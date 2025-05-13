from corvus.structures import Handler, Exchange, Loop, Update
import numpy as np
import random
from scipy.interpolate import RegularGridInterpolator as rgi
import corvutils.pyparsing as pp
import os, sys, subprocess, shutil #, resource
import math
import pprint
import copy
#import mltp
#try:
#    import ray.util.multiprocessing as mltp
#except ImportError:
import multiprocessing as mltp

#from matplotlib import pyplot as plt

#Define dictionary of implemented calculations
implemented = {}
strlistkey = lambda L:','.join(sorted(L))

implemented['cfavg'] = {'type':'Exchange','out':['cfavg'],'cost':1,'req':['cluster_array'],
'desc':'Average over an array of clusters and absorbing atoms.'}

#implemented['spectrum_set'] = {'type':'Exchange','out':['scan'],'cost':1,'req':['parameter_scan'],
#'desc':'Calculate an array of spectra from an array of input parameters, absorbing atoms, and inputs.'}


class helper(Handler):
    def __str__(self):
        return 'helper Handler'

    @staticmethod
    def canProduce(output):
        if isinstance(output, list) and output and isinstance(output[0], str):
            return strlistkey(output) in implemented
        elif isinstance(output, str):
            return output in implemented
        else:
            raise TypeError('Output should be token or list of tokens')

    @staticmethod
    def requiredInputFor(output):
        if isinstance(output, list) and output and isinstance(output[0], str):
            unresolved = {o for o in output if not Feff.canProduce(o)}
            canProduce  = (o for o in output if Feff.canProduce(o))
            additionalInput = (set(implemented[o]['req']) for o in canProduce)
            return list(set.union(unresolved,*additionalInput))
        elif isinstance(output, str):
            if output in implemented:
                return implemented[output]['req']
            else:
                return [output]
        else:
            raise TypeError('Output should be token or list of tokens')

    @staticmethod
    def cost(output):
        if isinstance(output, list) and output and isinstance(output[0], str):
            key = strlistkey(output)
        elif isinstance(output, str):
            key = output
        else:
            raise TypeError('Output should be token or list of tokens')
        if key not in implemented:
            raise LookupError('Corvus cannot currently produce ' + key + ' using helper')
        return implemented[key]['cost']

    @staticmethod
    def prep(config):
       if 'xcLabel' in config:
           subdir = config['pathprefix'] + str(config['xcIndex']) + config['xcLabel'] + '_helper'
       else:
           subdir = config['pathprefix'] + str(config['xcIndex']) + '_helper'

       xcDir = os.path.join(config['cwd'], subdir)
       # Make new output directory if it doesn't exist
       if not os.path.exists(xcDir):
           os.mkdir(xcDir)
       # Store current Exchange directory in configuration
       config['xcDir'] = xcDir
    
    @staticmethod
    def sequenceFor(output,inp=None):
        from corvus.controls import availableHandlers
        
        if isinstance(output, list) and output and isinstance(output[0], str):
            key = strlistkey(output)
        elif isinstance(output, str):
            key = output
        else:
            raise TypeError('Output should be token or list of tokens')
        if key not in implemented:
            raise LookupError('Corvus cannot currently produce ' + key + ' using helper')
        f = lambda subkey : implemented[key][subkey]
        required = f('req')
        
        if f('type') == 'Exchange':
            return Exchange(helper, required, f('out'), cost=f('cost'), desc=f('desc'))

    @staticmethod
    def run(config, input, output):
        #print('Inside helper')
        from corvus.controls import generateAndRunWorkflow
        dir = config['xcDir']
        if "multiprocessing.ncpu" in input:
           ncpu = input["multiprocessing.ncpu"][0][0]
        else:
           ncpu = mltp.cpu_count()

        print("Number of cpu for multiprocessing: ", ncpu)
        for target in output:
            #print('Inside helper', output)
            if (target == 'cfavg'):
                #Define the target of the average
                targetList = input['cfavg.target']

                cluster_array = input['cluster_array']
                print("Number of absorbers:", len(cluster_array))
                en = []
                mu = []
                dataNd = []
                step = 1.e30
                totalWeight = 0.0
                weights = []
                dirs=[]
                
                

                # Set total number of processes
                if 'cfavg.max_configurations' in input and len(cluster_array) > 1:
                    # Use nmax randomly chosen configurations
                    totprocs = min(input['cfavg.max_configurations'][0][0],len(cluster_array))
                    if 'cfavg.choose_random_absorbers' in input: 
                       if input['cfavg.chose_random_absorbers'][0][0]:
                          absorbers=random.sample(range(1, len(cluster_array)), totprocs)
                          for iabs in absorbers:
                             new_cluster_array = [cluster_array[i] for i in absorbers]

                          cluster_array = new_cluster_array
                else:
                    totprocs = len(cluster_array)
                
 
                outputs = []
                numdone=0
                while totprocs > 0:
                    poolSize = min(ncpu,totprocs)
                    #print("Pool size: ", poolSize)
                    print("processes left to run: ", totprocs)
                    inputs = []
                    configs = []
                    tLists = []
                    arguments = []
                    for i,clust_elem in enumerate(cluster_array[numdone:numdone+poolSize]):
                        inputs = inputs + [copy.copy(input)]
                        inputs[i]['absorbing_atom'] = [[clust_elem[0]]]
                        inputs[i]['cluster'] = clust_elem[2]
                        # Make sure we are working with absolute units.
                        inputs[i]['feff.absolute'] = [[True]]
                        configs = configs + [copy.copy(config)]
                        configs[i]['cwd'] = config['xcDir']
                        configs[i]['xcIndexStart'] = i+numdone+1
                        configs[i]['xcLabel'] = clust_elem[3]
                        tLists = tLists + [targetList]
                        arguments = arguments + [(configs[i],inputs[i],targetList)]
                        #targetList = [['xanes']]
                    pool = mltp.Pool(processes=poolSize)
                    outputs = outputs + pool.starmap(multiproc_genAndRun,arguments)
                    numdone = numdone + poolSize
                   
                    #print('Check: ', len(outputs), poolSize, totprocs)
                    pool.close()
                    totprocs = totprocs - poolSize

                #print(len(outputs),len(cluster_array))
                    #generateAndRunWorkflow(config2, input, targetList)
                    
                    #Prcs = mltp.Process(target=generateAndRunWorkflow,args=(config2,inputs[i],targetList))
                    #Prcs.start()
                    #tasks.append(Prcs)

                #for Prcss in tasks:
                #    Prcss.join()

                for i,clust_elem in enumerate(cluster_array):
                    # get results from inputs.
                    #print(targetList[0][0])
                    #print(inputs[i])
                    weight = clust_elem[1]
                    weights = weights + [weight]
                    #mu0 = mu0
                    totalWeight = totalWeight + weight
                    data = np.array(outputs[i][targetList[0][0]])
                    if data.ndim == 2: # array of rows corresponding to set of 1d data 
                       en0,mu0=data
                       # Save in array of output.
                       en = en + [en0]
                       step = min(step,np.amin(en0[1:]-en0[:-1]))
                       mu = mu + [mu0]
                       #plt.plot(en,mu)
                    elif data.ndim == 3: # 2d data like RIXS. Assummes first two indices are x and y.
                       # Now reform data as a 2d ndarray
                       dataNd = dataNd + [data]
                       
                     
                   
                weights = np.array(weights)
                if data.ndim == 2:     
                   en = np.array(en)
                   mu = np.array(mu)

                   # Make the common grid.
                   emin=np.amin(en)
                   emax=np.amax(en)
                   egrid = np.arange(emin,emax,step)

                   # Interpolate onto common grid.
                   mu_interp = []

                for i,clust_elem in enumerate(cluster_array):
                    if data.ndim == 2:
                        # interpolate onto the common grid and add to total.
                       
                        mui = np.interp(egrid, en[i], mu[i], left = 0.0)
                        mu_interp = mu_interp + [mui]


                    elif data.ndim == 3:
                        # interpolate onto common 2d grid - just use the first grid
                        print('Adding contribution from absorber ', i)
                        datai = rgi((dataNd[i][1,:,0],dataNd[i][0,0,:]),dataNd[i][2],method='linear', bounds_error=False,fill_value=0.0)
                        if i == 0:
                            data_tot = datai((dataNd[0][1],dataNd[0][0])).flatten()*weights[i]
                        else:
                            data_tot = data_tot + datai((dataNd[0][1],dataNd[0][0])).flatten()*weights[i]

                if data.ndim == 2:                
                    # Get average and standard deviation.
                    weights = weights/totalWeight
                    mu_avg,mu_stdev = weighted_avg_and_std(mu_interp, weights)
                    np.savetxt(config['pathprefix'] + '.cfavg.' + targetList[0][0] + '.out',np.array([egrid,mu_avg,mu_stdev]).T)
                    output['cfavg'] = np.array([egrid,mu_avg,mu_stdev]).tolist()

                elif data.ndim == 3:
                    # Transform data for output
                    data_out = np.array([dataNd[0][0].flatten(),dataNd[0][1].flatten(),data_tot])
                    output['cfavg'] = np.array(data_out).tolist()                    
                    
                    f = open(config['pathprefix'] + '.cfavg.' + targetList[0][0] + '.out', 'w')
                    i=0
                    nd = dataNd[0][0].shape[0]
                    print(nd)
                    for row in data_out.T:
                        if i == nd:
                            f.write('\n')
                            i = 0

                        f.write('    '.join(map(str,row)) + '\n')
                        i += 1
            #elif(target == 'spectrum_set'):
                # Loop through set of parameters, create and run the
                # Set the target of spectrum_set - XANES, XES, RIXS, ...
                #targetList = input['spectrum_set.target']

                #parameter_set = input['parameter_set']
                #en = []
                #mu = []
                #step = 1.e30
                #totalWeight = 0.0
                #weights = []
                #dirs=[]




    @staticmethod
    def cleanup(config):
        pass

def multiproc_genAndRun(conf,inp,tList):
    from corvus.controls import generateAndRunWorkflow
    generateAndRunWorkflow(conf, inp, tList)
    return inp

def weighted_avg_and_std(values, weights):
    """
        Return the weighted average and standard deviation.
    
        values, weights -- Numpy ndarrays with the same shape.
        """
    average = np.average(values, axis=0, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, axis=0, weights=weights)
    return (average, np.sqrt(variance))
