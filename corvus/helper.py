from corvus.structures import Handler, Exchange, Loop, Update
import re
from time import sleep
import numpy as np
import random
from scipy.interpolate import RegularGridInterpolator as rgi
import corvutils.pyparsing as pp
import os, sys, subprocess, shutil #, resource
import math
import pprint
import copy
from corvus.abort import check_abort
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

implemented['loop'] = {'type':'Exchange','out':['loop'],'cost':1,'req':['loop_target'],
'desc':'Run multiple calculations with different settings.'}
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
       return
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
        if check_abort(None,'helper.run'): return
        #print('Inside helper')
        from corvus.controls import generateAndRunWorkflow
        #dir = config['xcDir']
        for target in output:
            #print('Inside helper', output)
            if (target == 'cfavg'):
                if "multiprocessing_ncpu" in input:
                    ncpu = input["multiprocessing_ncpu"][0][0]
                else:
                    ncpu = mltp.cpu_count()

                print("Number of threads available for multiprocessing: ", ncpu)
                #Define the target of the average
                targetList = input['cfavg_target']
                for targ in targetList[0]:
                    # Set the current working directory to the correct directory. 
                    if 'xcLabel' in config:
                        subdir = config['pathprefix'] + str(config['xcIndex']) + config['xcLabel'] + '_cfavg_' + targ
                    else:
                        subdir = config['pathprefix'] + str(config['xcIndex']) + '_cfavg_' + targ
                    config['xcIndex'] = config['xcIndex']+1

                    xcDir = os.path.join(config['cwd'], subdir)
                    # Make new output directory if it doesn't exist
                    if not os.path.exists(xcDir):
                        os.mkdir(xcDir)

                    cluster_array = input['cluster_array']
                    print("Number of absorbers:", len(cluster_array))
                    dirs=[]
                    
                    

                    # Set total number of processes
                    if 'cfavg_max_configurations' in input and len(cluster_array) > 1:
                        # Use nmax randomly chosen configurations
                        totprocs = min(input['cfavg_max_configurations'][0][0],len(cluster_array))
                        nclust=totprocs
                        if 'cfavg_choose_random_absorbers' in input: 
                             if input['cfavg_choose_random_absorbers'][0][0]:
                                absorbers=random.sample(range(1, len(cluster_array)), totprocs)
                                for iabs in absorbers:
                                   new_cluster_array = [cluster_array[i] for i in absorbers]

                                cluster_array = new_cluster_array
                    else:
                        totprocs = len(cluster_array)
                        nclust=totprocs
                    
     
                    outputs = []
                    numdone=0
                    while totprocs > 0:
                        poolSize = min(ncpu,totprocs)
                        print("Using ", poolSize, ' processors.')
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
                            configs[i]['cwd'] = xcDir
                            configs[i]['xcIndexStart'] = i+numdone+1
                            
                            # Create a safe folder name from input
                            invalid_chars = r'[<>:"/\|?*`]'
                            folder_name = re.sub(invalid_chars, '', clust_elem[3])
                            folder_name = re.sub(r'[\x00-\x1f\x7f]', '', folder_name)
                            configs[i]['xcLabel'] = folder_name

                            tLists = tLists + [[targ]]
                            arguments = arguments + [(configs[i],inputs[i],[[targ]])]
                            #targetList = [['xanes']]
                        with mltp.Pool(processes=poolSize) as pool:
                            poolout =  pool.starmap_async(multiproc_genAndRun,arguments)
                            while not poolout.ready():
                                sleep(1)
                        outputs = outputs + poolout.get()
                        numdone = numdone + poolSize
                       
                        #print('Check: ', len(outputs), poolSize, totprocs)
                        pool.close()
                        totprocs = totprocs - poolSize
                    
                    if input['write_input_only'][0][0]: 
                        output[target] = None
                        continue
                    #print(len(outputs),len(cluster_array))
                        #generateAndRunWorkflow(config2, input, targetList)
                        
                        #Prcs = mltp.Process(target=generateAndRunWorkflow,args=(config2,inputs[i],targetList))
                        #Prcs.start()
                        #tasks.append(Prcs)

                    #for Prcss in tasks:
                    #    Prcss.join()
                    mu_pol = []
                    data_pol = []
                    ipol = 1
                    npol =  len(outputs[0][targ])
                    #print(outputs[0][targetList[0][0]])
                    #print('npol=',npol)
                    #exit()
                    UnicodeEncodeError = []
                    print('Summing XAS of all unique absorbers:')
                    while ipol <= npol-1:
                        en = []
                        mu = []
                        step = 1.e30
                        totalWeight = 0.0
                        weights = []
                        ndim=0
                        
                        for i,clust_elem in enumerate(cluster_array[0:nclust]):
                            # get results from inputs.
                            #print(targetList[0][0])
                            #print(inputs[i])
                            weight = clust_elem[1]
                            if ipol == 1:
                               print('Absorber: ', clust_elem[3])
                               print('weight ratio: ', clust_elem[1])
                            weights = weights + [weight]
                            #mu0 = mu0
                            totalWeight = totalWeight + weight
                            data = np.array(outputs[i][targ])
                            ndim = data.ndim
                            if data.ndim == 2: # array of rows corresponding to set of 1d data 
                                en0 = data[0]
                                mu0 = data[ipol]
                                # Save in array of output.
                                en = en + [en0]
                                step = min(step,np.amin(en0[1:]-en0[:-1]))
                                mu = mu + [mu0]
                                #plt.plot(en,mu)
                            elif data.ndim == 3: # 2d data like RIXS. Assummes first two indices are x and y.
                                # Now reform data as a 2d ndarray
                                dataNd = dataNd + [data]
                            
                            
                        
                        weights = np.array(weights)/totalWeight
                        if data.ndim == 2:
                            en = np.array(en)
                            mu = np.array(mu)

                            # Make the common grid.
                            emin=np.amin(en)
                            emax=np.amax(en)
                        if ipol == 1:
                            step = max(step,0.01)
                            egrid = np.arange(emin,emax,step)

                        # Interpolate onto common grid.
                        mu_interp = []

                        for i,clust_elem in enumerate(cluster_array[0:nclust]):
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
                 

                        if ndim == 2:
                            # Get average and standard deviation.
                            mu_avg,mu_stdev = weighted_avg_and_std(np.array(mu_interp), weights)
                            mu_pol = mu_pol + [mu_avg]
                            
                        elif ndim == 3:
                            data_pol = data_pol + [data_tot]
                            break
                        #print("mu_pol",mu_pol)
                        #mu_avg,mu_stdev = weighted_avg_and_std(mu_interp)
                        #mu_stdev = np.std(mu_interp,0)/totalWeight*len(cluster_array)
                        ipol = ipol + 1
                    dirName = os.path.join(config['cwd'], config['pathprefix'] + '.cfavg_' + targ)
                    if data.ndim == 2: 
                        mu_pol = [egrid] + mu_pol 
                        output['cfavg'] = np.array(mu_pol).tolist()
                        #print(mu_pol[0])
                        #print(mu_pol[1])              
                        np.savetxt(dirName + '.out',np.array(mu_pol).T)
                        

                    elif data.ndim == 3:
                        # Transform data for output
                        data_out = np.array([dataNd[0][0].flatten(),dataNd[0][1].flatten(),data_tot])
                        output['cfavg'] = np.array(data_out).tolist()
                        
                        f = open(dirName + '.1.out', 'w')
                        i=0
                        nd = dataNd[0][0].shape[0]
                        #print(nd)
                        for row in data_out.T:
                            if i == nd:
                                f.write('\n')
                                i = 0

                            f.write('    '.join(map(str,row)) + '\n')
                            i += 1

            if (target == 'loop'):
                #Define the target of the average
                targetList = input['loop_target']

                # Get loop parameter and loop values
                # Parameter is on first line.
                loop_parameter = input['loop_parameter'][0][0]
                # values are on next N lines.
                loop_values = input['loop_parameter'][1:]
                
                print("Number of calculations in loop:", len(loop_values))
                dirs=[]
                
                

                # Set total number of processes
                totprocs = len(loop_values)
                ncalc=totprocs
                
 
                outputs = []
                numdone=0
                #while totprocs > 0:
                #    poolSize = min(ncpu,totprocs)
                #    print("Using ", poolSize, ' processors.')
                #    print("processes left to run: ", totprocs)
                inputs = []
                configs = []
                tLists = []
                arguments = []
                for i,value in enumerate(loop_values):
                    inputs = inputs + [copy.copy(input)]
                    #print(loop_values[i])
                    inputs[i][loop_parameter] = [loop_values[i]]
                    #inputs[i]['multiprocessing_ncpu'] = [[1]]
                    del inputs[i]['target_list']
                    del inputs[i]['loop_parameter']
                    inputs[i]['target_list'] = targetList

                    # Set the current working directory to the correct directory. 
                    if 'xcLabel' in config:
                        subdir = config['pathprefix'] + str(i+1) + config['xcLabel'] + '_loop'
                    else:
                        subdir = config['pathprefix'] + str(i+1) + '_loop'

                    xcDir = os.path.join(config['cwd'], subdir)
                    # Make new output directory if it doesn't exist
                    if not os.path.exists(xcDir):
                        os.mkdir(xcDir)

                    dirs = dirs + [xcDir]
                    configs = configs + [copy.copy(config)]
                        
                    configs[i]['xcIndexStart'] = 1
                    configs[i]['cwd'] = xcDir
                    #configs[i]['pathprefix'] = xcDir
                        
                    tLists = tLists + [targetList]
                    arguments = arguments + [(configs[i],inputs[i],targetList)]
                        #targetList = [['xanes']]
                    #with mltp.Pool(processes=poolSize) as pool:
                    #    poolout =  pool.starmap_async(multiproc_genAndRun,arguments)
                    #    while not poolout.ready():
                    #        sleep(1)
                    generateAndRunWorkflow(configs[i], inputs[i], targetList) 
                    outputs = outputs + [output]
                
                if input['write_input_only'][0][0]: 
                    output[target] = None
                    return

                output[target] = [dirs]



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
