from corvus.structures import Handler, Exchange, Loop, Update
import numpy as np
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
                
                dirs=[]
                
                

                # Set total number of processes
                totprocs = len(cluster_array)
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
                        configs[i]['cwd'] = config['xcDir']
                        configs[i]['xcIndexStart'] = i+numdone+1
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
                mu_pol = []
                ipol = 1
                UnicodeEncodeError = []
                while ipol <= 4:
                    en = []
                    mu = []
                    step = 1.e30
                    totalWeight = 0.0
                    weights = []
                    for i,clust_elem in enumerate(cluster_array):
                        # get results from inputs.
                        #print(targetList[0][0])
                        #print(inputs[i])
                        xns=np.array(outputs[i][targetList[0][0]])
                        en0=xns[0]
                        mu0=xns[ipol]
                        if False:
                           weight = clust_elem[1]
                        else:
                           weight = 1.0 # Make weight 1 (if polarization is requested).
                        
                        weights = weights + [weight]
                        #mu0 = mu0
                        totalWeight = totalWeight + weight
                        
                        # Save in array of XANES output.
                        en = en + [en0]
                        step = min(step,np.amin(en0[1:]-en0[:-1]))
                        mu = mu + [mu0]
                        #plt.plot(en,mu)

                    en = np.array(en)
                    mu = np.array(mu)
                    weights = np.array(weights)
                    # Make the common grid.
                    emin=np.amin(en)
                    emax=np.amax(en)
                    step = max(step,0.01)
                    if ipol == 1:
                        egrid = np.arange(emin,emax,step)
                    
                    # Interpolate onto common grid.
                    mu_interp = []
                    for i,clust_elem in enumerate(cluster_array):
                        # interpolate onto the common grid and add to total.
                        mui = np.interp(egrid, en[i], mu[i], left = 0.0)
                        mu_interp = mu_interp + [mui]

                    # Get average and standard deviation.
                    mu_avg,mu_stdev = weighted_avg_and_std(np.array(mu_interp), weights)
                    mu_pol = mu_pol + [mu_avg]
                    #print("mu_pol",mu_pol)
                    #mu_avg,mu_stdev = weighted_avg_and_std(mu_interp)
                    #mu_stdev = np.std(mu_interp,0)/totalWeight*len(cluster_array)
                    ipol = ipol + 1
                
                mu_pol = [egrid] + mu_pol
                
                output['cfavg'] = np.array(mu_pol).tolist()

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
