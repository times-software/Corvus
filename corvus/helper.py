from corvus.structures import Handler, Exchange, Loop, Update
import numpy as np
import corvutils.pyparsing as pp
import os, sys, subprocess, shutil #, resource
import math
import pprint
import copy
#from matplotlib import pyplot as plt

#Define dictionary of implemented calculations
implemented = {}
strlistkey = lambda L:','.join(sorted(L))

implemented['xanes_cfavg'] = {'type':'Exchange','out':['xanes_cfavg'],'cost':1,'req':['cluster_array'],
'desc':'Average over an array of clusters and absorbing atoms.'}


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
        
        if f('type') is 'Exchange':
            return Exchange(helper, required, f('out'), cost=f('cost'), desc=f('desc'))

    @staticmethod
    def run(config, input, output):
        #print('Inside helper')
        from corvus.controls import generateAndRunWorkflow
        dir = config['xcDir']
        for target in output:
            #print('Inside helper', output)
            if (target == 'xanes_cfavg'):
                #Define the target of the average
                targetList = [['xanes']]

                cluster_array = input['cluster_array']
                print("Number of absorbers:", len(cluster_array))
                en = []
                mu = []
                step = 1.e30
                totalWeight = 0.0
                for i,clust_elem in enumerate(cluster_array):
                    #print(i, clust_elem[0])
                    input['absorbing_atom'] = [[clust_elem[0]]]
                    weight = clust_elem[1]
                    print('weight=',weight)
                    input['cluster'] = clust_elem[2]
                    # Make sure we are working with absolute units.
                    input['feff.absolute'] = [[True]]
                    config2 = copy.deepcopy(config)
                    config2['cwd'] = config['xcDir']
                    config2['xcIndexStart'] = i+1
                    targetList = [['xanes']]
                    generateAndRunWorkflow(config2, input, targetList)
            
                    # get results from input.
                    en0,mu0=np.array(input['xanes'])
                    mu0 = mu0*weight
                    totalWeight = totalWeight + weight
                    
                    # Save in array of XANES output.
                    en = en + [en0]
                    step = min(step,np.amin(en0[1:]-en0[:-1]))
                    mu = mu + [mu0]
                    #plt.plot(en,mu)

                en = np.array(en)
                mu = np.array(mu)
                # Make the common grid.
                emin=np.amin(en)
                emax=np.amax(en)
                egrid = np.arange(emin,emax,step)

                # Interpolate onto common grid.
                mu_interp = []
                for i,clust_elem in enumerate(cluster_array):
                    # interpolate onto the common grid and add to total.
                    mui = np.interp(egrid, en[i], mu[i], left = 0.0)
                    mu_interp = mu_interp + [mui]

                # Get average and standard deviation.
                mu_avg = np.average(mu_interp,0)/totalWeight*len(cluster_array)
                mu_stdev = np.std(mu_interp,0)/totalWeight*len(cluster_array)

                output['xanes_cfavg'] = np.array([egrid,mu_avg,mu_stdev]).tolist()

    @staticmethod
    def cleanup(config):
        pass


