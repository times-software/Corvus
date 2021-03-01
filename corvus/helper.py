from corvus.structures import Handler, Exchange, Loop, Update
import numpy as np
import corvutils.pyparsing as pp
import os, sys, subprocess, shutil, resource
import math
import pprint
import copy

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
        print('Inside helper')
        from corvus.controls import generateAndRunWorkflow
        dir = config['xcDir']
        for target in output:
            print('Inside helper', output)
            if (target == 'xanes_cfavg'):
                #Define the target of the average
                targetList = [['xanes']]

                cluster_array = input['cluster_array']
                for i,clust_elem in enumerate(cluster_array):
                    print(i, clust_elem[0])
                    input['absorbing_atom'] = [[clust_elem[0]]]
                    input['cluster'] = clust_elem[1]
                    # Make sure we are working with absolute units.
                    input['feff.absolute'] = [[True]]
                    config2 = copy.deepcopy(config)
                    config2['cwd'] = config['xcDir']
                    config2['xcIndexStart'] = i+1
                    targetList = [['xanes']]
                    generateAndRunWorkflow(config2, input, targetList)

                    # get results from input.
                    en,mu=np.array(input['xanes'])

                    # If first run, make the common grid.
                    if i == 0:
                        egrid = en
                        mu_total = mu
                    else:
                        # interpolate onto the common grid and add to total.
                        mu_total = mu_total + np.interp(egrid, en, mu, left = 0.0, right = 0.0)


                output['xanes_cfavg'] = [egrid.tolist(),(mu_total/len(cluster_array)).tolist()]

    @staticmethod
    def cleanup(config):
        pass


