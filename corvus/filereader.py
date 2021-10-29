from corvus.structures import Handler, Exchange, Loop, Update
import corvutils.pyparsing as pp
import os, sys, subprocess, shutil #, resource
import re
from scipy.interpolate import CubicSpline
from scipy.integrate import quad
import numpy as np
# Debug: FDV
import pprint

pp_debug = pprint.PrettyPrinter(indent=4)


# Define dictionary of implemented calculations
implemented = {}
strlistkey = lambda L:','.join(sorted(L))
subs = lambda L:[{L[j] for j in range(len(L)) if 1<<j&k} for k in range(1,1<<len(L))]
#for s in subs(['cell_vectors', 'cell_struct_xyz_red', 'cell_scaling_iso', 'cell_scaling_abc', 'number_density']):
#    key = strlistkey(s)
#    autodesc = 'Get ' + ', '.join(s) + ' using cif2cell'
#    cost = 10
#    implemented[key] = {'type':'Exchange','out':list(s),'req':['cif_input'],
#                        'desc':autodesc,'cost':cost}

implemented['xanes'] = {'type':'Exchange','out':['xanes'],'cost':0,
                        'req':['xanes_file'],'desc':'Read xanes from file.'}

implemented['spectralFunction'] = {'type':'Exchange','out':['spectralFunction'],'cost':0,
                        'req':['spectralFunction_file'],'desc':'Read spectral function from file.'}


class filereader(Handler):
    def __str__(self):
        return 'filereader Handler'

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
            unresolved = {o for o in output if not filereader.canProduce(o)}
            canProduce = (o for o in output if filereader.canProduce(o))
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
            raise LookupError('Corvus cannot currently produce ' + key + ' using FEFF')
        return implemented[key]['cost']

    @staticmethod
    def sequenceFor(output,inp=None):
        if isinstance(output, list) and output and isinstance(output[0], str):
            key = strlistkey(output)
        elif isinstance(output, str):
            key = output
        else:
            raise TypeError('Output should be token of list of tokens')
        if key not in implemented:
            raise LookupError('Corvus cannot currently produce ' + key + ' using FEFF')
        f = lambda subkey : implemented[key][subkey]
        if f('type') == 'Exchange':
            return Exchange(filereader, f('req'), f('out'), cost=f('cost'), desc=f('desc'))

    @staticmethod
    def prep(config):
        subdir = config['pathprefix'] + str(config['xcIndex']) + '_CIF2CELL'
        xcDir = os.path.join(config['cwd'], subdir)
        # Make new output directory if if doesn't exist
        if not os.path.exists(xcDir):
            os.mkdir(xcDir)
        # Store current Exchange directory in configuration
        config['xcDir'] = xcDir

    #@staticmethod
    #def setDefaults(input,target):

    @staticmethod
    def run(config, input, output):


          
        # Loop over targets in output.
        if 'xanes' in output:
            fl = input['xanes_file'][0][0]
            w, mu0 = np.loadtxt(fl,usecols = (0,1)).T
            output['xanes'] = [w, mu0]

        if 'spectralFunction' in output:
            fl = input['spectralfunction_file'][0][0]
            w, spfcn = np.loadtxt(fl,usecols = (0,1)).T
            output['spectralFunction'] = [w, spfcn]
            

    @staticmethod
    def cleanup(config):
        pass





