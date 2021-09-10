from corvus.structures import Handler, Exchange, Loop, Update
import corvutils.pyparsing as pp
import os, sys, subprocess, shutil #, resource
import re
from scipy.interpolate import CubicSpline
from scipy.integrate import quad
from scipy.signal import convolve
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

implemented['mbxanes'] = {'type':'Exchange','out':['mbxanes'],'cost':0,
                        'req':['xanes_cfavg','spectralFunction'],'desc':'Calculate many-body xanes from xanes and spectral function.'}
#'req':['xanes','spectal_function'],'desc':'Calculate supercell from cif input.'}



class mbconv(Handler):
    def __str__(self):
        return 'mbconv Handler'

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
            unresolved = {o for o in output if not mbconv.canProduce(o)}
            canProduce = (o for o in output if mbconv.canProduce(o))
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
        required = f('req')
        # JJK - Need to add requirements of internal workflow here.
        if 'mbconv' in list(inp.keys()):
            required.extend()

        if f('type') is 'Exchange':
            return Exchange(mbconv, f('req'), f('out'), cost=f('cost'), desc=f('desc'))

    @staticmethod
    def prep(config):
        subdir = config['pathprefix'] + str(config['xcIndex']) + '_MBXANES'
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
        if 'mbxanes' in output:
            # In future use file_reader handler to read in XANES and spectral function if already calculated.
            w  = np.array(input.get('xanes_cfavg')[0])
            mu0= np.array(input.get('xanes_cfavg')[1])
            wsf= np.flip(-1.0*np.array(input.get('spectralFunction')[0]))
            sf = np.flip(np.array(input.get('spectralFunction')[1]))
            # Interpolate both XANES and spectral function onto an even grid
            #w, mu0 = np.loadtxt('xanes.dat',usecols = (0,1)).T
            #wsf,sf = np.loadtxt('spfcn.dat',usecols = (0,1)).T
            min_diff = np.amin(np.ediff1d(w))
            min_diff = min(min_diff,np.amin(np.ediff1d(wsf)))
        
            mu0_cs = CubicSpline(w,mu0)
            spfcn_cs = CubicSpline(wsf,sf)
            # Use larger of two ranges to specify range
            w_terp = np.arange(w[0],w[-1],min_diff)
            wsf_terp = np.arange(wsf[0],wsf[-1],min_diff)
            mu0_terp = mu0_cs(w_terp)
            spfcn_terp = spfcn_cs(wsf_terp)
 
            mu_mb = convolve(mu0_terp,spfcn_terp,mode='full')*min_diff

            # If extra broadening is requested, perform a convolution of that as well.
            if 'mbconv.extra_broadening' in input:
                gam = input['mbconv.extra_broadening'][0][0]
                A_br = gam/np.pi*1.0/(wsf_terp**2 + gam**2)
                mu_mb = np.convolve(mu_mb,A_br,mode='same')*min_diff
              
            scale=w_terp[-1] - w_terp[0] + wsf_terp[-1] - wsf_terp[0]
            first = w_terp[0] + wsf_terp[0]
            w_terp = np.linspace(0.0,scale,mu_mb.size) 
            w_terp = w_terp + first
            mu0_terp = mu0_cs(w_terp)
            output['mbxanes'] = [w_terp,mu_mb]
            np.savetxt('mbxanes.dat',np.array([w_terp, mu_mb, mu0_terp]).transpose())



    @staticmethod
    def cleanup(config):
        pass





