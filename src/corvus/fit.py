from structures import Handler, Exchange, Loop, Update
import numpy as np
import corvutils.pyparsing as pp
import os, sys, subprocess, shutil, resource
from lmfit import Minimizer, Parameters, report_fit, fit_report 
import re
import math
# Debug: FDV
import pprint
pp_debug = pprint.PrettyPrinter(indent=4)

# Define dictionary of implemented calculations
implemented = {}
fitconvfile=None
firstcall=True
strlistkey = lambda L:','.join(sorted(L))
implemented['fit'] = {'type':'Exchange','out':['fit'],'cost':3,
                            'req':['fit.target','fit.datafile','fit.parameters'],'desc':'Calculate a property using fit handler.'}


class fit(Handler):
    def __str__(self):
        return 'fit Handler'

    @staticmethod
    def canProduce(output):
        if isinstance(output, list) and output and isinstance(output[0], basestring):
            return strlistkey(output) in implemented
        elif isinstance(output, basestring):
            return output in implemented
        else:
            raise TypeError('Output should be token or list of tokens')

    @staticmethod
    def requiredInputFor(output):
        if isinstance(output, list) and output and isinstance(output[0], basestring):
            unresolved = {o for o in output if not Feff.canProduce(o)}
            canProduce = (o for o in output if Feff.canProduce(o))
            additionalInput = (set(implemented[o]['req']) for o in canProduce)
            return list(set.union(unresolved,*additionalInput))
        elif isinstance(output, basestring):
            if output in implemented:
                return implemented[output]['req']
            else:
                return [output]
        else:
            raise TypeError('Output should be token or list of tokens')

    @staticmethod
    def cost(output):
        if isinstance(output, list) and output and isinstance(output[0], basestring):
            key = strlistkey(output)
        elif isinstance(output, basestring):
            key = output
        else:
            raise TypeError('Output should be token or list of tokens')
        if key not in implemented:
            raise LookupError('Corvus cannot currently produce ' + key + ' using fit')
        return implemented[key]['cost']

    @staticmethod
    def sequenceFor(output,inp=None):
        from controls import availableHandlers

        if isinstance(output, list) and output and isinstance(output[0], basestring):
            key = strlistkey(output)
        elif isinstance(output, basestring):
            key = output
        else:
            raise TypeError('Output should be token of list of tokens')
        if key not in implemented:
            raise LookupError('Corvus cannot currently produce ' + key + ' using fit')
        f = lambda subkey : implemented[key][subkey]
        required = f('req')
        # JJK - Need to add requirements of internal workflow here.
        if 'fit.target' in inp.keys():
            target = inp['fit.target'][0]
            for h in availableHandlers():
                if h.canProduce(target):
                   required.extend(h.requiredInputFor(target))
                   break
        else:
            print("Missing user input: fit.target is required to run a fit.")
            sys.exit()
        
        if f('type') is 'Exchange':
            return Exchange(fit, required, f('out'), cost=f('cost'), desc=f('desc'))

    @staticmethod
    def prep(config):
        subdir = config['pathprefix'] + str(config['xcIndex']) + '_fit'
        xcDir = os.path.join(config['cwd'], subdir)
        # Make new output directory if if doesn't exist
        if not os.path.exists(xcDir):
            os.mkdir(xcDir)
        # Store current Exchange directory in configuration
        config['xcDir'] = xcDir

    #@staticmethod
    #def setDefaults(input,target):

    # JJ Kas - run now performs all 3 methods, i.e., generateInput, run, translateOutput
    # Maybe we should also include prep here. Is there a reason that we want to limit the directory names
    # to automated Corvus_FEFFNN? Also if we have prep included here, we can decide on making a new directory
    # or not. 
    @staticmethod
    def run(config, input, output):
        dataFile = os.path.join(config['cwd'], input['fit.datafile'][0][0])
        x,data = np.loadtxt(dataFile,usecols=(0,1),unpack=True) 
        params = Parameters()
        # Loop over parameters defined in input. Later, we might loop over
        # all possible parameters and set which ones vary instead.
        for param in input['fit.parameters']:
            params.add(param[0],value=param[1])

        #open('fitconvergence.dat', 'ab')
        # do fit, here with the default leastsq algorithm
        minner = Minimizer(Xanes2Min, params, fcn_args=(x, data, input, config, output))
        result = minner.minimize(epsfcn=0.0001)
        final  = data + result.residual

        report_fit(result)
        with open('fit_result.txt', 'w') as fh:
            fh.write(fit_report(result))

        output['fit'] = [x.tolist(),final.tolist(), data.tolist()] # For now set this as fit. Later we may want to make a statement about what is implemented and what is not for fit. 

    @staticmethod
    def cleanup(config):
        pass


# Global functions
def Xanes2Min(params, x, data, input, config, output):
    from controls import generateAndRunWorkflow
    import copy
    Xanes2Min.count
    #lastParams
    Xanes2Min.count += 1
    input2 = copy.deepcopy(input)

    energy_shift = 0.0
    atoms = input['cluster']
    amp = 1.0
    input2['feff.corrections'] = [[0.0,0.0]]
    # Set controls based on what has changed since last call
    # to function.
    if Xanes2Min.count == 1:
        control = [1,1,1,1,1,1]
    else:
        control = [0, 0, 0, 0, 0, 0]
    ipar = 0
    for param in params.values():
        if Xanes2Min.lastParams is not None:            
            diff = param != Xanes2Min.lastParams.values()[ipar] 
        else: 
            diff = True 

        # Use case insensitive equal.
        if param.name.lower() == 'expansion':
            # Uniform expansion of coordinates in cluster.
            if diff:
                control = [1, 1, 1, 1, 1, 1] 

            expansion = param.value
        
            atoms = [[f[0],expansion*f[1],expansion*f[2],expansion*f[3]] for f in atoms]

        elif param.name.lower() == 'broadening':
            # Lorentzian broadening applied to spectrum.
            if diff:
                control[5] = 1

            broadening = param.value
            #input2['spectral_broadening'] = [[broadening]]
            input2['feff.corrections'][0][1] = broadening

        elif param.name.lower() == 'delta_e0':
            # Shift in absolute edge energy (shift of energy grid of spectrum).
            energy_shift = param.value
        elif param.name.lower() == 'bond':
            # Move a set of atoms away from absorber along a bond.
            # Find vector to move along (r_2 - r_1)/r12
            # Get the two atoms defining the bond vector.
            if diff:
                control = [1, 1, 1, 1, 1, 1] 

            bond = param.value
            bond_atoms = [item-1 for sublist in input2['fit.bond'] for item in sublist]
            vec = [ input2['cluster'][bond_atoms[1]][i] - input2['cluster'][bond_atoms[0]][i] for i in [1,2,3]]
            vecSquared = [ vec[i]**2 for i in [0,1,2] ] 
            norm = math.sqrt(sum(vecSquared))
            vec = [ vec[i]/norm*bond for i in [0,1,2] ]
            for atom in bond_atoms[1:]:
                for i in [1,2,3]:
                    atoms[atom][i] += vec[i-1]

        elif param.name.lower() == 'delta_efermi':
            #input2['fermi_shift'] = [[param.value]]
            input2['feff.corrections'][0][0] = param.value
            if diff:
                control[5] = 1
        elif param.name.lower() == 'amplitude':
            amp = param.value
        else:
            print('WARNING: UNKOWN PARAMETER ' + param.name + '!')
            print('STOPPING NOW!!!')
            exit()
        
        ipar += 1

    input2['cluster'] = atoms
    input2['feff.control'] = [control]
    # Need a copy of config to start wf over
    config2 = copy.deepcopy(config)
    
    # Set current working directory to xCDir, so that internal wf
    # will run inside of outer wf directory.
    config2['cwd'] = config['xcDir']

    if False: # Save all runs of underlying handler in separate directories.
        config2['xcIndexStart'] = Xanes2Min.count
    else:
        config2['xcIndexStart'] = 1

    dir = config['xcDir']
    # Loop over targets in output. Not sure if there will ever be more than one output target here.
    # Set output and error files
    for target in output:
        with open(os.path.join(dir, 'corvus.fit.stdout'), 'w') as out, open(os.path.join(dir, 'corvus.fit.stderr'), 'w') as err:

            # Set the tagetList according to fit.target
            #xanes2MinIterator += 1
            targetList = input['fit.target']

            # generate and run the workflow for target, unless no run is necessary.
            generateAndRunWorkflow(config2, input2,targetList)
    
            x0,y=np.array(input2[input['fit.target'][0][0]])
            y = y*amp
            # If there is an energy shift, shift the x-axis before
            # interpolating onto data grid
            x0 = x0 + energy_shift
            # On first call, check if experimental data is outside calculated
            # data, and redefine experimental data within range.
            global firstcall
            global fitconvfile
            if firstcall:
               print('Opening convergence file')
               try:
                   os.remove('fitconvergence.dat')
               except OSError:
                   pass
               fitconvfile=open('fitconvergence.dat', 'a')
               np.savetxt(fitconvfile,np.array([x,data]).transpose())
               fitconvfile.write('\n')   
               firstcall = False

            yterp = np.interp(x, x0, y, left = 0.0, right = 0.0)
             
            np.savetxt(fitconvfile,np.array([x,yterp]).transpose())
            fitconvfile.write('\n')

            i = 0
            residual = np.zeros(yterp.size)
            for yi in yterp:
               if (x[i] >= x0[0]) and (x[i] <= x0[-1]):
                 residual[i] = yi - data[i] 
               else:
                 residual[i] = 0.0
               i = i + 1

            Xanes2Min.lastParams = copy.copy(params)

            return residual

Xanes2Min.count=0
Xanes2Min.lastParams=None
