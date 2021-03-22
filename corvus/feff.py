from corvus.structures import Handler, Exchange, Loop, Update
import corvutils.pyparsing as pp
import os, sys, subprocess, shutil #, resource
import re
import math
import numpy as np
import scipy.integrate
import glob
from scipy.signal import argrelextrema
#from CifFile import ReadCif
#from cif2cell.uctools import *

# Added by FDV
import time
#import threading as thrd
import multiprocessing as mltp
from pymatgen.io.cif import CifParser
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import more_itertools as mit

# Debug: FDV
import pprint
pp_debug = pprint.PrettyPrinter(indent=4)

# Do a linear interpolation of the YY data for point ii in XX
def Linear_Interp(ii,XX,YY):

  x_ii = XX[ii]
  if   ii==0:
    x0 = XX[1]
    x1 = XX[2]
    y0 = YY[1]
    y1 = YY[2]
  elif ii==len(XX):
    x0 = XX[-2]
    x1 = XX[-3]
    y0 = YY[-2]
    y1 = YY[-3]
  else:
    x0 = XX[ii-1]
    x1 = XX[ii+1]
    y0 = YY[ii-1]
    y1 = YY[ii+1]

  y_ii = y0 + (x_ii-x0)*(y1-y0)/(x1-x0)

  return y_ii

# Do a linear interpolation of the YY data for point ii in XX, with information
# about the points to use for reference
def Linear_Interp_2(ii,iim,iip,XX,YY):

  x_ii = XX[ii]
  x0   = XX[iim]
  x1   = XX[iip]
  y0   = YY[iim]
  y1   = YY[iip]

  y_ii = y0 + (x_ii-x0)*(y1-y0)/(x1-x0)

  return y_ii

# Small routine to temporarily fix the NaN issues
# A simpler routine assume the NaNs were isolated, but they can appear in pairs
# This new implementation takes care of that too, hpefully.
def Temp_Fix_NaN(k2,exafs,mu0):
  nEne = len(k2)
# Make a list of all points that have NaNs
  iNaNs = [ ik for ik,(abs1,abs2) in enumerate(zip(exafs,mu0)) if math.isnan(abs1) or math.isnan(abs1) ]
  if len(iNaNs) >0:
# Make a fake iNaNs to testing
# iNaNs = [ 0,1,2, 8,9, 20,21, nEne-2,nEne-1 ]
# Debug
    print('FDV iNaNs:',iNaNs)
# Group them into consecutive sets to simplify fixing
    iNaNs_Groups = [ list(grp) for grp in mit.consecutive_groups(iNaNs)]
# Debug
    print('FDV iNaNs_Groups:',iNaNs_Groups)
# Make sure that the groups of NaNs are not too big, since that would indicate
# bigger problems
    mxGroup = max([ len(grp) for grp in iNaNs_Groups ])
    if mxGroup > 3:
      print('Found long group of NaNs while trying to fix.')
      print('There is probably a bigger problem, stopping')
      sys.exit()
# Now we make a list of what points to use for the linear interpolation
    iNaNs_Intp = []
    for Grp in iNaNs_Groups:
      ilw = min(Grp)-1
      ihg = max(Grp)+1
      if ilw < 0:
        ilw = ihg
        ihg = ihg+1
      if ihg > nEne-1:
        ihg = ilw
        ilw = ilw-1
      for Ind in Grp:
        iNaNs_Intp.append([Ind,[ilw,ihg]])
# Debug
    print('FDV iNaNs_Intp:',iNaNs_Intp)
# Old code, saving for reference
# exafs_new = []
# mu0_new = []
# for ik,(abs1,abs2) in enumerate(zip(exafs,mu0)):
#   if math.isnan(abs1) or math.isnan(abs1):
#     print('Found NaN fixing:',ik)
#     abs1_new = Linear_Interp(ik,k2,exafs)
#     abs2_new = Linear_Interp(ik,k2,mu0)
#     exafs_new.append(abs1_new)
#     mu0_new.append(abs2_new)
#   else:
#     exafs_new.append(abs1)
#     mu0_new.append(abs2)
# Now we can interpolate in a smarter way, we only fix points that need fixing
    for Intp in iNaNs_Intp:
      ik = Intp[0]
      ilw = Intp[1][0]
      ihg = Intp[1][1]
      exafs[ik] = Linear_Interp_2(ik,ilw,ihg,k2,exafs)
      mu0[ik]   = Linear_Interp_2(ik,ilw,ihg,k2,mu0)

  return (exafs, mu0)

# Added by FDV
# List of feff modules that should be run in parallel if requested
Parallel_Exes = [ 'feff_timer', 'ldos', 'fms', 'pot' ]

# Define dictionary of implemented calculations
implemented = {}
strlistkey = lambda L:','.join(sorted(L))
subs = lambda L:[{L[j] for j in range(len(L)) if 1<<j&k} for k in range(1,1<<len(L))]
for s in subs(['potlist','atomlist']):
    key = strlistkey(s)
    autodesc = 'Basic FEFF ' + ', '.join(s) + ' from ABINIT cell definition'
    input = ['acell','znucl','xred','rprim','natom']
    cost = 1
    implemented[key] = {'type':'Exchange','out':list(s),'req':input,
                        'desc':autodesc,'cost':cost}

implemented['feffAtomicData'] = {'type':'Exchange','out':['feffAtomicData'],'cost':1,
                        'req':['cluster','absorbing_atom'],'desc':'Calculate atomic data using FEFF.'}

implemented['feffSCFPotentials'] = {'type':'Exchange','out':['feffSCFPotentials'],'cost':1,
                        'req':['cluster','absorbing_atom','feffAtomicData'],'desc':'Calculate SCF potentials using FEFF.'}

implemented['feffCrossSectionsAndPhases'] = {'type':'Exchange','out':['feffCrossSectionsAndPhases'],'cost':1,
                        'req':['cluster','absorbing_atom','feffSCFPotentials'],'desc':'Calculate atomic cross sections and phases using FEFF.'}

implemented['feffGreensFunction'] = {'type':'Exchange','out':['feffGreensFunction'],'cost':1,
                        'req':['cluster','absorbing_atom','feffCrossSectionsAndPhases'],'desc':'Calculate Greens function using FEFF.'}

implemented['feffPaths'] = {'type':'Exchange','out':['feffPaths'],'cost':1,
                        'req':['cluster','absorbing_atom','feffGreensFunction'],'desc':'Calculate paths using FEFF.'}

implemented['feffFMatrices'] = {'type':'Exchange','out':['feffFMatrices'],'cost':1,
                        'req':['cluster','absorbing_atom','feffPaths'],'desc':'Calculate scattering matrices using FEFF.'}

implemented['feffXANES'] = {'type':'Exchange','out':['feffXANES'],'cost':1,
                        'req':['cluster','absorbing_atom'],'desc':'Calculate XANES using FEFF.'}


implemented['feffXES'] = {'type':'Exchange','out':['feffXES'],'cost':1,
                        'req':['cluster','absorbing_atom'],'desc':'Calculate XANES using FEFF.'}

implemented['feffRIXS'] = {'type':'Exchange','out':['feffRIXS'],'cost':1,
                        'req':['cluster','absorbing_atom'],'desc':'Calculate XANES using FEFF.'}

implemented['opcons'] = {'type':'Exchange','out':['opcons'],'cost':1,
                        'req':['cif_input'],'desc':'Calculate optical constants using FEFF.'}

# Added by FDV
# Trying to implement and EXAFS with optimized geometry and ab initio DW factors
implemented['opt_dynmat_s2_exafs'] = {'type':'Exchange',
   'out':['opt_dynmat_s2_exafs'], 'cost':3,
   'req':['opt_dynmat','absorbing_atom'],
   'desc':'Calculate EXAFS with optimized geometry and ab initio DW factors from a dynamical matrix using FEFF.'}

class Feff(Handler):
    def __str__(self):
        return 'FEFF Handler'

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
            canProduce = (o for o in output if Feff.canProduce(o))
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
        if f('type') is 'Exchange':
            return Exchange(Feff, f('req'), f('out'), cost=f('cost'), desc=f('desc'))

    @staticmethod
    def prep(config):
        if 'xcIndexStart' in config:
            if config['xcIndexStart'] > 0:
                subdir = config['pathprefix'] + str(config['xcIndex']) + '_FEFF'
                xcDir = os.path.join(config['cwd'], subdir)
            else:
                xcDir = config['xcDir']
        else:
            subdir = config['pathprefix'] + str(config['xcIndex']) + '_FEFF'
            xcDir = os.path.join(config['cwd'], subdir)
            
        # Make new output directory if if doesn't exist
        if not os.path.exists(xcDir):
            os.makedirs(xcDir)
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
        # set atoms and potentials

        # Set directory to feff executables.
# Debug: FDV
#       pp_debug.pprint(config)
        feffdir = config['feff']
# Debug: FDV
#       sys.exit()
        
        # Copy feff related input to feffinput here. Later we will be overriding some settings,
        # so we want to keep the original input intact.
        feffInput = {key:input[key] for key in input if key.startswith('feff.')}

        # Generate any data that is needed from generic input and populate feffInput with
        # global data (needed for all feff runs.)
        if 'cif_input' in input: # Prefer using cif for input, but still use REAL space
            # Replace path with absolute path
            feffInput['feff.cif'] = [[os.path.abspath(input['cif_input'][0][0])]]
                
            if 'feff.reciprocal' not in input:
                feffInput['feff.real'] = [[True]]

        if 'cluster' in input:
            atoms = getFeffAtomsFromCluster(input)
            setInput(feffInput,'feff.atoms',atoms)
            potentials = getFeffPotentialsFromCluster(input)
            setInput(feffInput,'feff.potentials',potentials)
            
        debyeOpts = getFeffDebyeOptions(input)
        
        if 'feff.exchange' in feffInput:
            exch = feffInput['feff.exchange']
        else:
            exch = [[0, 0.0, 0.0, 2]] 

        if 'spectral_broadening' in input:
            exch[0][2] = input['spectral_broadening'][0][0]
        
        if 'fermi_shift' in input:
            exch[0][1] = input['fermi_shift'][0][0]

        feffInput['feff.exchange'] = exch
        if debyeOpts is not None:
            setInput(feffInput,'feff.debye',debyeOpts)

        # Set directory for this exchange
        exdir = config['xcDir']

        # Set input file
        inpf = os.path.join(exdir, 'feff.inp')
        # Loop over targets in output. Not sure if there will ever be more than one output target here.

        for target in output:
            if (target == 'feffAtomicData'):
                # Set output and error files
                with open(os.path.join(exdir, 'corvus.FEFF.stdout'), 'w') as out, open(os.path.join(exdir, 'corvus.FEFF.stderr'), 'w') as err:

                    # Write input file for FEFF.
                    writeAtomicInput(feffInput,inpf)

                    # Loop over executable: This is specific to feff. Other codes
                    # will more likely have only one executable.
                    # Run rdinp and atomic part of calculation
                    execs = ['rdinp','atomic','screen']
                    for exe in execs:
                        if 'feff.MPI.CMD' in feffInput:
                            executable = feffInput.get('feff.MPI.CMD')[0]
                            args = feffInput.get('feff.MPI.ARGS',[['']])[0] + [os.path.join(feffdir,exe)]
                        else:
                            executable = [os.path.join(feffdir,exe)]
                            args = ['']

                        runExecutable('',exdir,executable,args,out,err)

                # For this case, I am only passing the directory for now so
                # that other executables in FEFF can use the atomic data.

                output[target] = exdir

            elif (target == 'feffSCFPotentials'):
                # Set output and error files
                with open(os.path.join(exdir, 'corvus.FEFF.stdout'), 'w') as out, open(os.path.join(exdir, 'corvus.FEFF.stderr'), 'w') as err:

                    # Write input file for FEFF.
                    writeSCFInput(feffInput,inpf)

                    # Loop over executable: This is specific to feff. Other codes
                    # will more likely have only one executable. Here I am running 
                    # rdinp again since writeSCFInput may have different cards than

                    # Run rdinp and atomic part of calculation
                    execs = ['rdinp','atomic', 'pot', 'screen']
                    for exe in execs:
                        if 'feff.MPI.CMD' in feffInput:
                            executable = feffInput.get('feff.MPI.CMD')[0]
                            args = feffInput.get('feff.MPI.ARGS',[['']])[0] + [os.path.join(feffdir,exe)]
                        else:
                            executable = [os.path.join(feffdir,exe)]
                            args = ['']

                        runExecutable('',exdir,executable,args,out,err)

                # For this case, I am only passing the directory for now so
                # that other executables in FEFF can use the atomic data.

                output[target] = exdir
            elif (target == 'feffCrossSectionsAndPhases'):
                # Set output and error files
                with open(os.path.join(exdir, 'corvus.FEFF.stdout'), 'w') as out, open(os.path.join(exdir, 'corvus.FEFF.stderr'), 'w') as err:

                    # Write input file for FEFF.
                    writeCrossSectionsInput(feffInput,inpf)

                    # Loop over executable: This is specific to feff. Other codes
                    # will more likely have only one executable. Here I am running 
                    # rdinp again since writeSCFInput may have different cards than

                    # Run rdinp and atomic part of calculation
                    execs = ['rdinp','atomic','screen', 'pot', 'xsph']
                    for exe in execs:
                        if 'feff.MPI.CMD' in feffInput:
                            executable = feffInput.get('feff.MPI.CMD')[0]
                            args = feffInput.get('feff.MPI.ARGS',[['']])[0] + [os.path.join(feffdir,exe)]
                        else:
                            executable = [os.path.join(feffdir,exe)]
                            args = ['']

                        runExecutable('',exdir,executable,args,out,err)
                output[target] = exdir

 
            elif (target == 'feffGreensFunction'):

                # Set output and error files
                with open(os.path.join(exdir, 'corvus.FEFF.stdout'), 'w') as out, open(os.path.join(exdir, 'corvus.FEFF.stderr'), 'w') as err:


                    # Write input file for FEFF.
                    writeGreensFunctionInput(feffInput,inpf)

                    # Loop over executable: This is specific to feff. Other codes
                    # will more likely have only one executable. Here I am running 
                    # rdinp again since writeSCFInput may have different cards than

                    execs = ['rdinp','atomic','pot','screen','xsph','fms','mkgtr']
                    for exe in execs:
                        if 'feff.MPI.CMD' in feffInput:
                            executable = feffInput.get('feff.MPI.CMD')[0]
                            args = feffInput.get('feff.MPI.ARGS',[['']])[0] + [os.path.join(feffdir,exe)]
                        else:
                            executable = [os.path.join(feffdir,exe)]
                            args = ['']

                        runExecutable('',exdir,executable,args,out,err)

                # For this case, I am only passing the directory for now so
                # that other executables in FEFF can use the atomic data.

                output[target] = exdir


            elif (target == 'feffPaths'):
                # Set output and error files
                with open(os.path.join(exdir, 'corvus.FEFF.stdout'), 'w') as out, open(os.path.join(exdir, 'corvus.FEFF.stderr'), 'w') as err:

                    # Write input file for FEFF.
                    writePathsInput(feffInput,inpf)

                    # Loop over executable: This is specific to feff. Other codes
                    # will more likely have only one executable. Here I am running 
                    # rdinp again since writeSCFInput may have different cards than

                    execs = ['rdinp','atomic','pot','screen','xsph','fms','mkgtr','path']
                    for exe in execs:
                        if 'feff.MPI.CMD' in feffInput:
                            executable = feffInput.get('feff.MPI.CMD')[0]
                            args = feffInput.get('feff.MPI.ARGS',[['']])[0] + [os.path.join(feffdir,exe)]
                        else:
                            executable = [os.path.join(feffdir,exe)]
                            args = ['']

                        runExecutable('',exdir,executable,args,out,err)

                # For this case, I am only passing the directory for now so
                # that other executables in FEFF can use the atomic data.

                output[target] = exdir

            elif (target == 'feffFMatrices'):
                # Set output and error files
                with open(os.path.join(exdir, 'corvus.FEFF.stdout'), 'w') as out, open(os.path.join(exdir, 'corvus.FEFF.stderr'), 'w') as err:

                    # Write input file for FEFF.
                    writeFMatricesInput(feffInput,inpf)

                    # Loop over executable: This is specific to feff. Other codes
                    # will more likely have only one executable. Here I am running 
                    # rdinp again since writeSCFInput may have different cards than

                    execs = ['rdinp','atomic','pot','screen','xsph','fms','mkgtr','path','genfmt']
                    for exe in execs:
                        if 'feff.MPI.CMD' in feffInput:
                            executable = feffInput.get('feff.MPI.CMD')[0]
                            args = feffInput.get('feff.MPI.ARGS',[['']])[0] + [os.path.join(feffdir,exe)]
                        else:
                            executable = [os.path.join(feffdir,exe)]
                            args = ['']

                        runExecutable('',exdir,executable,args,out,err)
                        
                output[target] = exdir

            elif (target == 'feffXANES'):
                # Set output and error files
                with open(os.path.join(exdir, 'corvus.FEFF.stdout'), 'w') as out, open(os.path.join(exdir, 'corvus.FEFF.stderr'), 'w') as err:

                    # Write input file for FEFF.
                    writeXANESInput(feffInput,inpf)

                    # Loop over executable: This is specific to feff. Other codes
                    # will more likely have only one executable. Here I am running 
                    # rdinp again since writeSCFInput may have different cards than

#                   execs = ['rdinp','atomic','pot','screen','opconsat','xsph','fms','mkgtr','path','genfmt','ff2x','sfconv']
                    execs = ['feff_timer','rdinp','atomic','pot','screen','ldos','opconsat','xsph','fms','mkgtr','path','genfmt','ff2x','sfconv','feff_timer']
                    for exe in execs:
                        if 'feff.MPI.CMD' in feffInput:
# Modified by FDV
                            executable = feffInput.get('feff.MPI.CMD')[0]
                            npflag     = feffInput.get('feff.MPI.NPFLAG')[0]
                            nnflag     = feffInput.get('feff.MPI.NNFLAG')[0]
# Check the input we have here
#                           print('XXX feffInput XXX')
#                           pp_debug.pprint(feffInput)
                            nproc      = feffInput.get('feff.MPI.NP')[0]
                            print('feff.MPI.PPN',feffInput.get('feff.MPI.PPN')[0][0])
                            PPN_Val = feffInput.get('feff.MPI.PPN')[0][0]
                            nnodes = max(int(nproc[0]/PPN_Val),1)
                            print('nnodes',nnodes)
#                           sys.exit()
# Adjust nproc to get best efficiency, only run certain modules in parallel
                            if exe not in Parallel_Exes:
# Testing for bug
                              nproc = 1
#                             nproc = nproc[0]
                            else:
                              nproc = nproc[0]
                            otherflags = feffInput.get('feff.MPI.OTHER')[0]
#                           print('executable',executable)
                            print('nnflag',nnflag)
                            print('nnodes',nnodes)
                            print('npflag',npflag)
                            print('nproc',nproc)
                            print('otherflag',otherflags)
                            print(os.path.join(feffdir,exe))
# Create a different version of args
                            args = nnflag + [str(nnodes)] + npflag + [str(nproc)] + otherflags + [os.path.join(feffdir,exe)]
                            print('executable',executable)
                            print('args',args)
#                           executable = feffInput.get('feff.MPI.CMD')[0]
#                           args = feffInput.get('feff.MPI.ARGS',[['']])[0] + [os.path.join(feffdir,exe)]
#                           print('executable',executable)
#                           print('args',args)
#                           sys.exit()
                        else:
                            executable = [os.path.join(feffdir,exe)]
                            args = ['']

                        runExecutable('',exdir,executable,args,out,err)

                outFile=os.path.join(exdir,'xmu.dat')
                output[target] = np.loadtxt(outFile,usecols = (0,3)).T.tolist()
                #print output[target]


            elif (target == 'feffXES'):
                # Set output and error files
                with open(os.path.join(exdir, 'corvus.FEFF.stdout'), 'w') as out, open(os.path.join(exdir, 'corvus.FEFF.stderr'), 'w') as err:

                    # Write input file for FEFF.
                    writeXESInput(feffInput,inpf)

                    # Loop over executable: This is specific to feff. Other codes
                    # will more likely have only one executable. 
                    execs = ['rdinp','atomic','pot','ldos','screen','opconsat','xsph','fms','mkgtr','path','genfmt','ff2x','sfconv']
                    for exe in execs:
                        if 'feff.MPI.CMD' in feffInput:
                            executable = feffInput.get('feff.MPI.CMD')[0]
                            args = feffInput.get('feff.MPI.ARGS',[['']])[0] + [os.path.join(feffdir,exe)]
                        else:
                            executable = [os.path.join(feffdir,exe)]
                            args = ['']

                        runExecutable('',exdir,executable,args,out,err)


                # For this case, I am only passing the directory for now so
                # that other executables in FEFF can use the data.

                outFile=os.path.join(exdir,'xmu.dat')
                output[target] = np.loadtxt(outFile,usecols = (0,3)).T.tolist()

            elif (target == 'feffRIXS'):
                # For RIXS, need to run multiple times as follows.
                
                # Core-Core RIXS
                # 1. Run for the deep core-level.
                # 2. Run for the shallow core-level.
                # 3. Collect files.
                # 4. Run rixs executable.
                
                # Core-valence RIXS
                # 1. Run for the deep core-level.
                # 2. Run for the valence level.
                # 3. Run and XES calculation.
                # 4. Run rixs executable.

                # Set global settings for all runs.
                # Set default energy grid
                setInput(feffInput,'feff.egrid',[['e_grid', -10, 10, 0.05],['k_grid', 'last', 4, 0.025]])
                setInput(feffInput,'feff.exchange',[[0, 0.0, -20.0, 0]])
                setInput(feffInput,'feff.corehole',[['RPA']],Force=True) # maybe not this one. Need 'NONE' for valence

                setInput(feffInput,'feff.edge',[['K','VAL']])
                edges = feffInput['feff.edge'][0]

                # Save original state of input
                savedInput = dict(feffInput)

                # Loop over edges and run XANES for each edge:
                # Save deep edge
                edge0 = edges[0]

                nEdge=0
                for edge in edges:
                    nEdge = nEdge + 1

                    # Make directory
                    dirname = os.path.join(exdir,edge)
                    if not os.path.exists(dirname):
                        os.mkdir(dirname)

                    if edge.upper() != "VAL":
                        outFileName='rixsET.dat'
                        # Delete XES input key
                        if 'feff.xes' in feffInput:
                            del feffInput['feff.xes']
                            
                        # Set edge.
                        setInput(feffInput,'feff.edge',[[edge]],Force=True)

                        # Set icore to other edge
                        setInput(feffInput,'feff.icore',[[getICore(edge0)]],Force=True)

                        # Set corehole RPA
                        setInput(feffInput,'feff.corehole',[['RPA']],Force=True)

                        # Set default energy grid
                        setInput(feffInput,'feff.egrid',[['e_grid', -10, 10, 0.05],['k_grid', 'last', 4, 0.025]])
                        
                        # Set RLPRINT
                        setInput(feffInput,'feff.rlprint',[[True]],Force=True)

                        feffinp = os.path.join(dirname, 'feff.inp')
                        # Write XANES input for this run
                        writeXANESInput(feffInput,feffinp)

                    else: # This is a valence calculation. Calculate NOHOLE and XES
                        # XANES calculation
                        outFileName='rixsET-sat.dat'
                        # Find out if we are using a valence hole for valence calculation.
                        # Set edge.
                        setInput(feffInput,'feff.edge',[[edge0]],Force=True)
                        if len(edges) == nEdge+1:
                            # Set corehole RPA
                            setInput(feffInput,'feff.corehole',[['RPA']],Force=True)

                            # We want to use this core-state as the core hole in place of the valence
                            # Set screen parameters                            
                            setInput(feffInput,'feff.screen',[['icore', getICore(edges[nEdge])]],Force=True)
                        else:
                            # Set corehole NONE
                            setInput(feffInput,'feff.corehole',[['NONE']],Force=True)
                            # Write wscrn.dat to VAL directory
                            wscrnFileW = os.path.join(dirname,'wscrn.dat')
                            writeList(wscrnLines,wscrnFileW)
                            
                        # Set icore to deep edge
                        setInput(feffInput,'feff.icore',[[getICore(edge0)]],Force=True)

                        # Set default energy grid
                        setInput(feffInput,'feff.egrid',[['e_grid', -10, 10, 0.05],['k_grid', 'last', 4, 0.025]])
                        
                        # Set RLPRINT
                        setInput(feffInput,'feff.rlprint',[[True]],Force=True)

                        #
                        # Run XES for the deep level
                        # Save XANES card, but delete from input
                        xanesInput = {}
                        if 'feff.xanes' in feffInput:
                            setInput(xanesInput, 'feff.xanes', feffInput['feff.xanes'])
                            del feffInput['feff.xanes']
                        
                        # Set XES options
                        # Set default energy grid
                        del feffInput['feff.egrid']
                        setInput(feffInput,'feff.egrid',[['e_grid', -40, 10, 0.1]])
                        setInput(feffInput,'feff.xes', [[-20, 10, 0.1]])
                        xesdir = os.path.join(exdir,'XES')
                        if not os.path.exists(xesdir):
                            os.mkdir(xesdir)

                        feffinp = os.path.join(xesdir,'feff.inp')
                        # Write XES input.
                        writeXESInput(feffInput,feffinp)

                        # Run executables to get XES
                        # Set output and error files
                        with open(os.path.join(xesdir, 'corvus.FEFF.stdout'), 'w') as out, open(os.path.join(xesdir, 'corvus.FEFF.stderr'), 'w') as err:
                            execs = ['rdinp','atomic','pot','screen','opconsat','xsph','fms','mkgtr','path','genfmt','ff2x','sfconv']
                            for exe in execs:
                                if 'feff.MPI.CMD' in feffInput:
                                    executable = feffInput.get('feff.MPI.CMD')[0]
                                    args = feffInput.get('feff.MPI.ARGS',[['']])[0] + [os.path.join(feffdir,exe)]
                                else:
                                    executable = [os.path.join(feffdir,exe)]
                                    args = ['']

                                runExecutable('',xesdir,executable,args,out,err)

                        # Make xes.dat from xmu.dat
                        xmuFile = open(os.path.join(xesdir,'xmu.dat'))
                        xesFile = os.path.join(dir,'xes.dat')
                        xesLines=[]
                        for line in xmuFile:
                            if line.lstrip()[0] != '#':
                                fields=line.split()
                                xesLines = xesLines + [str(xmu - float(fields[1])) + '   ' + fields[3]]
                            elif 'Mu' in line:
                                fields = line.strip().split()
                                xmu = float(fields[2][3:len(fields[2])-3])

                        # Write lines in reverse order so that column 1 is sorted correctly.
                        writeList(xesLines[::-1],xesFile)
                        
                        # Now make input file for XANES calculation.
                        if 'feff.xanes' in xanesInput:
                            feffInput['feff.xanes'] = xanesInput['feff.xanes']
                            
                        del feffInput['feff.xes']
                        # Set default energy grid
                        setInput(feffInput,'feff.egrid',[['e_grid', -10, 10, 0.05],['k_grid', 'last', 4, 0.025]],Force=True)
                        feffinp = os.path.join(dirname, 'feff.inp')
                        # Write XANES input for this run
                        writeXANESInput(feffInput,feffinp)

                    
                        
                    # Run XANES for this edge
                    # Set output and error files
                    with open(os.path.join(dirname, 'corvus.FEFF.stdout'), 'w') as out, open(os.path.join(dirname, 'corvus.FEFF.stderr'), 'w') as err:
                        execs = ['rdinp','atomic','pot','screen','opconsat','xsph','fms','mkgtr','path','genfmt','ff2x','sfconv']
                        for exe in execs:
                            if 'feff.MPI.CMD' in feffInput:
                                executable = feffInput.get('feff.MPI.CMD')[0]
                                args = feffInput.get('feff.MPI.ARGS',[['']])[0] + [os.path.join(feffdir,exe)]
                            else:
                                executable = [os.path.join(feffdir,exe)]
                                args = ['']

                            runExecutable('',dirname,executable,args,out,err)


                    # Now copy files from this edge to main directory
                    shutil.copyfile(os.path.join(dirname,'wscrn.dat'), os.path.join(exdir,'wscrn_' + str(nEdge) + '.dat'))
                    shutil.copyfile(os.path.join(dirname,'phase.bin'), os.path.join(exdir,'phase_' + str(nEdge) + '.bin'))
                    shutil.copyfile(os.path.join(dirname,'gg.bin'), os.path.join(exdir,'gg_' + str(nEdge) + '.bin'))
                    shutil.copyfile(os.path.join(dirname,'xsect.dat'), os.path.join(exdir,'xsect_' + str(nEdge) + '.dat'))
                    shutil.copyfile(os.path.join(dirname,'xmu.dat'), os.path.join(exdir,'xmu_' + str(nEdge) + '.dat'))
                    shutil.copyfile(os.path.join(dirname,'rl.dat'), os.path.join(exdir,'rl_' + str(nEdge) + '.dat'))
                    shutil.copyfile(os.path.join(dirname,'.dimensions.dat'), os.path.join(exdir,'.dimensions.dat'))
                    # If this is the first edge, get the screened potential.
                    if nEdge == 1:
                        wscrnLines = []
                        with open(os.path.join(dirname,'wscrn.dat'),'r') as wscrnFileR:
                            for wscrnLine in wscrnFileR.readlines():
                                if wscrnLine.lstrip()[0] == '#':
                                    wscrnLines = wscrnLines + [wscrnLine.strip()]
                                else:
                                    wscrnFields = wscrnLine.strip().split()
                                    wscrnLines = wscrnLines + [wscrnFields[0] + ' 0.0  0.0']
                                
                                
                # Finally, run rixs executable
                feffInput = savedInput
                setInput(feffInput,'feff.rixs', [[0.1, 0.1]])
                
                feffinp = os.path.join(exdir, 'feff.inp')
                # Write XANES input for this run
                writeXANESInput(feffInput,feffinp)

                # Set output and error files                
                with open(os.path.join(exdir, 'corvus.FEFF.stdout'), 'w') as out, open(os.path.join(exdir, 'corvus.FEFF.stderr'), 'w') as err:
                    execs = ['rdinp','atomic','rixs']
                    for exe in execs:
                        if 'feff.MPI.CMD' in feffInput:
                            executable = feffInput.get('feff.MPI.CMD')[0]
                            args = feffInput.get('feff.MPI.ARGS',[['']])[0] + [os.path.join(feffdir,exe)]
                        else:
                            executable = [os.path.join(feffdir,exe)]
                            args = ['']

                        runExecutable('',exdir,executable,args,out,err)

                outFile=os.path.join(exdir,outFileName)
                output[target] = np.loadtxt(outFile).T.tolist()
    
## OPCONS BEGIN
            elif (target == 'opcons'):
        
# Opcons imports
                import copy

                from corvus.controls import generateAndRunWorkflow
                # Define some constants

                hart = 2*13.605698
                alpinv = 137.03598956
                bohr = 0.529177249

# Used in fixing element symbols
                only_alpha = re.compile('[^a-zA-Z]')
                
# Set prefix for sdtout of feff runs.
                runExecutable.prefix = '\t\t\t'
                
# Copy general input to local one
                input2 = copy.deepcopy(input)

# Modify the common values of local input
                input2['feff.setedge'] = input.get('feff.setedge',[[True]])
                input2['feff.absolute'] = [[True]]
                input2['feff.rgrid'] = [[0.01]]

# Copy general config to local one
                config2 = copy.deepcopy(config)

# Set directory to run in.
                config2['cwd'] = config['xcDir']

# Set xcIndexStart to -1 so that xcDir will be set below rather than in prep.
                config2['xcIndexStart'] = -1

# Use absolute units for everything.
                config2['feff.absolute'] = [[True]]

# Initialize variables that collect results (?)
                NumberDensity = []
                vtot = 0.0
                xas_arr = []
                xas_conv_arr = []
                xas0_arr = []
                xas0_conv_arr = []
                en_arr = []
                component_labels = []

# The logic of the lines below is weird: In opcons calculations the absorber is chosen on the fly by looping over all unique atoms
                if 'absorbing_atom' not in input:
                    absorbers = []
                else:
                    absorbers = input['absorbing_atom'][0]

# Build a list of absorbers for the system
# I think this also build a fake cluster to go in the input
                if 'cif_input' in input2:
# Added by FDV
# Trying to replicate the data reading with pymatgen. This is meant to ensure
# consistency between the space group naming in the CIF input generation
# and usage
                    Sys_Str_Parser = CifParser(os.path.abspath(input2['cif_input'][0][0]))
# There should only be 1 structure in the cif file
                    if len(Sys_Str_Parser.get_structures()) > 1:
                      print('More than one structure in CIF file. Check input')
                      sys.exit()
                    Sys_Str = Sys_Str_Parser.get_structures()[0]
# There is no way to generate a sym structure from a read CIF in pymatgen as far
# as I can tell. So we have to symmetrize again
                    spgAna = SpacegroupAnalyzer(Sys_Str)
                    Sys_Str_Sym = spgAna.get_symmetrized_structure()
                    formula = Sys_Str_Sym.formula

                    VTot_fdv = Sys_Str_Sym.volume
# Debug
#                   pp_debug.pprint(dir(Sys_Str_Sym))
#                   pp_debug.pprint(Sys_Str_Sym.volume)
#                   sys.exit()
#                   pp_debug.pprint(Sys_Str_Sym)
#                   pp_debug.pprint(Sys_Str_Sym.equivalent_sites)
#                   print(Sys_Str.get_space_group_info())
#                   print(Sys_Str.get_primitive_structure())
#                   print(Sys_Str.equivalent_sites)

# Now we should have the info to replicate the code below
                    symmult_fdv = [ len(EqSite) for EqSite in Sys_Str_Sym.equivalent_sites ]

# Debug
#                   print(dir(Sys_Str_Sym.equivalent_sites[0][0]))
#                   print(Sys_Str_Sym.equivalent_sites[0][0].specie)
#                   print(Sys_Str_Sym.equivalent_sites[0][0].species_string)
#                   sys.exit()

                    elements_fdv = [ EqSite[0].species_string for EqSite in Sys_Str_Sym.equivalent_sites ]
                    component_labels_fdv = [ Elem+str(iElem+1) for (iElem,Elem) in enumerate(elements_fdv) ]
                    cluster_fdv = [ ['Cu', 0.0, 0.0, iSite*2.0 ] for iSite in range(len(Sys_Str_Sym.equivalent_sites)) ]

                    if 'absorbing_atom' not in input:
                      absorbers_fdv = list(range(1,len(Sys_Str_Sym.equivalent_sites)+1))

# Debug
                    print(symmult_fdv)
                    print(cluster_fdv)
                    print(elements_fdv)
                    print(component_labels_fdv)
                    print(absorbers_fdv)
#                   sys.exit()

#                   cifFile = ReadCif(os.path.abspath(input2['cif_input'][0][0]))
#                   cif_dict = cifFile[list(cifFile.keys())[0]]
#                   cell_data = CellData()
#                   cell_data.getFromCIF(cif_dict)
#                   print(cell_data.atomdata)
#                   cell_data.primitive()
#                   print(VTot_fdv/(bohr**3))
#                   print(cell_data.volume()*(cell_data.lengthscale/bohr)**3)
#                   sys.exit()
#                   print(cell_data.atomdata)
#                   symmult = []
#                   cluster = []
                
#                   for ia,a in enumerate(cell_data.atomdata): # This loops over sites in the original cif
#                       symmult = symmult + [len(a)]
#                       element = list(a[0].species.keys())[0]
#                       component_labels = component_labels + [element + str(ia+1)]
#                       if 'absorbing_atom' not in input:
#                           absorbers = absorbers + [ia+1]
#                       cluster = cluster + [['Cu', 0.0, 0.0, ia*2.0 ]]

# Debug
#                   pp_debug.pprint(symmult)
#                   pp_debug.pprint(cluster)
#                   pp_debug.pprint(component_labels)
#                   pp_debug.pprint(absorbers)
#                   sys.exit()

# Replace the original cif2cell generated data with my own generated with
# pymatgen
                    symmult = symmult_fdv
                    cluster = cluster_fdv
                    elements = elements_fdv
                    component_labels = component_labels_fdv
                    absorbers = absorbers_fdv
#                   print('element',element)

                    if 'cluster' not in input2:    
                        input2['cluster'] = cluster
 
# Debug: FDV
#               print('ABSORBERS')
#               pp_debug.pprint(absorbers)

# OPCONS LOOP SETUP BEGIN -------------------------------------------------------------------------------------
# Added by FDV
# Creating a list to collect the inputs for delayed execution
                WF_Params_Dict = {}

# Determine how many processors we can dedicate to each edge run
                OC_NP_Tot = input2['feff.MPI.NPTOT'][0][0]
                OC_Tot_Runs = 0
#               for (element,absorber) in zip(elements,absorbers):
                for (element,absorber) in zip(elements_fdv,absorbers):
                  for edge in feff_edge_dict[only_alpha.sub('',element)]:
                    OC_Tot_Runs = OC_Tot_Runs + 1

# Debug
                print(' FDV: OC_NP_Tot, OC_Tot_Runs')
                print(OC_NP_Tot, OC_Tot_Runs)
#               sys.exit()

                NP_Base = OC_NP_Tot//OC_Tot_Runs
                NP_Extra = OC_NP_Tot%OC_Tot_Runs
                NP_Run_Partition = OC_Tot_Runs*[NP_Base]
                for iRun in range(NP_Extra):
                  NP_Run_Partition[iRun] += 1

                print('NP_Run_Partition',NP_Run_Partition)
#               sys.exit()

# For each atom in absorbing_atoms run a full-spectrum calculation (all edges, XANES + EXAFS)
                iRun_Count = 0
#               for absorber in absorbers:
                for (element,absorber) in zip(elements,absorbers):

                    print('')
                    print("##########################################################")
                    print("       Component: " + component_labels[absorber-1])
                    print("##########################################################")
                    print('')

### BEGIN INPUT GEN --------------------------------------------------------------------------------------------
                    input2['absorbing_atom'] = [[absorber]]
                    input2['feff.target'] = [[absorber]]

                    if 'cif_input' in input2:
                        input2['feff.target'] = [[absorber]]
#                       element = list(cell_data.atomdata[absorber-1][0].species.keys())[0]
                        if 'number_density' not in input:
                            NumberDensity = NumberDensity + [symmult[absorber - 1]]

                    else:
# This only works if all elements are treated as the same for our calculation
                        element = input['cluster'][absorber-1][0]
                        if 'number_density' not in input:
                            n_element = 0
                            for atm in input['cluster']:
                                if element in atm:
                                    n_element += 1
                                    
                            NumberDensity = NumberDensity + [n_element]
                    
                    print('Number in unit cell: ' + str(NumberDensity[-1]))

### END INPUT GEN --------------------------------------------------------------------------------------------

                    # For each edge for this atom, run a XANES and EXAFS run

# Commented out by FDV, unused, simplifying
#                   FirstEdge = True
                    Item_Absorber = {}
                    for iedge,edge in enumerate(feff_edge_dict[only_alpha.sub('',element)]):

                        Item_Edge = {}
                        print("\t" + edge)
                        print("\t\t" + 'XANES')

### BEGIN INPUT GEN --------------------------------------------------------------------------------------------
                        input2['feff.edge'] = [[edge]]
                        
# Run XANES 
                        input2['target_list'] = [['feffXANES']]

# Set energy grid for XANES.
                        input2['feff.egrid'] = [['e_grid', -10, 10, 0.1], ['k_grid','last',5,0.07]]
                        input2['feff.control'] = [[1,1,1,1,1,1]]

                        config2['xcDir'] = os.path.join(config2['cwd'],component_labels[absorber-1],edge,'XANES')

                        # NOTE: This will unbalance the load. Need to figure out best way to do this.
                        if iedge == 0 and element == elements[0]:
                          # Run LDOS only for first edge of first absorber. Others will be copied.
                          input2['feff.ldos'] = [[-40.0, 0.0, 0.5,160]]
                          dos_dir = config2['xcDir']
                            
                        targetList = [['feffXANES']]
                        if 'feff.scf' in input:
                            input2['feff.scf'] = input['feff.scf']
                        else:
                            input2['feff.scf'] = [[4.0,0,100,0.1,0]]
                        if 'feff.fms' in input:
                            input2['feff.fms'] = input['feff.fms']
                        else:
                            input2['feff.fms'] = [[6.0]]
                            
                        input2['feff.rpath'] = [[0.1]]
### END INPUT GEN --------------------------------------------------------------------------------------------

# Added by FDV
                        print(' FDV: XANES INPUT')
# Adjust the number of processors defined in the input so we only use those
# assigned to this run
                        input2['feff.MPI.NP'] = [[NP_Run_Partition[iRun_Count]]]
                        print('XXX input2 XXX')
                        pp_debug.pprint(input2)
                        Item_xanes = { 'config2':copy.deepcopy(config2),
                                       'input2':copy.deepcopy(input2),
                                       'targetList':copy.deepcopy(targetList) }

# Commented out by FDV, unused, simplifying
#                       FirstEdge = False

### BEGIN INPUT GEN --------------------------------------------------------------------------------------------
                        print("\t\t" + 'EXAFS')

                        xanesDir = config2['xcDir']
                        exafsDir = os.path.join(config2['cwd'],component_labels[absorber-1],edge,'EXAFS')
                        config2['xcDir'] = exafsDir
                        
                        # Delete ldos card if it exists.
                        if 'feff.ldos' in input2:
                          del input2['feff.ldos']
                        
                        input2['feff.control'] = [[0, 1, 1, 1, 1, 1]]
                        input2['feff.egrid'] = [['k_grid', -20, -2, 1], ['k_grid',-2,0,0.07], ['k_grid', 0, 40, 0.07],['exp_grid', 'last', 500000.0, 10.0]]
                        if 'feff.fms' in input:
                            input2['feff.rpath'] = [[max(input['feff.fms'][0][0],0.1)]]
                        else:
                            input2['feff.rpath'] = [[6.0]]
                        input2['feff.fms'] = [[0.0]]
### END INPUT GEN --------------------------------------------------------------------------------------------

# Added by FDV
                        print(' FDV: EXAFS INPUT')
# Adjust the number of processors defined in the input so we only use those
# assigned to this run
# For the EXAFS part we force a single processor
                        input2['feff.MPI.NP'] = [[1]]
                        pp_debug.pprint(input2)
                        iRun_Count += 1
                        Item_exafs = { 'config2':copy.deepcopy(config2),
                                       'input2':copy.deepcopy(input2),
                                       'targetList':copy.deepcopy(targetList) }

                        Item_Absorber[edge] = { 'xanes':Item_xanes,
                                                'exafs':Item_exafs }
                    print('')
                    WF_Params_Dict[absorber] = Item_Absorber
                print('')
# OPCONS LOOP SETUP END ---------------------------------------------------------------------------------------

# Debug: FDV
                print('#### FDV ####')
                print('#### All WF Params ####')
                pp_debug.pprint(WF_Params_Dict)
# Monty has issue on 2.7, so will just use pickle
                import pickle
                pickle.dump(WF_Params_Dict,open('WF_Params_Dict.pickle','wb'))
# Debug
#               sys.exit()

# OPCONS LOOP RUN BEGIN ---------------------------------------------------------------------------------------

# For each atom in absorbing_atoms run a full-spectrum calculation (all edges, XANES + EXAFS)
# Splitting the run loop into two, run all the XANES calcs first, then do the EXAFS. This might
# simplify things a bit since the EXAFS runs wil only be done in 1 core.
# XANES LOOP

# Attempting to multithread the calls
# For now we launch all threads to run with whatever number of cores are in the input and not try to do
# anything smart. We will just put enough cores available in the overall run.
                nTasks_Est = len(absorbers)*len(list(WF_Params_Dict[absorber].keys()))
                nTasks = 0
                Tasks = []
                for absorber in absorbers:
                  for edge in WF_Params_Dict[absorber].keys():

                    config2    = WF_Params_Dict[absorber][edge]['xanes']['config2']
                    input2     = WF_Params_Dict[absorber][edge]['xanes']['input2']
                    targetList = WF_Params_Dict[absorber][edge]['xanes']['targetList']
                    if 'opcons.usesaved' not in input:
#                     generateAndRunWorkflow(config2,input2,targetList)
#                     thrd.Thread(target=generateAndRunWorkflow,args=(config2,input2,targetList)).start()
                      Prcs = mltp.Process(target=generateAndRunWorkflow,args=(config2,input2,targetList))
                      Tasks.append(Prcs)
                      nTasks += 1
                    else:
                    # Run if xmu.dat doesn't exist.
                      if not os.path.exists(os.path.join(config2['xcDir'],'xmu.dat')):
#                       generateAndRunWorkflow(config2,input2,targetList)
#                       thrd.Thread(target=generateAndRunWorkflow,args=(config2,input2,targetList)).start()
                        Prcs = mltp.Process(target=generateAndRunWorkflow,args=(config2,input2,targetList))
                        Tasks.append(Prcs)
                        nTasks += 1
                      else:
                        print("\t\t\txmu.dat already calculated. Skipping.")

# Launch all the tasks
                for Tsk in Tasks:
                    Tsk.start()

# Ensure all of the tasks have finished
                for Tsk in Tasks:
                    Tsk.join()

# Monitor the number of threads until they are all done..
#               while thrd.active_count()>1:
#                 print('Waiting on',thrd.active_count()-1,'thread(s) of',nThreads)
#                 time.sleep(10)

# EXAFS LOOP
                nTasks = 0
                Tasks = []
                for absorber in absorbers:
                  for edge in WF_Params_Dict[absorber].keys():

                    config2    = WF_Params_Dict[absorber][edge]['xanes']['config2']
                    xanesDir = config2['xcDir']
                    exafsDir = os.path.join(config2['cwd'],component_labels[absorber-1],edge,'EXAFS')
                    if not os.path.exists(exafsDir):
                      os.makedirs(exafsDir)
                    shutil.copyfile(os.path.join(xanesDir,'apot.bin'), os.path.join(exafsDir,'apot.bin'))
                    shutil.copyfile(os.path.join(xanesDir,'pot.bin'), os.path.join(exafsDir,'pot.bin'))

# Modified by FDV
# Commented out and moved to an independent loop
                    config2    = WF_Params_Dict[absorber][edge]['exafs']['config2']
                    input2     = WF_Params_Dict[absorber][edge]['exafs']['input2']
                    targetList = WF_Params_Dict[absorber][edge]['exafs']['targetList']
                    if 'opcons.usesaved' not in input2:
#                     generateAndRunWorkflow(config2,input2,targetList)
#                     thrd.Thread(target=generateAndRunWorkflow,args=(config2,input2,targetList)).start()
                      Prcs = mltp.Process(target=generateAndRunWorkflow,args=(config2,input2,targetList))
                      Tasks.append(Prcs)
                      nTasks += 1
                    else:
                      # Run if xmu.dat doesn't exist.
                      if not os.path.exists(os.path.join(config2['xcDir'],'xmu.dat')):
#                       generateAndRunWorkflow(config2,input2,targetList)
#                       thrd.Thread(target=generateAndRunWorkflow,args=(config2,input2,targetList)).start()
                        Prcs = mltp.Process(target=generateAndRunWorkflow,args=(config2,input2,targetList))
                        Tasks.append(Prcs)
                        nTasks += 1

# Launch all the tasks
                for Tsk in Tasks:
                    Tsk.start()

# Ensure all of the tasks have finished
                for Tsk in Tasks:
                    Tsk.join()

# Monitor the number of threads until they are all done..
#               while thrd.active_count()>1:
#                 print('Waiting on',thrd.active_count()-1,'thread(s) of',nThreads)
#                 time.sleep(10)

# OPCONS LOOP RUN END -----------------------------------------------------------------------------------------

# OPCONS LOOP ANA BEGIN ---------------------------------------------------------------------------------------
# For each atom in absorbing_atoms run a full-spectrum calculation (all edges, XANES + EXAFS)
                for iabs,absorber in enumerate(absorbers):

                    print('')
                    print("##########################################################")
                    print("       Component: " + component_labels[absorber-1])
                    print("##########################################################")
                    print('')

# Commented out by FDV, unused, simplifying
#                   FirstEdge = True
                    for iedge,edge in enumerate(WF_Params_Dict[absorber].keys()):

                        print("\t" + edge)
                        print("\t\t" + 'XANES')

                        
# Added by FDV
                        config2    = WF_Params_Dict[absorber][edge]['xanes']['config2']
                        input2     = WF_Params_Dict[absorber][edge]['xanes']['input2']
                        targetList = WF_Params_Dict[absorber][edge]['xanes']['targetList']
### BEGIN OUTPUT ANA --------------------------------------------------------------------------------------------
                        if 'cif_input' in input:
                            # Get total volume from cif in atomic units. 
#                           vtot = cell_data.volume()*(cell_data.lengthscale/bohr)**3
                            vtot = VTot_fdv/(bohr**3)
                        else:
                            # Get norman radii from xmu.dat
                            with open(os.path.join(config2['xcDir'],'xmu.dat')) as f:
                                for line in f: # Go through the lines one at a time
                                    words = line.split()
                                    if 'Rnm=' in words:
                                        vtot = vtot + (float(words[words.index('Rnm=')+1])/bohr)**3*4.0/3.0*np.pi
                                        break
                                  
                            f.close()

                        outFile = os.path.join(config2['xcDir'],'xmu.dat')
                        e1,k1,xanes = np.loadtxt(outFile,usecols = (0,2,3)).T
#                       print('k1',k1)
#                       print('xanes fdv',xanes)
                        xanes = np.maximum(xanes,0.0)
### END OUTPUT ANA --------------------------------------------------------------------------------------------

# Added by FDV
                        config2    = WF_Params_Dict[absorber][edge]['exafs']['config2']
                        input2     = WF_Params_Dict[absorber][edge]['exafs']['input2']
                        targetList = WF_Params_Dict[absorber][edge]['exafs']['targetList']
### BEGIN OUTPUT ANA --------------------------------------------------------------------------------------------
                        outFile = os.path.join(config2['xcDir'],'xmu.dat')
                        e2,ep,k2,exafs,mu0 = np.loadtxt(outFile,usecols = (0,1,2,3,4)).T
# Debug: FDV
# Adding a little call to temporarily fix the NaN issues
                        (exafs, mu0) = Temp_Fix_NaN(k2,exafs,mu0)
                        print('fdv e2',e2)
                        print('fdv k2',k2)
                        print('xcDir',config2['xcDir'])
                        print('fdv exafs',exafs)
                        print('fdv mu0',mu0)
                        exafs = np.maximum(exafs,0.0)
                        mu0 = np.maximum(mu0,0.0)
                        e0 = e2[100] - (k2[100]*bohr)**2/2.0*hart
                        
                        # Interpolate onto a union of the two energy-grids and smoothly go from one to the other between  
                        e_tot = np.unique(np.append(e1,e2))
                        k_tot = np.where(e_tot > e0, np.sqrt(2.0*np.abs(e_tot-e0)/hart), -np.sqrt(2.0*np.abs(e0 - e_tot)/hart))/bohr
                        kstart = 3.0
                        kfin = 4.0
                        weight1 = np.cos((np.minimum(np.maximum(k_tot,kstart),kfin)-kstart)/(kfin-kstart)*np.pi/2)**2
                        weight2 = 1.0 - weight1
                        xas_element = NumberDensity[iabs]*(np.interp(e_tot,e1,xanes)*weight1 + np.interp(e_tot,e2,exafs)*weight2)
                        # Zero out portions below 1/5th of the edge energy.
                        xas_element[np.where(e_tot < 0.2*e0)] = 0.0
                        xas0_element = NumberDensity[iabs]*np.interp(e_tot,e2,mu0)
                        xas0_element[np.where(e_tot < 0.2*e0)] = 0.0
                        
                        # Find EFermi:
                        # Find mu, where k = 0
                        ind_Fermi = np.argmin(k2**2)
                        EFermi = ep[ind_Fermi]
                        
                        
                        # Load pot.inp and read to get info on ldos files 
                        infl = open(os.path.join(config2['xcDir'],'pot.inp'), 'r')
                        Lines = infl.readlines()
                        infl.close()

                        # second field of second line is the number of unique potentials. 
                        second = Lines[1].split()
                        npot = int(second[1])

                        # Which dos to use?  Will take this from ihole in pot.inp.
                        ihole = int(second[3])

                        # We also need the target of this particular run, since the ldos will not necessarily have come
                        # from this run. NOTE: This can only really be used by opcons at present. 
                        infl = open(os.path.join(config2['xcDir'],'feff.inp'), 'r')
                        feffInpLines = infl.readlines()
                        infl.close()
                        for line in feffInpLines:
                          if 'target' in line.lower():
                            target_atom = int(line.split()[1])

                        # If this is the edge and absorber with the LDOS in it, load LDOS.
                        if iedge == 0 and iabs == 0:
                          ipot = 0
                          ldos_files = []
                          dos_array = []
                          # Get the ldos.
                          while ipot <= npot:
                            ldos_files = ldos_files + ['ldos' + '{:02}'.format(ipot) + '.dat']
                            dos_array = dos_array + [np.loadtxt(os.path.join(dos_dir,ldos_files[ipot])).T]
                            ipot += 1

                        # Use the dos from this target absorber, and ldos corresponding to the correct symmetry of the core-level.
                        idos = getHoleSymm(ihole-1)+1
                        w = dos_array[target_atom][0]
                        dos = dos_array[target_atom][idos]
                        print('Edge and absorber:', edge, absorber)
                        print('LDOS file:', ldos_files[target_atom])
                        print('Angular momentum', idos - 1)
                        # Below the first point of the XANES grid, use the background of the EXAFS calculation.
                        xas_element[np.where(e_tot < e1[0])] = xas0_element[np.where(e_tot < e1[0])] 


                        e_cv = 50.0 # This should be adjustable in general.
                        # If edge energy is less than e_cv, run LDOS convolution on xas_element and xas0_element
                        if e0 < e_cv:
                          xas_conv_element = dos_conv(e_tot, EFermi, e0, k_tot, xas_element, w, dos)
                          xas0_conv_element = dos_conv(e_tot, EFermi, e0, k_tot, xas0_element, w, dos)
                        else:
                          xas_conv_element = xas_element
                          xas0_conv_element = xas0_element

                        np.savetxt(os.path.join(config2['xcDir'],'xmu_conv.dat'),np.array([e_tot,xas_conv_element]).T)
                        xas_arr = xas_arr + [xas_element]
                        xas_conv_arr = xas_conv_arr + [xas_conv_element]
                        xas0_arr = xas0_arr + [xas0_element]
                        xas0_conv_arr = xas0_conv_arr + [xas0_conv_element]
                        en_arr = en_arr + [e_tot]
                        #plt.plot(e_tot, xas_element)
                        #plt.show()
### END OUTPUT ANA --------------------------------------------------------------------------------------------
                    print('')
                print('')
# OPCONS LOOP ANA END -----------------------------------------------------------------------------------------

# POST LOOP ANALYSYS: If everything is correct we should not have to change anything below
                # Interpolate onto common grid from 0 to 500000 eV
                # Make common grid as union of all grids.
                energy_grid = np.unique(np.concatenate(en_arr))

                # Now loop through all elements and add xas from each element
                xas_tot = np.zeros_like(energy_grid)
                xas_conv_tot = np.zeros_like(energy_grid)
                xas0_tot = np.zeros_like(energy_grid)
                xas0_conv_tot = np.zeros_like(energy_grid)
                for i,en in enumerate(en_arr):
                    xas_tot = xas_tot + np.interp(energy_grid,en,xas_arr[i],left=0.0,right=0.0)
                    xas_conv_tot = xas_conv_tot + np.interp(energy_grid,en,xas_conv_arr[i],left=0.0,right=0.0)
                    xas0_tot = xas0_tot + np.interp(energy_grid,en,xas0_arr[i],left=0.0,right=0.0)
                    xas0_conv_tot = xas0_conv_tot + np.interp(energy_grid,en,xas0_conv_arr[i],left=0.0,right=0.0)
                    
                xas_tot = xas_tot/vtot
                xas_conv_tot = xas_conv_tot/vtot
                xas0_tot = xas0_tot/vtot
                xas0_conv_tot = xas0_conv_tot/vtot

                # transform to eps2. xas_tot*-4pi/apha/\omega*bohr**2
                energy_grid = energy_grid/hart
                eps2 = xas_tot*4*np.pi*alpinv*bohr**2/energy_grid
                eps2 = eps2[np.where(energy_grid > 0)]
                eps2_conv = xas_conv_tot*4*np.pi*alpinv*bohr**2/energy_grid
                eps2_conv = eps2_conv[np.where(energy_grid > 0)]
                eps2_bg = xas0_tot*4*np.pi*alpinv*bohr**2/energy_grid
                eps2_bg = eps2_bg[np.where(energy_grid > 0)]
                eps2_conv_bg = xas0_conv_tot*4*np.pi*alpinv*bohr**2/energy_grid
                eps2_conv_bg = eps2_conv_bg[np.where(energy_grid > 0)]
 
                energy_grid = energy_grid[np.where(energy_grid > 0)]
                #plt.plot(energy_grid,eps2)
                #plt.show()

                if False:
                    # Test with Lorentzian
                    eps2 = -5.0/((energy_grid - 277.0)**2 + 5.0**2) + 5.0/((energy_grid + 277.0)**2 + 5.0**2) 
                
                # Perform KK-transform
                print('Performaing KK-transform of eps2:')
                print('')

                # Background
                w,eps1_bg = kk_transform(energy_grid, eps2_bg)
                # Add Drude term if requested
                if 'opcons.drude' in input:
                    wp = input['opcons.drude'][0][0]/hart
                    Gamma = input['opcons.drude'][0][1]/hart
                    print('wp,gamma',wp, Gamma)
                    eps_drude = - wp**2/(w**2 + 1j*Gamma*w)
                    eps1_drude = np.real(eps_drude)
                    eps2_drude = np.imag(eps_drude)

                eps2_bg = np.interp(w,energy_grid,eps2_bg) + eps2_drude
                eps1_bg = eps1_bg + 1.0 + eps1_drude
                eps_bg = eps1_bg + 1j*eps2_bg

                # Background with LDOS convolution.
                w,eps1_conv_bg = kk_transform(energy_grid, eps2_conv_bg)
                eps2_conv_bg = np.interp(w,energy_grid,eps2_conv_bg) + eps2_drude
                eps1_conv_bg = eps1_conv_bg + 1.0 + eps1_drude
                eps_conv_bg = eps1_conv_bg + 1j*eps2_conv_bg

                # With fine-structure
                w,eps1 = kk_transform(energy_grid, eps2)
                eps2 = np.interp(w,energy_grid,eps2) + eps2_drude
                eps1 = eps1 + 1.0 + eps1_drude
                eps = eps1 + 1j*eps2

                # With fine-structure and LDOS convolution.
                w,eps1_conv = kk_transform(energy_grid, eps2_conv)
                eps2_conv = np.interp(w,energy_grid,eps2_conv) + eps2_drude
                eps1_conv = eps1_conv + 1.0 + eps1_drude
                eps_conv = eps1_conv + 1j*eps2_conv

                
                # Transform to optical constants
                index_of_refraction = np.sqrt(eps)
                index_of_refraction_bg = np.sqrt(eps_bg)
                reflectance = np.abs((index_of_refraction-1)/(index_of_refraction+1))**2
                reflectance_bg = np.abs((index_of_refraction_bg-1)/(index_of_refraction_bg+1))**2
                absorption = 2*w*1.0/alpinv*np.imag(index_of_refraction)/bohr*1000
                absorption_bg = 2*w*1.0/alpinv*np.imag(index_of_refraction_bg)/bohr*1000
                energy_loss = -1.0*np.imag(1.0/eps)
                energy_loss_bg = -1.0*np.imag(eps_bg**(-1))

                # With LDOS convolution
                index_of_refraction_conv = np.sqrt(eps_conv)
                index_of_refraction_conv_bg = np.sqrt(eps_conv_bg)
                reflectance_conv = np.abs((index_of_refraction-1)/(index_of_refraction+1))**2
                reflectance_conv_bg = np.abs((index_of_refraction_conv_bg-1)/(index_of_refraction_conv_bg+1))**2
                absorption_conv = 2*w*1.0/alpinv*np.imag(index_of_refraction_conv)/bohr*1000
                absorption_conv_bg = 2*w*1.0/alpinv*np.imag(index_of_refraction_conv_bg)/bohr*1000
                energy_loss_conv = -1.0*np.imag(1.0/eps_conv)
                energy_loss_conv_bg = -1.0*np.imag(eps_conv_bg**(-1))

                # Calculate sumrules
                # eps2 sumrule: V/(2pi^2 N) * \int_0^{\omega} d\omega' \omega' eps2(\omega')
                # Find the number of chemical units
                iUnit = 2
                Number_Of_Units = 1
                while iUnit <= min(NumberDensity):
                    formula = np.array(NumberDensity)/float(iUnit)
                    print('formula', formula)
                    print(np.all(np.mod(formula,1) == 0))
                    if int(min(NumberDensity)/iUnit) == 0:
                        print('break')
                        break
                    elif np.all(np.mod(formula,1) == 0):
                        Number_Of_Units = iUnit
                    iUnit += iUnit

                print("Number of Chemical Formula Units in Unit Cell: ", Number_Of_Units)
                print("Number Densities:", NumberDensity, np.sum(np.array(NumberDensity)))
                #Total_NumberDensity = np.sum(np.array(NumberDensity))/vtot
                Total_NumberDensity = Number_Of_Units/vtot
                n_eff = 1.0/Total_NumberDensity*scipy.integrate.cumulative_trapezoid(w*eps2, x=w, initial=0.0)
                n_eff_conv = 1.0/Total_NumberDensity*scipy.integrate.cumulative_trapezoid(w*eps2_conv, x=w, initial=0.0)
                
                print('Sumrule gives: ',1.0/Total_NumberDensity*np.trapz(w*eps2,w)/(2.0*np.pi**2))
                print('Sumrule with convolution gives: ', 1.0/Total_NumberDensity*np.trapz(w*eps2_conv,w)/(2.0*np.pi**2))
                n_eff = n_eff/(2.0*np.pi**2)
                n_eff_conv = n_eff_conv/(2.0*np.pi**2)
                print("Number Density: ",Total_NumberDensity)

                w = w*hart
                np.savetxt('epsilon.dat',np.array([w,eps1,eps2,eps1_bg,eps2_bg]).transpose())
                np.savetxt('epsilon_drude.dat',np.array([w,eps1_drude,eps2_drude]).transpose())
                np.savetxt('index.dat',np.array([w,np.real(index_of_refraction),np.imag(index_of_refraction),np.real(index_of_refraction_bg),np.imag(index_of_refraction_bg)]).transpose())
                np.savetxt('reflectance.dat',np.array([w,reflectance,reflectance_bg]).transpose())
                np.savetxt('absorption.dat',np.array([w,absorption,absorption_bg]).transpose())
                np.savetxt('loss.dat', np.array([w,energy_loss,energy_loss_bg]).transpose())
                np.savetxt('sumrules.dat', np.array([w,n_eff]).transpose())
                np.savetxt('epsilon_conv.dat',np.array([w,eps1_conv,eps2_conv,eps1_conv_bg,eps2_conv_bg]).transpose())
                np.savetxt('index_conv.dat',np.array([w,np.real(index_of_refraction_conv),np.imag(index_of_refraction_conv),np.real(index_of_refraction_conv_bg),np.imag(index_of_refraction_conv_bg)]).transpose())
                np.savetxt('reflectance_conv.dat',np.array([w,reflectance_conv,reflectance_conv_bg]).transpose())
                np.savetxt('absorption_conv.dat',np.array([w,absorption_conv,absorption_conv_bg]).transpose())
                np.savetxt('loss_conv.dat', np.array([w,energy_loss_conv,energy_loss_conv_bg]).transpose())
                np.savetxt('sumrules_conv.dat', np.array([w,n_eff_conv]).transpose())

                output[target] = (w,eps)

                print('Finished with calculation of optical constants.')
## OPCONS END

# Added by FDV
            elif (target == 'opt_dynmat_s2_exafs'):

# Before we can call the rest of the input preparation, we need to get the
# optmized structure and put it in the "cluster" input, so the stuff below
# can work OK.
# JK - in future, would like to replace this with a "helper" handler that translates
# between generic properties. Thus cluster would be produced by the "helper" handler 
# from cell_struc_xyz_red. 
                input['cluster'] = input['cell_struc_xyz_red']

# First generate any data that is needed from input
                feffInput['feff.atoms']      = getFeffAtomsFromCluster(input)
                feffInput['feff.potentials'] = getFeffPotentialsFromCluster(input)


# Set output and error files
                with open(os.path.join(exdir, 'corvus.FEFF.stdout'), 'w') as out, open(os.path.join(exdir, 'corvus.FEFF.stderr'), 'w') as err:

# Write input file for FEFF.
                  writeEXAFSDMDWInput(feffInput,inpf)

# Before running, we need to write the dym file to be used in this run
                  dymFilename = 'corvus.dym'
                  from corvus.dmdw import writeDym
                  writeDym(input['opt_dynmat'], os.path.join(exdir, dymFilename))
#                 sys.exit()

# Loop over executable: This is specific to feff. Other codes
# will more likely have only one executable. Here I am running 
# rdinp again since writeSCFInput may have different cards than

# Run rdinp and atomic part of calculation
                  execs = ['rdinp','atomic','pot','screen','opconsat','xsph','fms','mkgtr','path','genfmt','ff2x','sfconv']
                  for exe in execs:
                    if 'feff.MPI.CMD' in feffInput:
                        executable = feffInput.get('feff.MPI.CMD')[0]
                        args = feffInput.get('feff.MPI.ARGS',[['']])[0] + [os.path.join(feffdir,exe)]
                    else:
                        executable = [os.path.join(feffdir,exe)]
                        args = ['']
                    runExecutable('',exdir,executable,args,out,err)

                output[target] = exdir

    @staticmethod
    def cleanup(config):
        pass



##### Generic Helper Methods ##########

def check(input, token, default=None):
    if token in input:
        return input[token]
    else:
        return default

def setInput(input, token, default, Force=False):
    # If token is already defined, leave it unless force, otherwise define with default.
    if token not in input or Force:
        input[token] = default

def writeInput(input,inpfile):
    lines=[]
    for key in input:
        if key != 'feff.end' and not key.startswith('feff.MPI') and not key.startswith('feff.lfms'):
            lines = lines + getInpLines(input,key) 
        
    lines = lines + ['END']
    
    # Print feff input file
    writeList(lines, inpfile)
    
def getInpLines(input,token):
    lines=[]
    key = token[len('feff.'):]

    if key == 'end':
        return lines  # Don't use this to write END card, since that might be
                      # written before end of input file.
    
    if token in input:
        # If the first element is not a boolean, this contains values
        # to be stored after keyword.
        if not isinstance(input[token][0][0],bool): # Keywords that have no arguments are written if bool is true
            if token in input: # Use pre-defined values
                for element in input[token]: # Takes care of single and multi-line input.
                    lines.append(' '.join([str(value) for value in element])) 
            else: # Use default
                lines.append(' '.join(default))
          	
            if key in ['atoms','potentials','egrid']: # Some keywords only have arguments starting on the 
                lines.insert(0,key.upper())           # next line.
                
            else:                                     # Most have arguments on the same line as keyword.
                lines[0] = key.upper() + ' ' + lines[0]
                
        elif input[token][0][0]: 
            lines.append(key.upper())
    else:
        if input[token][0][0]:
            lines.append(key.upper()) # Keywords with no arguments.

    # Add a blank line after each line
    lines.append('')

    return lines

def writeList(lines, filename):
    with open(filename, 'w') as f:
        f.write('\n'.join(lines))

def getICore(edge):
    iCore = {'K':1, 'L1':2, 'L2':3, 'L3':4, 'M1':5, 'M2':6, 'M3':7, 'M4':8, 'M5':9}
    if edge in iCore:
        return iCore[edge]
    else:
        print("###########################################")
        print("###########################################")
        print("Error: Unknown edge name.")
        print("###########################################")
        print("###########################################")
        sys.exit()


def runExecutable(execDir,workDir,executable, args,out,err):
    # Runs executable located in execDir from working directory workDir.
    # Tees stdout to file out in real-time, and stderr to file err.
    print('--- FDV ---',args)
    print((runExecutable.prefix + 'Running exectuable: ' + executable[0] + ' ' + ' '.join(args)))
    # Modified by FDV:
    # Adding the / to make the config more generic
    # Modified by JJK to use os.path.join (even safer than above).
    execList = [os.path.join(execDir,executable[0])] + args
    p = subprocess.Popen(execList, bufsize=0, cwd=workDir, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8')
    while True:
        output = p.stdout.readline()
        error = p.stderr.readline()
        if output == '' and p.poll() is not None:
            break
        if output:
            print(output.strip())
            out.write(output.strip() + os.linesep)
        if error:
            print(error.strip())
            err.write(error.strip() + os.linesep)
    rc = p.poll()
    #while True:
    #    pout = p.stdout.readline()
    #    if pout == '' and p.poll() is not None:
    #        break
    #    if pout:
    #        print((runExecutable.prefix + pout.strip()))
    #        out.write(pout)

    #while True:
    #    perr = p.stderr.readline()
    #    if perr == '' and p.poll() is not None:
    #        break
    #    if perr:
    #        print((runExecutable.prefix + '###################################################'))
    #        print((runExecutable.prefix + '###################################################'))
    #        print((runExecutable.prefix + 'Error in executable: ' + executable[0]))
    #        print((runExecutable.prefix + perr.strip()))
    #        print((runExecutable.prefix + '###################################################'))
    #        print((runExecutable.prefix + '###################################################'))
    #        err.write(perr)


    #p.wait()
    
runExecutable.prefix = ''    
    
def readColumns(filename, columns=[1,2]):
    # Read file and clear out comments
    with open(filename, 'r') as file:
        cleanStr = file.read()
    comments = pp.ZeroOrMore(pp.pythonStyleComment).setParseAction(pp.replaceWith(''))
    try:
        cleanStr = comments.transformString(cleanStr)
    except pp.ParseException as pe:
        print(('Parsing Error using pyparsing: invalid input:', pe))
        sys.exit()
    # Define grammar for ncols of data based on number of entries in first row
    floating = pp.Word(pp.nums + ".+-E").setParseAction(lambda t: float(t[0]))
    EOL = pp.LineEnd().suppress()
    row1entry = floating.copy().setWhitespaceChars(" \t")
    row1 = pp.Group(pp.ZeroOrMore(EOL) + pp.OneOrMore(row1entry) + EOL)
    row = pp.Forward()
    def defineTotalCols(toks):
        ncols = len(toks[0])
        row << pp.Group(floating * ncols)
        return None
    row1.addParseAction(defineTotalCols)
    text = row1 + pp.ZeroOrMore(row)
    try:
        data = text.parseString(cleanStr).asList()
    except pp.ParseException as pe:
        print(('Parsing Error using pyparsing: invalid input:', pe))
        sys.exit()
    cols = list(map(list, list(zip(*data))))
    return [cols[i-1] for i in columns]

#### Specific Helper Methods

# Periodic Table information from FEFF
ptable = {}
atomicMasses = [
     1.00798,   4.0026,    6.941,     9.012182, 10.811,    12.011,  14.00672,  15.9994,   
    18.9984,   20.18,     22.989767, 24.3051,   26.981538, 28.0855, 30.973762, 32.064,
    35.4527,   39.948,    39.0983,   40.078,    44.95591,  47.88,   50.9415,   51.9961,
    54.93805,  55.844,    58.9332,   58.69,     63.546,    65.4,    69.723,    72.59,     
    74.92159,  78.96,     79.904,    83.8,      85.4678,   87.62,   88.90585,  91.224, 
    92.90638,  95.93,     98.0,     101.07,    102.9055,  106.42,  107.8681,  112.412,
   114.82,    118.71,    121.76,    127.6,     126.90447, 131.29,  132.90544, 137.327,
   138.9054,  140.115,   140.90765, 144.24,    145.0,     150.36,  151.965,   157.25,    
   158.92534, 1612.5,     164.93032, 167.26,    168.93421, 173.04,  174.967,   178.49,
   180.9479,  183.85,    186.207,   190.2,     192.22,    195.08,  196.96654, 200.6, 
   204.3833,  207.2,     208.98037, 209.0,     210.0,     222.0,   223.0,     226.0,
   227.0,     232.0381,  231.04,    238.0289,  237.0,     244.0,   243.0,     247.0,
   247.0,     251.0,     252.0,     257.0,     258.0,     259.0,   262.0 ]
atomicSymbols = [
    'H' ,    'He',    'Li',    'Be',    'B' ,    'C' ,    'N' ,    'O' ,  
    'F' ,    'Ne',    'Na',    'Mg',    'Al',    'Si',    'P' ,    'S' ,  
    'Cl',    'Ar',    'K' ,    'Ca',    'Sc',    'Ti',    'V' ,    'Cr',  
    'Mn',    'Fe',    'Co',    'Ni',    'Cu',    'Zn',    'Ga',    'Ge',  
    'As',    'Se',    'Br',    'Kr',    'Rb',    'Sr',    'Y' ,    'Zr',  
    'Nb',    'Mo',    'Tc',    'Ru',    'Rh',    'Pd',    'Ag',    'Cd',  
    'In',    'Sn',    'Sb',    'Te',    'I' ,    'Xe',    'Cs',    'Ba',  
    'La',    'Ce',    'Pr',    'Nd',    'Pm',    'Sm',    'Eu',    'Gd',  
    'Tb',    'Dy',    'Ho',    'Er',    'Tm',    'Yb',    'Lu',    'Hf',  
    'Ta',    'W' ,    'Re',    'Os',    'Ir',    'Pt',    'Au',    'Hg',  
    'Tl',    'Pb',    'Bi',    'Po',    'At',    'Rn',    'Fr',    'Ra',  
    'Ac',    'Th',    'Pa',    'U' ,    'Np',    'Pu',    'Am',    'Cm',  
    'Bk',    'Cf',    'Es',    'Fm',    'Md',    'No',    'Lr' ]  
assert len(atomicMasses) == len(atomicSymbols), "FEFF Handler: Mismatch in periodic table!"
nElem = len(atomicSymbols)
for i in range(nElem):
    num = i + 1
    sym = atomicSymbols[i]
    mass = atomicMasses[i]
    ptable[num] = {'mass':mass, 'symbol':sym, 'number':num}
    ptable[sym] = {'mass':mass, 'symbol':sym, 'number':num}

def getFeffAtomsFromCluster(input):
    # Find absorbing atom.
    absorber = input['absorbing_atom'][0][0] - 1
    atoms = [x for i,x in enumerate(input['cluster']) if i!=absorber]
    uniqueAtoms = list(set([ x[0] for x in atoms]))
    #print absorber
    #print input['absorbing_atom']
    #print uniqueAtoms
    #print ''.join(''.join([str(e),' ']) for e in input['cluster'][absorber][1:]) + ' 0'
    #print input
    feffAtoms = []
    feffAtoms.append([0.0, 0.0, 0.0, 0])
    for line in atoms: 
        feffAtom = [ e - input['cluster'][absorber][i+1] for i,e in enumerate(line[1:4]) ]
        feffAtom.append(uniqueAtoms.index(line[0])+1)
        #print feffAtom
        #print ''.join(''.join([str(e),' ']) for e in line[1:]) + str(uniqueAtoms.index(line[0])+1)
        #print ptable[x[0]]['number']
        feffAtoms.append(feffAtom)

    return feffAtoms

def getFeffPotentialsFromCluster(input):
    absorber = input['absorbing_atom'][0][0] - 1
    atoms = [x for i,x in enumerate(input['cluster']) if i!=absorber]
    uniqueAtoms = list(set([ x[0] for x in atoms]))
    lfms1 = input.get('feff.lfms1')[0][0]
    lfms2 = input.get('feff.lfms2')[0][0]
    feffPots = [[]]
    feffPots[0] = [0, ptable[input['cluster'][absorber][0]]['number'], input['cluster'][absorber][0], lfms1, lfms2, 1.0 ]
    for i,atm in enumerate(uniqueAtoms):
       xnat = [ x[0] for x in input['cluster'] ].count(atm)
       feffPots.append([i+1, int(ptable[atm]['number']), atm, lfms1, lfms2, xnat ])

    return feffPots

def getFeffDebyeOptions(input):
# Now we have to make the DEBYE input line depending on what input we have
# for the DMDW approach.
# JK - For now options for debye are all defined as just a single string.
#      Maybe split these up in future, since some are optional.
#   pp_debug.pprint(input)
#   sys.exit()
    Debye_Def = None
    # JK - have to have nuctemp in input to run debye waller factors.
    if 'nuctemp' in input:
        Debye_Def  = str(input['nuctemp'][0][0]) + ' '
       
    # If a dynamical matrix is present, can use a made up debye temp if
    # not debyetemp is not present.
    if 'opt_dynmat' in input:
        if 'debyetemp' in input:
            Debye_Def += str(input['debyetemp'][0][0]) + ' '
        else:
            Debye_Def += '300 '

        dymFilename = 'corvus.dym'
        Debye_Def += '5 '
        Debye_Def += dymFilename + ' '
        # To get the right number of Lanczos iterations we should check the number of
        # atoms in the DM, but will leave that for later.
        Debye_Def += str(input['dmdw_nlanczos'][0][0]) + ' 0 1 '

    elif 'debyetemp' in input:
        Debye_Def += str(input['debyetemp'][0][0]) + ' '
        Debye_Def += '0 '
        Debye_Def = [[Debye_Def]]
    else:
        Debye_Def = None
    
#    print 'Debye_Def'
#    print Debye_Def
#   sys.exit()
    return Debye_Def

def cell2atoms(cellatoms, acell, rprim=None, angdeg=None, cutoff=None, nmax=1000):
    # cellatom = dict with 'coord', 'symbol'
    from math import sqrt
    from decimal import Decimal, ROUND_UP
    for at in cellatoms:
        assert ('symbol' in at) and ('coord' in at)
        at['coord'] = [scale * reduced for scale, reduced in zip(acell, at['coord'])]
    atoms = []
    nord = 3 
    iRange = list(range(-nord, nord+1))
    grid = ((i,j,k) for i in iRange for j in iRange for k in iRange)
    for ijk in grid:
        dx = [a*sum(i*c for i,c in zip(ijk,row)) for row, a in zip(rprim,acell)]
        for at in cellatoms:
            copyAt = {k:v for k,v in list(at.items())}
            copyAt['coord'] = list(map(sum, list(zip(at['coord'], dx))))
            copyAt['dist'] = sqrt(sum(x*x for x in copyAt['coord']))
            atoms.append(copyAt) 
    ru8 = lambda x: Decimal(x).quantize(Decimal('0.12345678'),rounding=ROUND_UP)
    atoms.sort(key=lambda at:(ru8(at['dist']), at['symbol']))
    if len(atoms) > nmax:
        return atoms[:nmax]
    else:
        return atoms

def abcell2atoms(input):
    from corvus.abinit import expandedList
    from corvus.conversions import bohr2angstrom
    for key in ['acell','znucl','xred','rprim','natom']:
        assert key in input
        # Expecting input parsed to strings for now
        assert isinstance(input[key], str) 
    # Parsing grammar 
    vector = pp.Group(pp.Word(pp.nums + ".+-E") * 3)
    listVectors = pp.Group(pp.OneOrMore(vector))
    # Translate strings into useable lists/numbers
    acell = bohr2angstrom(list(map(float, expandedList(input['acell']))))
    atnums = list(map(int, expandedList(input['znucl'], length=input['natom'])))
    coords = input['xred']
    rprim = input['rprim']
    try:
        coords = listVectors.parseString(coords).asList()[0]
        coords = [list(map(float, v)) for v in coords]
        rprim = listVectors.parseString(rprim).asList()[0]
        rprim = [list(map(float, v)) for v in rprim]
    except pp.ParseException as pe:
        print(("Parsing Error using pyparsing: invalid input:", pe))
    # Specify a single lmax for each atomic species
    if 'lmax' in input:
        lmax = list(map(int, expandedList(input['lmax'], length=input['natom'])))
    else:
        lmax = list(map(int, expandedList('*-1', length=input['natom'])))
    lmax_by_atnum = {}
    for i,num in enumerate(atnums):
        if num in lmax_by_atnum:
            lmax_by_atnum[num] = max(lmax_by_atnum[num], lmax[i])
        else:
            lmax_by_atnum[num] = lmax[i]
    # Prepare list of cluster atoms, core hole on first atom, as starting point
    # Start with creating cellatoms object to duplicate
    cellatoms = []
    atnumSet = set()
    shift = lambda x: list(map(lambda a,b: a-b, list(zip(x, coords[0]))))
    for i, num in enumerate(atnums):
        at = {'atnum':num, 'symbol':ptable[num]['symbol'], 'coord':shift(coords[i])}
        atnumSet.add(at['atnum'])
        cellatoms.append(at)
    # Core hole atom denoted by atomic symbol, remaining atoms denoted by atomic number
    #   (keys for 'ptable' lookup below)
    atTypeKeys = [cellatoms[0]['symbol']] + sorted(atnumSet)
    # Build list of potentials
    potList = []
    for ipot, key in enumerate(atTypeKeys):
        pot = {'ipot':ipot, 'atnum':ptable[key]['number'], 'tag':ptable[key]['symbol']}
        pot['lmax'] = lmax_by_atnum[pot['atnum']]
        potList.append(pot)
    # Build list of atoms
    atomList = cell2atoms(cellatoms, acell, rprim)
    for at in atomList:
        at['ipot'] = atTypeKeys.index(at['atnum'])
    atomList[0]['ipot'] = 0
    return (potList, atomList)

def dym2atoms(input, center=1):
    from operator import itemgetter
    from math import sqrt
    from corvus.conversions import bohr2angstrom
    assert 'dynmat' in input
    dym = input['dynmat']
    # Check that we have a valid center atom number
    nAt = dym['nAt']
    if center > nAt or center < 1:
        raise ValueError('Index for central atom is invalid')
    iCenter = center - 1 
    # Recenter and order by distance
    shift = lambda x: list(map(lambda a,b: a-b, list(zip(x,dym['atCoords'][iCenter]))))
    swapcoords = [[i,bohr2angstrom(shift(dym['atCoords'][i]))] for i in range(nAt)]
    for a in swapcoords:
        a.append(sqrt(sum(x*x for x in a[1])))
    swapcoords.sort(key=itemgetter(2))
    dym['printOrder'] = [swapcoords[iAt][0] for iAt in range(nAt)]
    # Core hole atom denoted by atomic symbol, remaining atoms denoted by atomic number
    #   (keys for 'ptable' lookup below)
    atTypeKeys = [ptable[dym['atNums'][iCenter]]['symbol']] + sorted(set(dym['atNums']))
    # Build list of potentials
    potList = []
    for ipot, key in enumerate(atTypeKeys):
        pot = {'ipot':ipot, 'atnum':ptable[key]['number'], 'tag':ptable[key]['symbol']}
        potList.append(pot)
    # Build list of atoms
    atomList = []
    for iAt in dym['printOrder']:
        at = {'coord':swapcoords[iAt][1], 'dist':swapcoords[iAt][2]}
        at['symbol'] = ptable[dym['atNums'][iAt]]['symbol']
        at['ipot'] = atTypeKeys.index(dym['atNums'][iAt])
        atomList.append(at)
    atomList[0]['ipot'] = 0
    return (potList, atomList)

def headerLines(input, lines):
    if 'title' in input:
        isStr = lambda x: isinstance(x, str)
        if isStr(input['title']):
            for t in input['title'].split('\n'):
                lines.append('TITLE ' + t)
        elif isinstance(input['title'], list) and all(map(isStr, input['title'])):
            for t in input['title']:
                lines.append('TITLE ' + t)
        lines.append('')

def writeAtomicInput(input, feffinp='feff.inp'):
    
    lines = []
    # Header lines
    lines.append('* This feff9 input file was generated by corvus.')
    
    # Set feff.KEY using setInput, which will use the
    # previously defined value if available, unless Force=True.
    setInput(input,'feff.control',[[1,0,0,0,0,0]],Force=True)
    setInput(input,'feff.print',[[5,0,0,0,0,0]],Force=True)

    writeInput(input,feffinp)


def writeSCFInput(input, feffinp='feff.inp'):

    lines = []
    lines.append('* This feff9 input file was generated by corvus.')
    lines.append('')

    setInput(input,'feff.control', [[1, 0, 0, 0, 0, 0]])
    setInput(input,'feff.print', [[5, 0, 0, 0, 0, 0]])
    setInput(input,'feff.scf',[[5.0, 0, 100, 0.1, 0]])
    setInput(input,'feff.edge',[['K']])
    
    writeInput(input,feffinp)

def writeCrossSectionsInput(input, feffinp='feff.inp'):
    lines = []
    lines.append('* This feff9 input file was generated by corvus.')
    lines.append('')

    setInput(input,'feff.control', [[1, 1, 0, 0, 0, 0]])
    setInput(input,'feff.print', [[5, 2, 0, 0, 0, 0]])
    setInput(input,'feff.exchange', [[0, 0.0, 0.0, 0]])
    setInput(input,'feff.edge',[['K']])

    writeInput(input,feffinp)


def writeGreensFunctionInput(input, feffinp='feff.inp'):
    lines = []
    lines.append('* This feff9 input file was generated by corvus.')
    lines.append('')

    setInput(input,'feff.control', [[1, 1, 1, 0, 0, 0]])
    
    writeInput(input,feffinp)


def writePathsInput(input, feffinp='feff.inp'):
    lines = [] 
    lines.append('* This feff9 input file was generated by corvus.')
    lines.append('')

    setInput(input,'feff.control', [[1, 1, 1, 1, 0, 0]])

    writInput(input,feffinp)


def writeFMatricesInput(input, feffinp='feff.inp'):
    lines = []
    lines.append('* This feff9 input file was generated by corvus.')
    lines.append('')

    setInput(input,'feff.control', [[1, 1, 1, 1, 1, 0]])


    writeInput(input, feffinp)


def writeXANESInput(input, feffinp='feff.inp'):
    lines = []
    lines.append('* This feff9 input file was generated by corvus.')
    lines.append('')

    setInput(input,'feff.control', [[1, 1, 1, 1, 1, 1]])
    setInput(input,'feff.print', [[0, 0, 0, 0, 0, 0]])
    setInput(input,'feff.exchange', [[0, 0.0, 0.0, 0]])
    setInput(input,'feff.xanes', [[10]])
    setInput(input,'feff.fms', [[6.0, 0]])
    setInput(input,'feff.scf', [[5.0, 0]])
    setInput(input,'feff.rpath', [[0.1]])
    setInput(input,'feff.edge', [['K']])

    writeInput(input,feffinp)


def writeXESInput(input, feffinp='feff.inp'):
    lines = []
    lines.append('* This feff9 input file was generated by corvus.')
    lines.append('')

    setInput(input,'feff.control', [[1, 1, 1, 1, 1, 1]])
    setInput(input,'feff.print', [[0, 0, 0, 0, 0, 0]])
    setInput(input,'feff.exchange', [[0, 0.0, 0.0, 0]])
    setInput(input,'feff.xes', [[-30, 0]])
    setInput(input,'feff.fms', [[6.0, 0]])
    setInput(input,'feff.scf', [[5.0, 0]])
    setInput(input,'feff.rpath', [[0.1]])
    setInput(input,'feff.edge', [['K']])

    writeInput(input,feffinp)


def writeELNESinp(input, feffinp='feff.inp'):
    for key in ['potlist', 'atomlist']:
        assert key in input
    lines = []
    lines.append(' * This feff9 input file was generated by corvus.feff.writeELNESinp')
    lines.append('')
    # Header part of feff.inp
    headerLines(input, lines)
    # Control and print cards. Default to run everything and no extra output
    lines.append(' *              pot    xsph     fms   paths  genfmt  ff2chi')
    lines.append(' CONTROL   ' + check(input, 'control', default=('{:8d}'*6).format(*[1]*6)))
    lines.append(' PRINT     ' + check(input, 'print', default=('{:8d}'*6).format(*[0]*6)))
    lines.append('')
    # Default to Cu example
    lines.append(' EDGE      ' + check(input, 'edge', default='K'))
    if 'corehole' in input:
        lines.append(' COREHOLE  ' + input['corehole'])
    lines.append(' SCF       ' + str(check(input, 'scf', default=4.0)))
    lines.append(' FMS       ' + str(check(input, 'fms', default=6.0)))
    lines.append('')
    defaultELNES = '\n'.join(['20.0 0.07 0.0','300','0 1 0','2.4 0.0','5 3','0.0 0.0'])
    lines.append(' ELNES     ' + check(input, 'elnes', default=defaultELNES))
    lines.append('')
    if 'egrid' in input:
        lines.append(' EGRID ')
        lines.append(input['egrid'])
    lines.append('')
    # List potentials and atoms
    lines.append('POTENTIALS')
    potstr = '{:5d}{:5d}   {:2}  {:3d} {:3d}'
    potstrShort = '{:5d}{:5d}   {:2}'
    for p in input['potlist']:
        if 'lmax' in p:
            lines.append(potstr.format(p['ipot'], p['atnum'], p['tag'], p['lmax'], p['lmax']))
        else:
            lines.append(potstrShort.format(p['ipot'], p['atnum'], p['tag']))
    lines.append('')
    # Input for the structure part of the feff.inp file    
    lines.append('ATOMS')
    atomstr = '{:11.5f}'*3 + '{:5d}   {:2}{:9.5f}{:5d}'
    for iAt, atom in enumerate(input['atomlist']):
        atLine = atom['coord'] + [atom['ipot'], atom['symbol'], atom['dist'], iAt]
        lines.append(atomstr.format(*atLine))
    lines.append('')
    lines.append('END')
    # Print feff input file
    writeList(lines, feffinp)

# Added by FDV
def writeEXAFSDMDWInput(input, feffinp='feff.inp'):

    lines = []
    lines.append('* This feff9 input file was generated by corvus.')
    lines.append('')


    setInput(input,'feff.control', [[1, 1, 1, 1, 1, 1]])
    setInput(input,'feff.print', [[0, 0, 0, 0, 0, 0]])
    setInput(input,'feff.exchange', [[0, 0.0, 0.0, 0]])
    setInput(input,'feff.exafs', [[18]])
    setInput(input,'feff.scf', [[4.0, 0]])
    setInput(input,'feff.rpath', [[8.0]])
    setInput(input,'feff.nleg', [[4]])
    setInput(input,'feff.edge', [['K']])

# JK -  Commented below: I'm making this a function getFeffDebyeInput that translates generic 
# input to feff input, which will be performed at the begining of the run 
# function.
# Now we have to make the DEBYE input line depending on what input we have
# for the DMDW approach.
#    dymFilename = 'corvus.dym'
#   pp_debug.pprint(input)
#   sys.exit()
#    Debye_Def  = str(input['nuctemp'][0][0]) + ' '
#    Debye_Def += str(input['debyetemp'][0][0]) + ' '
#    Debye_Def += '5 '
#    Debye_Def += dymFilename + ' '

# To get the right number of Lanczos iterations we should check the number of
# atoms in the DM, but will leave that for later.
#    Debye_Def += str(input['dmdw_nlanczos'][0][0]) + ' 0 1 '

#    print 'Debye_Def'
#    print Debye_Def
#   sys.exit()

    writeInput(input,feffinp)

def legacy_dym2feffinp(dym, center=1, feffinp='feff.inp', feffdym='feff.dym', spectrum='EXAFS', header=True, input={}):
    from operator import itemgetter
    from math import sqrt
    from corvus.dmdw import writeDym

    # Check that we have a valid center atom number
    nAt = dym['nAt']
    if center > nAt or center < 1:
        raise ValueError('Index for central atom is invalid')
    iCenter = center - 1 

    # Recenter and order by distance
    shift = lambda x: list(map(lambda a,b: a-b, list(zip(x,dym['atCoords'][iCenter]))))
    swapcoords = [[i,shift(dym['atCoords'][i])] for i in range(nAt)]
    for a in swapcoords:
        a.append(sqrt(sum(x*x for x in a[1])))
    swapcoords.sort(key=itemgetter(2))
    dym['printOrder'] = [swapcoords[j][0] for j in range(nAt)]
    
    # Atomic types as list of keys for periodic table
    #   center atom indicated by symbol instead of number for bijective map between
    #   index and key
    atTypeKeys = [ptable[dym['atNums'][iCenter]]['symbol']] + sorted(set(dym['atNums']))
    atTypeList = [atTypeKeys.index(dym['atNums'][i]) for i in dym['printOrder']]
    atTypeList[0] = 0

    # Gather input for feff.inp file
    lines = []
    # Header part of feff.inp
    if header:
        lines.append(' * This feff9 input file was generated by corvus.feff.dym2feffinp')
        lines.append('')
        lines.append(' TITLE dymfile name:  ' + feffdym)
        lines.append(' TITLE absorbing atom:' + '{:4d}'.format(0))
        lines.append('')
    # Edge and s0^2 information. Default to K and 1.0 respectively for now
    lines.append(' EDGE      ' + check(input, 'edge', default='K'))
    lines.append(' S02       ' + '{:6.4f}'.format(check(input, 's02', default=1.0)))
    lines.append('')
    # Control and print cards. Default to run everything and extra output for pot module
    lines.append(' *              pot    xsph     fms   paths  genfmt  ff2chi')
    lines.append(' CONTROL   ' + check(input, 'control', default=('{:8d}'*6).format(*[1]*6)))
    lines.append(' PRINT     ' + check(input, 'print', default=('{:8d}'*6).format(*[1]+[0]*5)))
    lines.append('')
    # Exchange controls. Default to 0 (Hedin-Lundqvist)
    lines.append(' *          ixc  [ Vr  Vi ]')
    lines.append(' EXCHANGE  ' + '{:4d}'.format(check(input,'ixc',default=0)))
    lines.append('')
    # SCF radius. Default to 4.0
    lines.append(' *            r_scf  [ l_scf   n_scf   ca ]')
    lines.append(' SCF       ' + '{:8.3f}'.format(check(input, 'scf', default=4.0)))
    lines.append('')
    # XANES info if applicable
    if spectrum == 'XANES':
        # XANES upper limit. Default to 4.0
        lines.append(' *             kmax   [ delta_k  delta_e ]')
        lines.append(' XANES     ' + '{:8.3f}'.format(check(input, 'kmax', default=8.0)))
        lines.append('')
        # FMS radius. Defaults to 6.0
        lines.append(' *            r_fms  [l_fms]')
        lines.append(' FMS       ' + '{:8.3f}'.format(check(input, 'r_fms', default=6.0)))
        lines.append(' RPATH     ' + '{:8.3f}'.format(check(input, 'rpath', default=0.1)))
        lines.append('')
    # EXAFS information if applicable
    if spectrum == 'EXAFS':
        lines.append(' EXAFS     ' + '{:8.3f}'.format(check(input, 'kmax', default=20.0)))
        lines.append(' RPATH     ' + '{:8.3f}'.format(check(input, 'rpath', default=8.0)))
        lines.append('')
    # DEBYE information
    lines.append(' *        Temp  Debye_Temp  DW_Opt' + ' '*9 + 'dymfile  DMDW_Order  DMDW_Type  DMDW_Route')
    defaultDebye = '{:6.1f}'.format(check(input, 'temp', default=450.0))
    defaultDebye += ' '*6 + '{:6.1f}'.format(check(input, 'debyetemp', default=315.0))
    defaultDebye += ' '*7 + '{:1d}'.format(check(input, 'dw_opt', default=5))
    defaultDebye += ' '*1 + '{:>15}'.format(feffdym)
    defaultDebye += ' '*10 + '{:2d}'.format(check(input, 'dmdw_order', default=6))
    defaultDebye += ' '*9 + '{:2d}'.format(check(input, 'dmdw_type', default=0))
    defaultDebye += ' '*10 + '{:2d}'.format(check(input, 'dmdw_route', default=1))
    lines.append(' DEBYE  ' + check(input, 'debye', default=defaultDebye)) 
    lines.append('')
    # Input for the potentials part of the feff.inp file
    lines.append('POTENTIALS')
    potstr = '{:5d}{:5d}   {:2}'
    for ipot, key in enumerate(atTypeKeys):
        lines.append(potstr.format(ipot, ptable[key]['number'], ptable[key]['symbol']))
    lines.append('')
    # Input for the structure part of the feff.inp file    
    lines.append('ATOMS')
    atomstr = '{:11.5f}'*3 + '{:5d}   {:2}{:9.5f}{:5d}'
    for iAt, atType in enumerate(atTypeList):
        atLine = swapcoords[iAt][1] + [atType, ptable[atTypeKeys[atType]]['symbol']]
        atLine += [swapcoords[iAt][2], iAt]
        lines.append(atomstr.format(*atLine))    
    lines.append('')
    lines.append('END')
    # Print feff input file
    writeList(lines, feffinp)
    # Print dym file
    writeDym(dym, feffdym)


feff_edge_dict={
 'H':
['K'],
 'He':
['K'],
 'Li':
['K', 'L1'],
 'Be':
['K', 'L1'],
 'B':
['K', 'L1', 'L2'],
 'C':
['K', 'L1', 'L2', 'L3'],
 'N':
['K', 'L1', 'L2', 'L3'],
 'O':
['K', 'L1', 'L2', 'L3'],
 'F':
['K', 'L1', 'L2', 'L3'],
 'Ne':
['K', 'L1', 'L2', 'L3'],
 'Na':
['K', 'L1', 'L2', 'L3', 'M1'],
 'Mg':
['K', 'L1', 'L2', 'L3', 'M1', 'M2'],
 'Al':
['K', 'L1', 'L2', 'L3', 'M1', 'M2'],
 'Si':
['K', 'L1', 'L2', 'L3', 'M1', 'M2'],
 'P':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3'],
 'S':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3'],
 'Cl':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3'],
 'Ar':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3'],
 'K':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'N1'],
 'Ca':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'N1', 'N2'],
 'Sc':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'N1', 'N2'],
 'Ti':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'N1', 'N2'],
 'V':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'N1', 'N2'],
 'Cr':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'N1', 'N2'],
 'Mn':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2'],
 'Fe':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2'],
 'Co':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2'],
 'Ni':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2'],
 'Cu':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1'],
 'Zn':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2'],
 'Ga':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2'],
 'Ge':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2'],
 'As':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3'],
 'Se':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3'],
 'Br':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3'],
 'Kr':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3'],
 'Rb':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'O1'],
 'Sr':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'O1'],
 'Y':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'O1'],
 'Zr':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'O1'],
 'Nb':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'O1'],
 'Mo':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'O1'],
 'Tc':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'O1'],
 'Ru':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'O1'],
 'Rh':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'O1'],
 'Pd':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5'],
 'Ag':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'O1'],
 'Cd':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'O1'],
 'In':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'O1', 'O2'],
 'Sn':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'O1', 'O2'],
 'Sb':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'O1', 'O2', 'O3'],
 'Te':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'O1', 'O2', 'O3'],
 'I':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'O1', 'O2', 'O3'],
 'Xe':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'O1', 'O2', 'O3'],
 'Cs':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'O1', 'O2', 'O3', 'P1'],
 'Ba':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'O1', 'O2', 'O3', 'P1'],
 'La':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'O1', 'O2', 'O3', 'O4', 'P1'],
 'Ce':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'O1', 'O2', 'O3', 'O4', 'P1'],
 'Pr':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'O1', 'O2', 'O3', 'O4', 'P1'],
 'Nd':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'O1', 'O2', 'O3', 'O4', 'P1'],
 'Pm':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'O1', 'O2', 'O3', 'O4', 'P1'],
 'Sm':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'O1', 'O2', 'O3', 'O4', 'P1'],
 'Eu':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'O1', 'O2', 'O3', 'O4', 'P1'],
 'Gd':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'P1'],
 'Tb':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'P1'],
 'Dy':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'P1'],
 'Ho':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'P1'],
 'Er':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'P1'],
 'Tm':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'P1'],
 'Yb':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'P1'],
 'Lu':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'P1'],
 'Hf':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'P1'],
 'Ta':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'P1'],
 'W':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'P1'],
 'Re':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'P1'],
 'Os':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'P1'],
 'Ir':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'P1'],
 'Pt':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'P1'],
 'Au':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'P1'],
 'Hg':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'P1'],
 'Tl':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'P1', 'P2'],
 'Pb':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'P1', 'P2'],
 'Bi':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'P1', 'P2', 'P3'],
 'Po':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'P1', 'P2', 'P3'],
 'At':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'P1', 'P2', 'P3'],
 'Rn':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'P1', 'P2', 'P3'],
 'Fr':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'P1', 'P2', 'P3', 'R1'],
 'Ra':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'P1', 'P2', 'P3', 'R1'],
 'Ac':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'P1', 'P2', 'P3', 'P4', 'R1'],
 'Th':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'P1', 'P2', 'P3', 'P4', 'R1'],
 'Pa':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'P1', 'P2', 'P3', 'P4', 'R1'],
 'U':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'P1', 'P2', 'P3', 'P4', 'R1'],
 'Np':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'P1', 'P2', 'P3', 'P4', 'R1'],
 'Pu':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'P1', 'P2', 'P3', 'R1'],
 'Am':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'R1'],
 'Cm':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'R1'],
 'Bk':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'R1'],
 'Cf':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'R1'],
 'Es':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'R1'],
 'Fm':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'R1'],
 'Md':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'R1'],
 'No':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'R1'],
 'Lr':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'R1', 'R2'],
 'Rf':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'R1'],
 'Db':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'R1'],
 'Sg':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'R1'],
 'Bh':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'P5', 'R1'],
 'Hs':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'P5', 'R1'],
 'Mt':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'P5', 'R1'],
 'Ds':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'P5', 'R1'],
 'Rg':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'P5', 'R1'],
 'Cn':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'P5', 'R1'],
 'Uut':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'P5', 'R1', 'R2'],
 'Fl':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'P5', 'R1', 'R2'],
 'UUp':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'P5', 'R1', 'R2', 'R3'],
 'Lv':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'P5', 'R1', 'R2', 'R3'],
 'Uus':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'P5', 'R1', 'R2', 'R3'],
 'Uuo':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'P5', 'R1', 'R2', 'R3'],
 'Uue':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'P5', 'R1', 'R2', 'R3', 'O8'],
 'Ubn':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'P5', 'R1', 'R2', 'R3', 'S1'],
 'Ubu':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'P5', 'R1', 'R2', 'R3', 'S1', 'S2'],
 'Ubb':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'P5', 'R1', 'R2', 'R3', 'S1', 'S2', 'R4'],
 'Ubt':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'P5', 'R1', 'R2', 'R3', 'S1', 'S2', 'R4', 'P6'],
 'Ubq':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'P5', 'R1', 'R2', 'R3', 'S1', 'S2', 'P6'],
 'Ubp':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'P5', 'R1', 'R2', 'R3', 'S1', 'S2', 'P6', 'O8'],
 'Ubh':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'P5', 'R1', 'R2', 'R3', 'S1', 'S2', 'P6', 'O8'],
 'Ubs':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'P5', 'R1', 'R2', 'R3', 'S1', 'S2', 'P6', 'O8'],
 'Ubo':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'P5', 'R1', 'R2', 'R3', 'S1', 'S2', 'P6', 'O8'],
 'Ube':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'P5', 'R1', 'R2', 'R3', 'S1', 'S2', 'P6', 'O8'],
 'Utn':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'P5', 'R1', 'R2', 'R3', 'S1', 'S2', 'P6', 'O8'],
 'Utu':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'P5', 'R1', 'R2', 'R3', 'S1', 'S2', 'P6', 'O8'],
 'Utb':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'P5', 'R1', 'R2', 'R3', 'S1', 'S2', 'P6', 'O8'],
 'Utt':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'P5', 'R1', 'R2', 'R3', 'S1', 'S2', 'P6', 'O8'],
 'Utq':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'P5', 'R1', 'R2', 'R3', 'S1', 'S2', 'P6', 'O8'],
 'Utp':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'P5', 'R1', 'R2', 'R3', 'S1', 'S2', 'P6', 'O8', 'O9'],
 'Uth':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'P5', 'R1', 'R2', 'R3', 'S1', 'S2', 'P6', 'O8', 'O9'],
 'Uts':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'P5', 'R1', 'R2', 'R3', 'S1', 'S2', 'P6', 'O8', 'O9'],
 'Uto':
['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'P5', 'R1', 'R2', 'R3', 'S1', 'S2', 'P6', 'O8', 'O9'],
 }


def kk_transform(w,eps2):
  
    """ Performs kk-transform on imaginary part of spectrum from 0 to infty
    input:
        w  - energy grid, does not need to be even
        eps2 - imaginary part of some spectrum

    output:
        f1 - real part of kk-tranform (will need to add constant if necessary
    """

    # Create new frequency grid that is centered between the points in w.
    wnew = (w[:-1] + w[1:])/2.0
    eps1 = np.zeros_like(wnew)
    for iw, w0 in enumerate(wnew):
        if False:
            for jw,w1 in enumerate(w[:-1]):
                integrand = w1*eps2[jw]/(w1**2-w0**2)
                delta = w[jw+1]-w1
                #eps1[iw] = np.trapz(integrand[:iw],w[:iw]) + np.trapz(integrand[iw+1:],w[iw+1:])
                a = (eps2[jw+1] - eps2[jw])/delta
                b = eps2[jw] - a*w[jw]
                g1 = (w0 + w1)/(w0 - w1)*(w0 - w1 - delta)/(w0 + w1 + delta)
                g1 = np.log(abs(g1))
                g2 = ((w1 + delta)**2 - w0**2)/(w1**2 - w0**2)
                g2 = np.log(abs(g2))
                eps1[iw] = eps1[iw] + a*delta + a*w0/2.0*g1 + b/2.0*g2


        integrand = w[:-1]*eps2[:-1]/(w[:-1]**2 - w0**2)
        delta = w[1:] - w[:-1]
        a = (eps2[1:] - eps2[:-1])/delta
        b = eps2[:-1] - a*w[:-1]
        g1 = (w0 + w[:-1])/(w0 - w[:-1])*(w0 - w[:-1] - delta)/(w0 + w[:-1] + delta)
        g1 = np.log(np.abs(g1))
        g2 = ((w[:-1] + delta)**2 - w0**2)/(w[:-1]**2 - w0**2)
        g2 = np.log(np.abs(g2))
        eps1[iw] = np.sum(a*delta + a*w0/2.0*g1 + b/2.0*g2)
    return wnew,eps1*2/np.pi

def dos_conv(e1,EFermi,E0,k,xanes,w_in, dos_in):
    ''' Convolve xmu.dat with appropriate LDOS to get low frequency optical spectrum. '''

    w = w_in[:]
    dos = dos_in[:]
    # Get all maxima in the dos
    ind_max = argrelextrema(dos,np.greater)[0]

    Etop = EFermi
    EGap = 0.0
    itop = -1
    # Now find maxima that are closest to Fermi level
    for i,ind in enumerate(ind_max[:-1]):
      if w[ind] < EFermi < w[ind_max[i+1]]:
        # Now run backward to find first point where dos < half max
        iw = ind_max[i+1]
        while iw > ind:
          if dos[iw] < dos[ind_max[i+1]]/2.0:
            Etop = w[iw]
            itop = iw
            break
          iw -= 1

        if itop < 0:
          break
        # Now run forward to find same from bottom of gap
        iw = ind
        while iw < ind_max[i+1]:
          if dos[iw] < dos[ind]/2.0:
            EGap = Etop - w[iw]
            ibottom = iw
            break
          iw += 1

        break
      
    print('Gap energy:', EGap)
    print('EFermi:', EFermi)
    print('ETop:', Etop)
    #import matplotlib.pyplot as plt
    #plt.plot(w-EFermi,dos)
    #plt.show()
    
    # Redefine DOS as occupied DOS.
    dos = dos[w<EFermi]
    w = w[w<EFermi]
    dos = dos/np.trapz(dos,w)

    # Make grids for dos (-100 to 100)
    e_step = 0.1
    e_grid = np.arange(0.1,100,0.1) 
    e_grid2 = np.flip(-e_grid)
    dos_terp = np.interp(e_grid2,w,dos,left=0.0,right=0.0)
    dos_terp = dos_terp/np.trapz(dos_terp)/0.1
    for i,en in enumerate(e_grid2):
      if en > EFermi:
        dos_terp[i] = 0.0
      mu_terp = np.interp(e1+E0-(Etop-en),e1,xanes*e1,left=0.0,right=xanes[-1])
      # Trapezoidal rule integration for mu. 
      if i == 0:
        mu = mu_terp*dos_terp[i]/e1
      else:
        mu = mu + mu_terp*dos_terp[i]*2.0/e1

    mu = mu*e_step/2.0

    # Cut mu off below EGap.
    mu[np.where(e1<EGap)] = 0.0
    return mu

def getHoleSymm(ihole):
    # K, L1, L2, L3, M1, M2, M3, M4, M5
    #['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'P1', 'P2', 'P3', 'P4', 'P5', 'R1', 'R2', 'R3', 'S1', 'S2', 'P6', 'O8', 'O9'],
    holeSymm = [0,     0,    1,    1,    0,    1,    1,    2,    2,    0,    1,    1,    2,    2,    3,    3,    0,    1,    1,    2,    2,    3,    3,    0,    1,    1,    2,    2,    0,    1,    1,    0,    1,    0,    4,    4]
    return holeSymm[ihole]
