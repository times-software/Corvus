from structures import Handler, Exchange, Loop, Update
import corvutils.pyparsing as pp
import os, sys, subprocess, shutil, resource
import re

# Debug: FDV
import pprint
pp_debug = pprint.PrettyPrinter(indent=4)

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
            raise LookupError('Corvus cannot currently produce ' + key + ' using FEFF')
        return implemented[key]['cost']

    @staticmethod
    def sequenceFor(output,inp=None):
        if isinstance(output, list) and output and isinstance(output[0], basestring):
            key = strlistkey(output)
        elif isinstance(output, basestring):
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
        subdir = config['pathprefix'] + str(config['xcIndex']) + '_FEFF'
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
        import numpy as np
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
        dir = config['xcDir']

        # Set input file
        inpf = os.path.join(dir, 'feff.inp')
        # Loop over targets in output. Not sure if there will ever be more than one output target here.
        for target in output:
            if (target == 'feffAtomicData'):
                # Set output and error files
                with open(os.path.join(dir, 'corvus.FEFF.stdout'), 'w') as out, open(os.path.join(dir, 'corvus.FEFF.stderr'), 'w') as err:

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

                        runExecutable('',dir,executable,args,out,err)

                # For this case, I am only passing the directory for now so
                # that other executables in FEFF can use the atomic data.

                output[target] = dir

            elif (target == 'feffSCFPotentials'):
                # Set output and error files
                with open(os.path.join(dir, 'corvus.FEFF.stdout'), 'w') as out, open(os.path.join(dir, 'corvus.FEFF.stderr'), 'w') as err:

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

                        runExecutable('',dir,executable,args,out,err)

                # For this case, I am only passing the directory for now so
                # that other executables in FEFF can use the atomic data.

                output[target] = dir
            elif (target == 'feffCrossSectionsAndPhases'):
                # Set output and error files
                with open(os.path.join(dir, 'corvus.FEFF.stdout'), 'w') as out, open(os.path.join(dir, 'corvus.FEFF.stderr'), 'w') as err:

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

                        runExecutable('',dir,executable,args,out,err)
                output[target] = dir

 
            elif (target == 'feffGreensFunction'):

                # Set output and error files
                with open(os.path.join(dir, 'corvus.FEFF.stdout'), 'w') as out, open(os.path.join(dir, 'corvus.FEFF.stderr'), 'w') as err:


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

                        runExecutable('',dir,executable,args,out,err)

                # For this case, I am only passing the directory for now so
                # that other executables in FEFF can use the atomic data.

                output[target] = dir


            elif (target == 'feffPaths'):
                # Set output and error files
                with open(os.path.join(dir, 'corvus.FEFF.stdout'), 'w') as out, open(os.path.join(dir, 'corvus.FEFF.stderr'), 'w') as err:

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

                        runExecutable('',dir,executable,args,out,err)

                # For this case, I am only passing the directory for now so
                # that other executables in FEFF can use the atomic data.

                output[target] = dir

            elif (target == 'feffFMatrices'):
                # Set output and error files
                with open(os.path.join(dir, 'corvus.FEFF.stdout'), 'w') as out, open(os.path.join(dir, 'corvus.FEFF.stderr'), 'w') as err:

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

                        runExecutable('',dir,executable,args,out,err)
                        
                output[target] = dir

            elif (target == 'feffXANES'):
                # Set output and error files
                with open(os.path.join(dir, 'corvus.FEFF.stdout'), 'w') as out, open(os.path.join(dir, 'corvus.FEFF.stderr'), 'w') as err:

                    # Write input file for FEFF.
                    writeXANESInput(feffInput,inpf)

                    # Loop over executable: This is specific to feff. Other codes
                    # will more likely have only one executable. Here I am running 
                    # rdinp again since writeSCFInput may have different cards than

                    execs = ['rdinp','atomic','pot','screen','opconsat','xsph','fms','mkgtr','path','genfmt','ff2x','sfconv']
                    for exe in execs:
                        if 'feff.MPI.CMD' in feffInput:
                            executable = feffInput.get('feff.MPI.CMD')[0]
                            args = feffInput.get('feff.MPI.ARGS',[['']])[0] + [os.path.join(feffdir,exe)]
                        else:
                            executable = [os.path.join(feffdir,exe)]
                            args = ['']

                        runExecutable('',dir,executable,args,out,err)

                outFile=os.path.join(dir,'xmu.dat')
                output[target] = np.loadtxt(outFile,usecols = (0,3)).T.tolist()
                #print output[target]


            elif (target == 'feffXES'):
                # Set output and error files
                with open(os.path.join(dir, 'corvus.FEFF.stdout'), 'w') as out, open(os.path.join(dir, 'corvus.FEFF.stderr'), 'w') as err:

                    # Write input file for FEFF.
                    writeXESInput(feffInput,inpf)

                    # Loop over executable: This is specific to feff. Other codes
                    # will more likely have only one executable. 
                    execs = ['rdinp','atomic','pot','screen','opconsat','xsph','fms','mkgtr','path','genfmt','ff2x','sfconv']
                    for exe in execs:
                        if 'feff.MPI.CMD' in feffInput:
                            executable = feffInput.get('feff.MPI.CMD')[0]
                            args = feffInput.get('feff.MPI.ARGS',[['']])[0] + [os.path.join(feffdir,exe)]
                        else:
                            executable = [os.path.join(feffdir,exe)]
                            args = ['']

                        runExecutable('',dir,executable,args,out,err)


                # For this case, I am only passing the directory for now so
                # that other executables in FEFF can use the data.

                outFile=os.path.join(dir,'xmu.dat')
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
                    dirname = os.path.join(dir,edge)
                    if not os.path.exists(dirname):
                        os.mkdir(dirname)

                    if edge.upper() != "VAL":

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
                        xesdir = os.path.join(dir,'XES')
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
                    shutil.copyfile(os.path.join(dirname,'wscrn.dat'), os.path.join(dir,'wscrn_' + str(nEdge) + '.dat'))
                    shutil.copyfile(os.path.join(dirname,'phase.bin'), os.path.join(dir,'phase_' + str(nEdge) + '.bin'))
                    shutil.copyfile(os.path.join(dirname,'gg.bin'), os.path.join(dir,'gg_' + str(nEdge) + '.bin'))
                    shutil.copyfile(os.path.join(dirname,'xsect.dat'), os.path.join(dir,'xsect_' + str(nEdge) + '.dat'))
                    shutil.copyfile(os.path.join(dirname,'xmu.dat'), os.path.join(dir,'xmu_' + str(nEdge) + '.dat'))
                    shutil.copyfile(os.path.join(dirname,'rl.dat'), os.path.join(dir,'rl_' + str(nEdge) + '.dat'))
                    shutil.copyfile(os.path.join(dirname,'.dimensions.dat'), os.path.join(dir,'.dimensions.dat'))
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
                
                feffinp = os.path.join(dir, 'feff.inp')
                # Write XANES input for this run
                writeXANESInput(feffInput,feffinp)

                # Set output and error files                
                with open(os.path.join(dir, 'corvus.FEFF.stdout'), 'w') as out, open(os.path.join(dir, 'corvus.FEFF.stderr'), 'w') as err:
                    execs = ['rdinp','atomic','rixs']
                    for exe in execs:
                        if 'feff.MPI.CMD' in feffInput:
                            executable = feffInput.get('feff.MPI.CMD')[0]
                            args = feffInput.get('feff.MPI.ARGS',[['']])[0] + [os.path.join(feffdir,exe)]
                        else:
                            executable = [os.path.join(feffdir,exe)]
                            args = ['']

                        runExecutable('',dir,executable,args,out,err)
                        
                outFile=os.path.join(dir,'rixsET.dat')
                output[target] = np.loadtxt(outFile).T.tolist()

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
                with open(os.path.join(dir, 'corvus.FEFF.stdout'), 'w') as out, open(os.path.join(dir, 'corvus.FEFF.stderr'), 'w') as err:

# Write input file for FEFF.
                  writeEXAFSDMDWInput(feffInput,inpf)

# Before running, we need to write the dym file to be used in this run
                  dymFilename = 'corvus.dym'
                  from dmdw import writeDym
                  writeDym(input['opt_dynmat'], os.path.join(dir, dymFilename))
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
                    runExecutable('',dir,executable,args,out,err)

                output[target] = dir

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
        if key != 'feff.end' and not key.startswith('feff.MPI'):
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
    print('Running exectuable: ' + executable[0] + ' ' + ' '.join(args))
    # Modified by FDV:
    # Adding the / to make the config more generic
    # Modified by JJK to use os.path.join (even safer than above).
    execList = [os.path.join(execDir,executable[0])] + args
    p = subprocess.Popen(execList, bufsize=0, cwd=workDir, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    while True:
        pout = p.stdout.readline()
        if pout == '' and p.poll() is not None:
            break
        if pout:
            print(pout.strip())
            out.write(pout)

    while True:
        perr = p.stderr.readline()
        if perr == '' and p.poll() is not None:
            break
        if perr:
            print('###################################################')
            print('###################################################')
            print('Error in executable: ' + executable[0])
            print(perr.strip())
            print('###################################################')
            print('###################################################')
            err.write(perr)


    p.wait()
    
    
def readColumns(filename, columns=[1,2]):
    # Read file and clear out comments
    with open(filename, 'r') as file:
        cleanStr = file.read()
    comments = pp.ZeroOrMore(pp.pythonStyleComment).setParseAction(pp.replaceWith(''))
    try:
        cleanStr = comments.transformString(cleanStr)
    except pp.ParseException as pe:
        print('Parsing Error using pyparsing: invalid input:', pe)
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
        print('Parsing Error using pyparsing: invalid input:', pe)
        sys.exit()
    cols = map(list, zip(*data))
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
for i in xrange(nElem):
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
    feffPots = [[]]
    feffPots[0] = [0, ptable[input['cluster'][absorber][0]]['number'], input['cluster'][absorber][0], -1, -1, 1.0 ]
    for i,atm in enumerate(uniqueAtoms):
       xnat = [ x[0] for x in input['cluster'] ].count(atm)
       feffPots.append([i+1, int(ptable[atm]['number']), atm, -1, -1, xnat ])

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
    iRange = range(-nord, nord+1)
    grid = ((i,j,k) for i in iRange for j in iRange for k in iRange)
    for ijk in grid:
        dx = [a*sum(i*c for i,c in zip(ijk,row)) for row, a in zip(rprim,acell)]
        for at in cellatoms:
            copyAt = {k:v for k,v in at.items()}
            copyAt['coord'] = map(sum, zip(at['coord'], dx))
            copyAt['dist'] = sqrt(sum(x*x for x in copyAt['coord']))
            atoms.append(copyAt) 
    ru8 = lambda x: Decimal(x).quantize(Decimal('0.12345678'),rounding=ROUND_UP)
    atoms.sort(key=lambda at:(ru8(at['dist']), at['symbol']))
    if len(atoms) > nmax:
        return atoms[:nmax]
    else:
        return atoms

def abcell2atoms(input):
    from abinit import expandedList
    from conversions import bohr2angstrom
    for key in ['acell','znucl','xred','rprim','natom']:
        assert key in input
        # Expecting input parsed to strings for now
        assert isinstance(input[key], basestring) 
    # Parsing grammar 
    vector = pp.Group(pp.Word(pp.nums + ".+-E") * 3)
    listVectors = pp.Group(pp.OneOrMore(vector))
    # Translate strings into useable lists/numbers
    acell = bohr2angstrom(map(float, expandedList(input['acell'])))
    atnums = map(int, expandedList(input['znucl'], length=input['natom']))
    coords = input['xred']
    rprim = input['rprim']
    try:
        coords = listVectors.parseString(coords).asList()[0]
        coords = [map(float, v) for v in coords]
        rprim = listVectors.parseString(rprim).asList()[0]
        rprim = [map(float, v) for v in rprim]
    except pp.ParseException as pe:
        print("Parsing Error using pyparsing: invalid input:", pe)
    # Specify a single lmax for each atomic species
    if 'lmax' in input:
        lmax = map(int, expandedList(input['lmax'], length=input['natom']))
    else:
        lmax = map(int, expandedList('*-1', length=input['natom']))
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
    shift = lambda x: map(lambda (a,b): a-b, zip(x, coords[0]))
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
    from conversions import bohr2angstrom
    assert 'dynmat' in input
    dym = input['dynmat']
    # Check that we have a valid center atom number
    nAt = dym['nAt']
    if center > nAt or center < 1:
        raise ValueError('Index for central atom is invalid')
    iCenter = center - 1 
    # Recenter and order by distance
    shift = lambda x: map(lambda (a,b): a-b, zip(x,dym['atCoords'][iCenter]))
    swapcoords = [[i,bohr2angstrom(shift(dym['atCoords'][i]))] for i in xrange(nAt)]
    for a in swapcoords:
        a.append(sqrt(sum(x*x for x in a[1])))
    swapcoords.sort(key=itemgetter(2))
    dym['printOrder'] = [swapcoords[iAt][0] for iAt in xrange(nAt)]
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
        isStr = lambda x: isinstance(x, basestring)
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
    from dmdw import writeDym

    # Check that we have a valid center atom number
    nAt = dym['nAt']
    if center > nAt or center < 1:
        raise ValueError('Index for central atom is invalid')
    iCenter = center - 1 

    # Recenter and order by distance
    shift = lambda x: map(lambda (a,b): a-b, zip(x,dym['atCoords'][iCenter]))
    swapcoords = [[i,shift(dym['atCoords'][i])] for i in xrange(nAt)]
    for a in swapcoords:
        a.append(sqrt(sum(x*x for x in a[1])))
    swapcoords.sort(key=itemgetter(2))
    dym['printOrder'] = [swapcoords[j][0] for j in xrange(nAt)]
    
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
