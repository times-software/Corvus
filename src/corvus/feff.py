from structures import Handler, Exchange, Loop, Update
import corvutils.pyparsing as pp
import os, sys, subprocess, shutil
import corvutils.corvus_misc_funcs as cu

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



class Feff(Handler):
    def __str__(self):
        return 'FEFF Handler'

    @staticmethod
    def Produces():
      return implemented

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
    def sequenceFor(output):
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

        import corvutils.corvus_misc_funcs as cu

        if config['verbose'] > 0:
          print 'Entering Handler {0}'.format(Feff.__name__)

        # set atoms and potentials

        # Set directory to feff executables.
# Debug: FDV
#       pp_debug.pprint(config)
        feffdir = config['feff']
# Debug: FDV
#       sys.exit()
        
        # Loop over targets in output. Not sure if there will ever be more than one output target here.
        for target in output:
            if (target == 'feffAtomicData'):
                # First generate any data that is needed from input
                input['feff.atoms']      = getFeffAtomsFromCluster(input)
                input['feff.potentials'] = getFeffPotentialsFromCluster(input)
               
                # Set directory for this exchange
                dir = config['xcDir']

                # Set output and error files
                with open(os.path.join(dir, 'corvus.FEFF.stdout'), 'w') as out, open(os.path.join(dir, 'corvus.FEFF.stderr'), 'w') as err:

                    # Set input file
                    inpf = os.path.join(dir, 'feff.inp')
                    
                    # Write input file for FEFF.
                    writeAtomicInput(input,inpf)

                    # Loop over executable: This is specific to feff. Other codes
                    # will more likely have only one executable.
                    # Run rdinp and atomic part of calculation
                    executables = ['rdinp','atomic','screen']
# Modified by FDV:
                    for executable in executables:
                        cu.subp(config,[feffdir+'/'+executable],dir,out,err)

                # Translate files produced by executable and update output
                # Normally this will look like the following example:
                             
                # filename = os.path.join(dir, 'xmu.dat')
                # output[target] = readColumns(filename, columns=[1,4])

                # For this case, I am only passing the directory for now so
                # that other executables in FEFF can use the atomic data.

                output[target] = dir

            elif (target == 'feffSCFPotentials'):
                # First generate any data that is needed from input
                input['feff.atoms']      = getFeffAtomsFromCluster(input)
                input['feff.potentials'] = getFeffPotentialsFromCluster(input)
                
                # Set directory for this exchange
                dir = config['xcDir']

                # Set output and error files
                with open(os.path.join(dir, 'corvus.FEFF.stdout'), 'w') as out, open(os.path.join(dir, 'corvus.FEFF.stderr'), 'w') as err:

                    # Set input file
                    inpf = os.path.join(dir, 'feff.inp')

                    # For now, copy all files from previous feff directory to current xcDir
                    previousDir = input['feffAtomicData']
                    src_files = os.listdir(previousDir)
                    for file_name in src_files:
                        full_file_name = os.path.join(previousDir, file_name)
                        if (os.path.isfile(full_file_name)):
                            shutil.copy(full_file_name, dir)

                    # Write input file for FEFF.
                    writeSCFInput(input,inpf)

                    # Loop over executable: This is specific to feff. Other codes
                    # will more likely have only one executable. Here I am running 
                    # rdinp again since writeSCFInput may have different cards than

                    # Run rdinp and atomic part of calculation
                    executables = ['rdinp','pot']
# Modified by FDV:
                    for executable in executables:
                        cu.subp(config,[feffdir+'/'+executable],dir,out,err)

                # Translate files produced by executable and update output
                # Normally this will look like the following example:

                # filename = os.path.join(dir, 'xmu.dat')
                # output[target] = readColumns(filename, columns=[1,4])

                # For this case, I am only passing the directory for now so
                # that other executables in FEFF can use the atomic data.

                output[target] = dir
            elif (target == 'feffCrossSectionsAndPhases'):
                # First generate any data that is needed from input
                input['feff.atoms']      = getFeffAtomsFromCluster(input)
                input['feff.potentials'] = getFeffPotentialsFromCluster(input)

                # Set directory for this exchange
                dir = config['xcDir']

                # Set output and error files
                with open(os.path.join(dir, 'corvus.FEFF.stdout'), 'w') as out, open(os.path.join(dir, 'corvus.FEFF.stderr'), 'w') as err:

                    # Set input file
                    inpf = os.path.join(dir, 'feff.inp')

                    # For now, copy all files from previous feff directory to current xcDir
                    previousDir = input['feffSCFPotentials']
                    src_files = os.listdir(previousDir)
                    for file_name in src_files:
                        full_file_name = os.path.join(previousDir, file_name)
                        if (os.path.isfile(full_file_name)):
                            shutil.copy(full_file_name, dir)

                    # Write input file for FEFF.
                    writeCrossSectionsInput(input,inpf)

                    # Loop over executable: This is specific to feff. Other codes
                    # will more likely have only one executable. Here I am running 
                    # rdinp again since writeSCFInput may have different cards than

                    # Run rdinp and atomic part of calculation
                    executables = ['rdinp','xsph']
# Modified by FDV:
                    for executable in executables:
                        cu.subp(config,[feffdir+'/'+executable],dir,out,err)

                # Translate files produced by executable and update output
                # Normally this will look like the following example:

                # filename = os.path.join(dir, 'xmu.dat')
                # output[target] = readColumns(filename, columns=[1,4])

                # For this case, I am only passing the directory for now so
                # that other executables in FEFF can use the atomic data.

                output[target] = dir

 
            elif (target == 'feffGreensFunction'):
                # First generate any data that is needed from input
                input['feff.atoms']      = getFeffAtomsFromCluster(input)
                input['feff.potentials'] = getFeffPotentialsFromCluster(input)

                # Set directory for this exchange
                dir = config['xcDir']

                # Set output and error files
                with open(os.path.join(dir, 'corvus.FEFF.stdout'), 'w') as out, open(os.path.join(dir, 'corvus.FEFF.stderr'), 'w') as err:

                    # Set input file
                    inpf = os.path.join(dir, 'feff.inp')

                    # For now, copy all files from previous feff directory to current xcDir
                    previousDir = input['feffCrossSectionsAndPhases']
                    src_files = os.listdir(previousDir)
                    for file_name in src_files:
                        full_file_name = os.path.join(previousDir, file_name)
                        if (os.path.isfile(full_file_name)):
                            shutil.copy(full_file_name, dir)

                    # Write input file for FEFF.
                    writeGreensFunctionInput(input,inpf)

                    # Loop over executable: This is specific to feff. Other codes
                    # will more likely have only one executable. Here I am running 
                    # rdinp again since writeSCFInput may have different cards than

                    # Run rdinp and atomic part of calculation
                    executables = ['rdinp','fms','mkgtr']
# Modified by FDV:
                    for executable in executables:
                        cu.subp(config,[feffdir+'/'+executable],dir,out,err)

                # Translate files produced by executable and update output
                # Normally this will look like the following example:

                # filename = os.path.join(dir, 'xmu.dat')
                # output[target] = readColumns(filename, columns=[1,4])

                # For this case, I am only passing the directory for now so
                # that other executables in FEFF can use the atomic data.

                output[target] = dir


            elif (target == 'feffPaths'):
                # First generate any data that is needed from input
                input['feff.atoms']      = getFeffAtomsFromCluster(input)
                input['feff.potentials'] = getFeffPotentialsFromCluster(input)

                # Set directory for this exchange
                dir = config['xcDir']

                # Set output and error files
                with open(os.path.join(dir, 'corvus.FEFF.stdout'), 'w') as out, open(os.path.join(dir, 'corvus.FEFF.stderr'), 'w') as err:

                    # Set input file
                    inpf = os.path.join(dir, 'feff.inp')

                    # For now, copy all files from previous feff directory to current xcDir
                    previousDir = input['feffGreensFunction']
                    src_files = os.listdir(previousDir)
                    for file_name in src_files:
                        full_file_name = os.path.join(previousDir, file_name)
                        if (os.path.isfile(full_file_name)):
                            shutil.copy(full_file_name, dir)

                    # Write input file for FEFF.
                    writePathsInput(input,inpf)

                    # Loop over executable: This is specific to feff. Other codes
                    # will more likely have only one executable. Here I am running 
                    # rdinp again since writeSCFInput may have different cards than

                    # Run rdinp and atomic part of calculation
                    executables = ['rdinp','path']
# Modified by FDV:
                    for executable in executables:
                        cu.subp(config,[feffdir+'/'+executable],dir,out,err)

                # Translate files produced by executable and update output
                # Normally this will look like the following example:

                # filename = os.path.join(dir, 'xmu.dat')
                # output[target] = readColumns(filename, columns=[1,4])

                # For this case, I am only passing the directory for now so
                # that other executables in FEFF can use the atomic data.

                output[target] = dir

            elif (target == 'feffFMatrices'):
                # First generate any data that is needed from input
                input['feff.atoms']      = getFeffAtomsFromCluster(input)
                input['feff.potentials'] = getFeffPotentialsFromCluster(input)

                # Set directory for this exchange
                dir = config['xcDir']

                # Set output and error files
                with open(os.path.join(dir, 'corvus.FEFF.stdout'), 'w') as out, open(os.path.join(dir, 'corvus.FEFF.stderr'), 'w') as err:

                    # Set input file
                    inpf = os.path.join(dir, 'feff.inp')

                    # For now, copy all files from previous feff directory to current xcDir
                    previousDir = input['feffPaths']
                    src_files = os.listdir(previousDir)
                    for file_name in src_files:
                        full_file_name = os.path.join(previousDir, file_name)
                        if (os.path.isfile(full_file_name)):
                            shutil.copy(full_file_name, dir)

                    # Write input file for FEFF.
                    writeFMatricesInput(input,inpf)

                    # Loop over executable: This is specific to feff. Other codes
                    # will more likely have only one executable. Here I am running 
                    # rdinp again since writeSCFInput may have different cards than

                    # Run rdinp and atomic part of calculation
                    executables = ['rdinp','genfmt']
# Modified by FDV:
                    for executable in executables:
                        cu.subp(config,[feffdir+'/'+executable],dir,out,err)

                # Translate files produced by executable and update output
                # Normally this will look like the following example:

                # filename = os.path.join(dir, 'xmu.dat')
                # output[target] = readColumns(filename, columns=[1,4])

                # For this case, I am only passing the directory for now so
                # that other executables in FEFF can use the atomic data.

                output[target] = dir

            elif (target == 'feffXANES'):
                # First generate any data that is needed from input
                input['feff.atoms']      = getFeffAtomsFromCluster(input)
                input['feff.potentials'] = getFeffPotentialsFromCluster(input)

                # Set directory for this exchange
                dir = config['xcDir']

                # Set output and error files
                with open(os.path.join(dir, 'corvus.FEFF.stdout'), 'w') as out, open(os.path.join(dir, 'corvus.FEFF.stderr'), 'w') as err:

# Debug FDV
#                   print os.path.join(dir, 'corvus.FEFF.stdout')
#                   print os.path.join(dir, 'corvus.FEFF.stderr')
                    # Set input file
                    inpf = os.path.join(dir, 'feff.inp')

                    # For now, copy all files from previous feff directory to current xcDir
                    # For now, we are going to run all modules for XANES.
                    #previousDir = input['feffFMatrices']
                    #src_files = os.listdir(previousDir)
                    #for file_name in src_files:
                    #    full_file_name = os.path.join(previousDir, file_name)
                    #    if (os.path.isfile(full_file_name)):
                    #        shutil.copy(full_file_name, dir)

                    # Write input file for FEFF.
                    writeXANESInput(input,inpf)

                    # Loop over executable: This is specific to feff. Other codes
                    # will more likely have only one executable. Here I am running 
                    # rdinp again since writeSCFInput may have different cards than

                    # Run rdinp and atomic part of calculation
                    executables = ['rdinp','atomic','pot','screen','opconsat','xsph','fms','mkgtr','path','genfmt','ff2x','sfconv']
# Modified by FDV:
                    for executable in executables:
                        cu.subp(config,[feffdir+'/'+executable],dir,out,err)

                # Translate files produced by executable and update output
                # Normally this will look like the following example:

                # filename = os.path.join(dir, 'xmu.dat')
                # output[target] = readColumns(filename, columns=[1,4])

                # For this case, I am only passing the directory for now so
                # that other executables in FEFF can use the atomic data.

# FDV: Testing the generation of gnuplot scripts out of xmu.dat files
                filename = os.path.join(dir, 'xmu.dat')
# Debug: FDV
#               pp_debug.pprint(input)
#               print input['mkgnuplot'][0][0]
                if input['mkgnuplot'][0][0]:
                  Ctrl_Dat = { 'Plot_mu'  : [ [1,4], ["Energy (eV)", "Mu"] ],
                               'Plot_mu0' : [ [1,5], ["Energy (eV)", "Mu0"] ],
                               'Plot_chi' : [ [1,6], ["Energy (eV)", "Chi"] ] }
                  cu.Make_Gnuplot(filename,Ctrl_Dat)

                output[target] = dir


            elif (target == 'feffXES'):
                # First generate any data that is needed from input
                input['feff.atoms']      = getFeffAtomsFromCluster(input)
                input['feff.potentials'] = getFeffPotentialsFromCluster(input)

                # Set directory for this exchange
                dir = config['xcDir']

                # Set output and error files
                with open(os.path.join(dir, 'corvus.FEFF.stdout'), 'w') as out, open(os.path.join(dir, 'corvus.FEFF.stderr'), 'w') as err:

                    # Set input file
                    inpf = os.path.join(dir, 'feff.inp')

                    # For now, copy all files from previous feff directory to current xcDir
                    # For now, we are going to run all modules for XANES.
                    #previousDir = input['feffFMatrices']
                    #src_files = os.listdir(previousDir)
                    #for file_name in src_files:
                    #    full_file_name = os.path.join(previousDir, file_name)
                    #    if (os.path.isfile(full_file_name)):
                    #        shutil.copy(full_file_name, dir)

                    # Write input file for FEFF.
                    writeXESInput(input,inpf)

                    # Loop over executable: This is specific to feff. Other codes
                    # will more likely have only one executable. Here I am running 
                    # rdinp again since writeSCFInput may have different cards than

                    # Run rdinp and atomic part of calculation
                    executables = ['rdinp','atomic','pot','screen','opconsat','xsph','fms','mkgtr','path','genfmt','ff2x','sfconv']
# Modified by FDV:
                    for executable in executables:
                        cu.subp(config,[feffdir+'/'+executable],dir,out,err)

                # Translate files produced by executable and update output
                # Normally this will look like the following example:

                # filename = os.path.join(dir, 'xmu.dat')
                # output[target] = readColumns(filename, columns=[1,4])

                # For this case, I am only passing the directory for now so
                # that other executables in FEFF can use the data.

                output[target] = dir

        if config['verbose'] > 0:
          print 'Done with Handler {0}'.format(Feff.__name__)

    @staticmethod
    def cleanup(config):
        pass



##### Generic Helper Methods ##########

def check(input, token, default=None):
    if token in input:
        return input[token]
    else:
        return default

def getInpLines(input,token,default=None):
    lines=[]
    key = token[len('feff.'):]

    if token in input:
        # If the first element is not a boolean, this contains values
        # to be stored after keyword.
        if not isinstance(input[token][0][0],bool):
            if token in input:
                for element in input[token]:
                    lines.append(' '.join([str(value) for value in element])) 
            else:
                lines.append(' '.join(default))
          	
            if key in ['atoms','potentials','egrid']:
                lines.insert(0,key.upper())
            else:
                lines[0] = key.upper() + ' ' + lines[0]
        elif input[token][0][0]:
            lines.append(key.upper())

    elif default:
        # Put default into input dictionary
   
        if key in ['egrid']:
           # for multi line input (fill with different targets later.
           lines = [x for x in default]
           lines.insert(0,key.upper())

        else:
           lines.append(key.upper() + ' ' + default)

    else:
        lines.append(key.upper())

    # Add a blank line after each line
    lines.append('')

    return lines

def writeList(lines, filename):
    with open(filename, 'w') as f:
        f.write('\n'.join(lines))

def readColumns(filename, columns=[1,2]):
    # Read file and clear out comments
    with open(filename, 'r') as file:
        cleanStr = file.read()
    comments = pp.ZeroOrMore(pp.pythonStyleComment).setParseAction(pp.replaceWith(''))
    try:
        cleanStr = comments.transformString(cleanStr)
    except pp.ParseException as pe:
        print 'Parsing Error using pyparsing: invalid input:', pe
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
        print 'Parsing Error using pyparsing: invalid input:', pe
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
        print "Parsing Error using pyparsing: invalid input:", pe
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
    
    # Set feff.KEY using getInpLine, which will use the
    # previously defined value if available.
    lines = lines + getInpLines(input,'feff.control', '1 0 0 0 0 0')
    lines = lines + getInpLines(input,'feff.print', '5 0 0 0 0 0')
    
    # loop over feff keys in input and write them to lines
    feffKeys = [key for key in input if key.startswith('feff.')]

    for key in feffKeys:
        lines = lines + getInpLines(input,key)

    lines = lines + getInpLines(input,'feff.end')
    
    # Print feff input file
    writeList(lines, feffinp)


def writeSCFInput(input, feffinp='feff.inp'):
    for key in ['feff.potentials', 'feff.atoms']:
        assert key in input
    lines = []
    lines.append('* This feff9 input file was generated by corvus.')
    lines.append('')

    lines = lines + getInpLines(input,'feff.control', '1 0 0 0 0 0')
    lines = lines + getInpLines(input,'feff.print', '5 0 0 0 0 0')

    # Header part of feff.inp
    #headerLines(input, lines)

    lines = lines + getInpLines(input,'feff.scf','5.0 0 100 0.1 0')
    lines = lines + getInpLines(input,'feff.edge','K')
    
    # loop over feff keys in input and write them to lines
    feffKeys = [key for key in input if key.startswith('feff.')]

    for key in feffKeys:
        lines = lines + getInpLines(input,key)

    lines = lines + getInpLines(input,'feff.end')

    # Print feff input file
    writeList(lines, feffinp)
    #print input

def writeCrossSectionsInput(input, feffinp='feff.inp'):
    for key in ['feff.potentials', 'feff.atoms']:
        assert key in input
    lines = []
    lines.append('* This feff9 input file was generated by corvus.')
    lines.append('')

    lines = lines + getInpLines(input,'feff.control', '0 1 0 0 0 0')
    lines = lines + getInpLines(input,'feff.print', '5 2 0 0 0 0')
    lines = lines + getInpLines(input,'feff.exchange', '0 0.0 0.0 0')

    # Header part of feff.inp
    #headerLines(input, lines)

    lines = lines + getInpLines(input,'feff.edge','K')

    # loop over feff keys in input and write them to lines
    feffKeys = [key for key in input if key.startswith('feff.')]

    for key in feffKeys:
        lines = lines + getInpLines(input,key)

    lines = lines + getInpLines(input,'feff.end')

    # Print feff input file
    writeList(lines, feffinp)


def writeGreensFunctionInput(input, feffinp='feff.inp'):
    for key in ['feff.potentials', 'feff.atoms']:
        assert key in input
    lines = []
    lines.append('* This feff9 input file was generated by corvus.')
    lines.append('')

    lines = lines + getInpLines(input,'feff.control', '0 0 1 0 0 0')
    
    # loop over feff keys in input and write them to lines
    feffKeys = [key for key in input if key.startswith('feff.')]

    for key in feffKeys:
        lines = lines + getInpLines(input,key)

    lines = lines + getInpLines(input,'feff.end')

    # Print feff input file
    writeList(lines, feffinp)


def writePathsInput(input, feffinp='feff.inp'):
    for key in ['feff.potentials', 'feff.atoms']:
        assert key in input
    lines = [] 
    lines.append('* This feff9 input file was generated by corvus.')
    lines.append('')

    lines = lines + getInpLines(input,'feff.control', '0 0 0 1 0 0')


    # loop over feff keys in input and write them to lines
    feffKeys = [key for key in input if key.startswith('feff.')]

    for key in feffKeys:
        lines = lines + getInpLines(input,key)

    lines = lines + getInpLines(input,'feff.end')

    # Print feff input file
    writeList(lines, feffinp)


def writeFMatricesInput(input, feffinp='feff.inp'):
    for key in ['feff.potentials', 'feff.atoms']:
        assert key in input
    lines = []
    lines.append('* This feff9 input file was generated by corvus.')
    lines.append('')

    lines = lines + getInpLines(input,'feff.control', '0 0 0 0 1 0')


    # loop over feff keys in input and write them to lines
    feffKeys = [key for key in input if key.startswith('feff.')]

    for key in feffKeys:
        lines = lines + getInpLines(input,key)

    lines = lines + getInpLines(input,'feff.end')

    # Print feff input file
    writeList(lines, feffinp)


def writeXANESInput(input, feffinp='feff.inp'):
    for key in ['feff.potentials', 'feff.atoms']:
        assert key in input
    lines = []
    lines.append('* This feff9 input file was generated by corvus.')
    lines.append('')

    lines = lines + getInpLines(input,'feff.control', '1 1 1 1 1 1')
    lines = lines + getInpLines(input,'feff.print', '0 0 0 0 0 0')
    lines = lines + getInpLines(input,'feff.exchange', '0 0.0 0.0 0')
    lines = lines + getInpLines(input,'feff.xanes', '10')
    lines = lines + getInpLines(input,'feff.fms', '6.0 0')
    lines = lines + getInpLines(input,'feff.scf', '5.0 0')
    lines = lines + getInpLines(input,'feff.rpath', '0.1')
    lines = lines + getInpLines(input,'feff.edge', 'K')


    # loop over feff keys in input and write them to lines
    feffKeys = [key for key in input if key.startswith('feff.')]

    for key in feffKeys:
        lines = lines + getInpLines(input,key)

    lines = lines + getInpLines(input,'feff.end')

    # Print feff input file
    writeList(lines, feffinp)


def writeXESInput(input, feffinp='feff.inp'):
    for key in ['feff.potentials', 'feff.atoms']:
        assert key in input
    lines = []
    lines.append('* This feff9 input file was generated by corvus.')
    lines.append('')

    lines = lines + getInpLines(input,'feff.control', '1 1 1 1 1 1')
    lines = lines + getInpLines(input,'feff.print', '0 0 0 0 0 0')
    lines = lines + getInpLines(input,'feff.exchange', '0 0.0 0.0 0')
    lines = lines + getInpLines(input,'feff.xes', '-30 0')
    lines = lines + getInpLines(input,'feff.fms', '6.0 0')
    lines = lines + getInpLines(input,'feff.scf', '5.0 0')
    lines = lines + getInpLines(input,'feff.rpath', '0.1')
    lines = lines + getInpLines(input,'feff.edge', 'K')


    # loop over feff keys in input and write them to lines
    feffKeys = [key for key in input if key.startswith('feff.')]

    for key in feffKeys:
        lines = lines + getInpLines(input,key)

    lines = lines + getInpLines(input,'feff.end')

    # Print feff input file
    writeList(lines, feffinp)


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
