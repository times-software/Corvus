import sys, os
import imp
# Debug: FDV
import pprint
pp_debug = pprint.PrettyPrinter(indent=4)

# Define the available handlers by hand here
# Note FDV: There should be a way to do this automatically, but can't think of
# one right now
def availableHandlers():
    # Need to check which of these is defined in config file befor importing
    # and adding to list of available handlers. 
    config = {}
    configure(config)
    handlers = []
    if config['abinit'] and config['anaddb'] and config['mrgddb'] and config['mrggkk']:
        from abinit import Abinit
        handlers = handlers + [Abinit]
    if config['dmdw']:
        from dmdw import Dmdw
        handlers = handlers + [Dmdw]
    if config['feff']:
        from feff import Feff
        handlers = handlers + [Feff]
    if config['vasp_gam'] and config['vasp_std']:
        from vasp import Vasp
        handlers = handlers + [Vasp]
    if config['nwchem']:
        from nwchem import Nwchem
        handlers = handlers + [Nwchem]
    if config['orca']: 
        # JK - Note that there is another command named orca on many linux systems.
        # Need to figure out how to check that this is indeed orca electronic structure code.
        from orca import Orca
        handlers = handlers + [Orca]
    if config['siesta']:
        from siesta import Siesta
        handlers = handlers + [Siesta]
    # import only if module lmfit exists (fit dependency). Should probably
    # do this with numpy and scipy as well. 
    try:
        import lmfit
        from fit import fit
        handlers = handlers + [fit]
    except ImportError:
        print("Warning: lmfit not found. fit handler will be disabled.")
        pass


    return handlers

# Commenting out Vasp for now since we have no content in the manual for it for
# now. I will readd once I include Scott's Handler.
#   return [Feff, FeffRixs, Dmdw, Abinit, Vasp, Nwchem, Orca]

def configure(config):
    from ConfigParser import RawConfigParser
    from platform import system
    
    # Store path to Corvus executables
    config['bin'] = os.path.dirname(os.path.abspath(__file__))
    # Store path to current working directory (where user is running Corvus)
    config['cwd'] = os.getcwd()

    # Load Corvus defaults
    rcp = RawConfigParser(allow_no_value=True)
    configFile = os.path.join(os.path.dirname(config['bin']), 'corvus.conf')
    if configFile not in rcp.read(configFile):
        printAndExit('Error reading corvus.conf')
    config['pathprefix'] = rcp.get('Defaults', 'prefix')
    config['inputsuffix'] = rcp.get('Defaults', 'inputsuffix')
    config['savesuffix'] = rcp.get('Defaults', 'savesuffix')
    config['inputFile'] = config['pathprefix'] + config['inputsuffix']
    config['saveFile'] = config['pathprefix'] + config['savesuffix']
    config['checkpoints'] = rcp.getboolean('Defaults', 'checkpoints') 
    config['parallelRun'] = rcp.get('Defaults', 'parallelrun')
    utilpath = os.path.join(os.path.dirname(config['bin']), 'corvutils')
# Modified by FDV
    config['parsnipConf'] = os.path.join(utilpath, 'parsnip.corvus.config')
    config['parsnipForm'] = os.path.join(utilpath, 'parsnip.corvus.formats')
    
# Here we add a check to see if we are running under Cygwin. In Cygwin
# executables have a ".exe" extension. We also define all the different
# executables in one spot to make things more clear, and generalize the code.

    Execs_Dict = { 'feff'  :['atomic','band','compton','crpa','dmdw','eels',
                             'ff2x','fms','fullspectrum','genfmt','ldos',
                             'mkgtr','opconsat','path','pot','rdinp','rhorrp',
                             'rixs','screen','sfconv','xsph'],
                   'dmdw'  :['dmdw'],
                   'orca'  :['orca'],
                   'abinit':['abinit','anaddb','mrgddb','mrggkk'],
                   'vasp'  :['vasp_gam','vasp_std'],
                   'nwchem':['nwchem'], 
                   'siesta':['siesta']}


    exe_ext = ''
    if system()[0:6]=='CYGWIN':
      exe_ext = '.exe'

    for code in Execs_Dict.keys():
      path2code = rcp.get('Executables', code)
      for cmd in Execs_Dict[code]:
        config[cmd] = which(cmd+exe_ext, path2code)

# NOTE FDV: Here I try to add code to make JJK's version work. We need the path
# to feff to be avilable in config.
    path2feff = rcp.get('Executables', 'feff')
    config['feff'] = path2feff

    path2siesta = rcp.get('Executables', 'siesta') 
    config['siesta'] = path2siesta
# Initialize system with user input
def initializeSystem(config, system):
    from corvutils import parsnip

    conf    = checkFile(config['parsnipConf'])
    inp     = config['inputFile']

# Debug: FDV
#print parsnip.parse(conf, inp, mode=['read','UseDefaults'])
#   sys.exit()

    if os.path.exists(inp):
        system.update(parsnip.parse(conf, inp, mode=['read','UseDefaults']))

# Debug: FDV
#   pp_debug.pprint(system)
#   sys.exit()

# Basic Workflow Generator
# JK - added system as input to generateWorkflow to accomodate handlers that call other 
# handlers or workflows such as "fit" or "average".
# JK - added config as input to allow use of paths to check for necessary executables
# before adding handlers to available handler list.
def generateWorkflow(target, handlers, system, config, desc=''):
    from structures import Workflow

# Note FDV: Moving to top so this is available throughout
# Define the available handlers by hand here
#   availableHandlers = [Feff, FeffRixs, Dmdw, Abinit, Vasp, Nwchem, Orca]
    #availableHandlers = [Nwchem, Dmdw]
    #availableHandlers = [Abinit, Dmdw]
# From the handlers list we generate a mapping dictionary to identify the 
# handlers requested by the user
    availHandlers = availableHandlers()
    availableHandlers_map = { h.__name__:h for h in availHandlers }
# Debug: FDV
    #pp_debug.pprint(availableHandlers_map)
#   sys.exit()

# Now we create the useHandlers list from the input and the available ones
    useHandlers = []
# Debug: FDV
    #print handlers
    #print availableHandlers_map.keys()
    for handler_name in handlers:
      if handler_name in availableHandlers_map.keys():
        useHandlers.append(availableHandlers_map[handler_name])
        #print useHandlers
      else:
        print("Handler %s not in list of available handlers:" % (handler_name))
        for s in availableHandlers_map.keys():
          print("%s" % s)
        sys.exit()

# Debuf: FDV
#   pp_debug.pprint(availableHandlers)
#   pp_debug.pprint(useHandlers)
#   sys.exit()
    workflow = Workflow(target, desc=desc)
    # subs creates a list of all possible subsets of the set of targets.
    subs = lambda L:[{L[j] for j in range(len(L)) if 1<<j&k} for k in range(1,1<<len(L))]
# Modified by FDV:
# Commenting original code and substituting with JJKs code, to see if it works
#-----------------------------------------------------------------------------
#    targets = set(workflow.getRequiredInput())
#    while len(targets) > 0:
#        noMatch = True
#        for h in useHandlers:
#            htargets = [t for t in targets if h.canProduce(t)]
#            if htargets:
#                noMatch = False
#                for subset in reversed(sorted(subs(htargets), key=len)):
#                    l = list(subset)
#                    if h.canProduce(l):
#                        workflow.addExchangeAt(0, h.sequenceFor(l))
#                        targets.difference_update(subset)
#                        targets.update(set(workflow.getRequiredInput()))
#                        break
#        if noMatch:
#            break
# Debug: FDV
#   berp = {}
#workflow print statements uncommented, Krsna
#   print type(workflow.sequence[0])
#   print workflow.sequence[0].handler
#   workflow.sequence[0]({1:'a'},{1:'a'},berp)
#   sys.exit()
#-----------------------------------------------------------------------------
    mainTargets = workflow.getRequiredInput()
    #print mainTargets
    #print "Done printing ri" # Debug JJK
    for t in reversed(mainTargets):
       targets = set([t])
       while len(targets) > 0:
           noMatch = True
           # JK - Make sure there is some handler that can produce each target.
           #for target in targets:
           #    handler_not_found = True
           #    for h in availHandlers:
           #        if h.canProduce(target):
           #            handler_not_found = False
           #    if handler_not_found:
           #        print("No handler to produce target:", target)
           #        print("Install the correct software if it exists.")
           #        sys.exit()



           for h in availHandlers:
               # In the below line we will want to check if target is a list, and if so, check
               # if h can produce any elements of target. Optionally, we can check 
               htargets = [target for target in targets if h.canProduce(target)]
               if htargets:
                   noMatch = False
                   for subset in reversed(sorted(subs(htargets), key=len)):
                       l = list(subset)
                       if h.canProduce(l):
                           workflow.addExchangeAt(0, h.sequenceFor(l,system))
                           # h.setDefaults(input,h)
                           # JK - We need to replace the following line with code that
                           # checks a list of lists rather than a list. If any requirement in an
                           # inner list is found, the entire inner list can be removed.
                           # For each target,
                           #for target in targets:
                           #    if isinstance(target,list):

                           targets.difference_update(subset)  
                           targets.update(set(workflow.getRequiredInput()))
                           break

           if noMatch:
               break

#-----------------------------------------------------------------------------

# Debug: FDV
#    print '#########################'
#    print '    WORKFLOW'
#    print '#########################'
#    print workflow
#    print '#########################'
#    print '  END WORKFLOW'
#    print '#########################'
#   sys.exit()

    return workflow

# Formats output strings
def prettyprint(token, config, system):
    import pprint
    # JK - Going to leave this for now, but output should probably
    # depend on the target. 
    f = open(config['pathprefix'] + '.' + token + '.out', 'w')
    data = system[token]

    if token.startswith('Table'): # Table:col,col,col
        labels = token[token.index(':')+1:].split(',')
        f.write('#' + ','.join(labels) + '\n')
        cols = []
        for l in labels:
            cols.append(data[l])
        for row in zip(*cols):
            f.write('    '.join(map(str,row)) + '\n')
    elif isinstance(data, list):
        # We have columns of data
        if isinstance(data[0], list) and len(data[0]) == len(data[1]):
            for row in zip(*data):
                f.write('    '.join(map(str,row)) + '\n')
        # Default to Python printing
        else: 
            pprint.pprint(data, stream=f)
    # Default to Python printing
    else:
        pprint.pprint(data, stream=f)

def which(program, search=None):
    import os
    is_exe = lambda f: os.path.isfile(f) and os.access(f, os.X_OK)
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        if search:
            for root, dirs, files in os.walk(os.path.expanduser(search)):
                if program in files:
                    filepath = os.path.join(root, program)
                    if is_exe(filepath):
                        return filepath
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            filepath = os.path.join(path, program)
            if is_exe(filepath):
                return filepath
    return None

def usage():
    print("Usage: corvus.py [options]")
    print("       Possible options:")
#uncommented print(line below, Krsna)
    print("         -t, --target      [comma-separated list]")
    print("         -w, --workflow    [filename]")
    print("         -i, --input       [filename]"  )
    print("         -c, --checkpoints"  )
    print("         -r, --resume"  )
    print("         -s, --save        [filename]"  )
    print("             --prefix      [alphanumeric string (-_. okay)]")
    print("             --parallelrun [command string]")

def checkFile(filename):
    if not os.path.exists(filename):
        printAndExit('File does not exist: ' + filename)
    else:
        return filename

def printAndExit(msg):
    sys.stderr.write(msg + '\n')
    sys.stderr.flush()
    sys.exit()

def generateAndRunWorkflow(config, system, targetList):
    import copy
    #if 'target_list' in system.keys():
    #  targetList = system['target_list']
    #else:
    #  print 'Provide target properties or Workflow'
    #  sys.exit()

    if 'usehandlers' in system.keys():
      handlerList = system['usehandlers'][0]
    else:
      print('Provide the handler list to be used in this calculation')
      sys.exit()

    autodesc = 'Calculate ' + ', '.join(targetList[0])
    workflow = generateWorkflow(targetList, handlerList, system, config, desc=autodesc)

    # Check for any missing user input or handlers.
    required = set(workflow.getRequiredInput())
    missing = list(required.difference(set(system.keys())))
    
    #config2=copy.deepcopy(config)
    # Create a copy of config
    if len(missing) > 0:
        printAndExit('Error: missing user input or handler for ' + str(missing))

    i = 0
    while i < len(workflow.sequence):
        #if config['checkpoints']:
        #    saveState['system'] = system
        #    saveState['index'] = i
        #    with open(config['saveFile'], 'w') as saveFile:
        #        pickle.dump(saveState, saveFile, pickle.HIGHEST_PROTOCOL)
        #print 'Workflow index'
        #print config['xcIndex']
        #print config['xcIndexStart']
        config['xcIndex'] = i + config['xcIndexStart']
        #print config['xcIndex']
        
        workflow.sequence[i].go(config, system)
        i += 1

# Handles command line arguments and runs through workflow
def oneshot(argv):
    import getopt, pickle

# Debug:FDV
# Test the writeDict function in abinit
#   from abinit import writeDict
# Create a test dictionary
#   dict={'aa':'12\n34','a':1,'b':[[1],[2]],'c':[[1]],'d':[1,2]}
#   writeDict(dict,'pepe')
#   sys.exit()

# Modified by FDV:
# Removing the target option from the cli. From now on we do it through the
# input.
#   shortopts = 'crt:i:w:s:j:'
#   longopts = ['target=','input=','workflow=','checkpoints','resume','save=',
#               'jump=','prefix=','parallelrun=']
    shortopts = 'cr:i:w:s:j:'
    longopts = ['input=','workflow=','checkpoints','resume','save=',
                'jump=','prefix=','parallelrun=']
    try:
        opts, args = getopt.getopt(argv[1:], shortopts, longopts)
    except getopt.GetoptError as err:
        print(str(err))
        usage()
        sys.exit(2)

    ### Read user options ###
    resume = checkpoints = False
    jumpPoint = 0
# Modified by FDV:
# Removing the target option from the cli. From now on we do it through the
# input.
#   inputFile = saveFile = workflowFile = targetList = pathPrefix = parallelRun = None
    inputFile = saveFile = workflowFile = pathPrefix = parallelRun = None
    for opt, arg in opts:
#       if opt in ('-t', '--target'):
#           targetList = arg.split(',')
#       elif opt in ('-i', '--input'):
#           inputFile = checkFile(arg)
        if opt in ('-i', '--input'):
            inputFile = checkFile(arg)
        elif opt in ('-w', '--workflow'):
            workflowFile = checkFile(arg) 
        elif opt in ('-c', '--checkpoints'):
            checkpoints = True
        elif opt in ('-s', '--save'):
            saveFile = checkFile(arg)
        elif opt in ('-r', '--resume'):
            resume = True
        elif opt in ('-j', '--jump'):
            if not arg.isdigit():
                printAndExit('Resume index should be a positive integer.')
            jumpPoint = int(arg) 
        elif opt in ('--prefix'):
            pathPrefix = arg 
        elif opt in ('--parallelrun'):
            parallelRun = arg

    # Set up the configuration and system database
    config = {}
    configure(config)
    # Overwrite defaults as requested
    if pathPrefix is not None:
        config['pathprefix'] = pathPrefix
        config['inputFile'] = config['pathprefix'] + config['inputsuffix']
        config['saveFile'] = config['pathprefix'] + config['savesuffix']
    if inputFile is not None:
        config['inputFile'] = inputFile
    if saveFile is not None:
        config['saveFile'] = saveFile
    if checkpoints:
        config['checkpoints'] = checkpoints
    if parallelRun is not None:
        config['parallelRun'] = parallelRun
#DASb
#line below uncommented by Krsna
#   print config
#DASe
# Debug: FDV
#   pp_debug.pprint(config)
#   sys.exit()

    system = {}
    workflowStart = 0

    # Read in saved state if restarting calculation midway
    if resume and os.path.exists(config['saveFile']):
        saveState = pickle.load(open(config['saveFile'], 'r'))
        system.update(saveState['system'])
        workflow = saveState['workflow']
        workflowStart = saveState['index']

    # Update System with any user input
    initializeSystem(config, system)
    #print(system)
#DASb
    #print 'system = ', system
#DASe
# Debug: FDV
#   print(system['target_list'])
#   sys.exit()

# Added by FDV
# At this point we set the target list based on the content of the input file,
# rather than the command line.
    if 'target_list' in system.keys():
      targetList = system['target_list']
    else:
      print('Provide target properties or Workflow')
      sys.exit()

# Added by FDV
# At this point we set the handler list to make the workflow generator work more
# generally.
    if 'usehandlers' in system.keys():
      handlerList = system['usehandlers'][0]
    else:
      print('Provide the handler list to be used in this calculation')
      sys.exit()

# Here we should probably add a check to see if the target list if empty, but
# I don't think it can happen. Leaving for future

# Debug: FDV
#print check targetList usehandlers by Krsna
    #print targetList
    #print handlerList
#   sys.exit()

    # Set up Workflow or read from file
    if workflowFile is not None:
        # Read Workflow from file
        workflow = pickle.load(open(workflowFile, 'r'))
    elif not resume:
        # Create Workflow
        autodesc = 'Calculate ' + ', '.join(targetList[0])
        workflow = generateWorkflow(targetList, handlerList, system, config, desc=autodesc)
        print('')
        print('')
        print('#########################')
        print('    WORKFLOW')
        print('#########################')
        print(workflow)
        print('#########################')
        print('  END WORKFLOW')
        print('#########################')
        print('')
        print('')

# Debug: FDV
#Krsna uncommented line below
#   print workflow
#   sys.exit()

    # Check for any missing user input
    required = set(workflow.getRequiredInput())
# Debug: FDV
#   print required
    missing = list(required.difference(set(system.keys())))
    if len(missing) > 0:
        printAndExit('Error: missing user input for ' + str(missing))

# Debug: FDV
#   sys.exit()

    # Skip to specified Workflow step (developer option)
    if jumpPoint is not 0:
        workflowStart = jumpPoint - 1

    # Save intial state
    saveState = {'system':system,'workflow':workflow}
    # Run Workflow
    i = max(0, workflowStart)
    #print workflow.sequence
    while i < len(workflow.sequence):
        if config['checkpoints']:
            saveState['system'] = system
            saveState['index'] = i
            with open(config['saveFile'], 'w') as saveFile:
                pickle.dump(saveState, saveFile, pickle.HIGHEST_PROTOCOL)
        
        config['xcIndex'] = i + 1
        workflow.sequence[i].go(config, system)
        i += 1
    # Save completed state
    saveState['system'] = system
    with open(config['saveFile'], 'w') as saveFile:
        pickle.dump(saveState, saveFile, pickle.HIGHEST_PROTOCOL)

    # Print output
    for token in workflow.target:
        if system[token] is not None:
            prettyprint(token, config, system)
        else: 
            printAndExit('Error: target [' + token + '] not produced')

