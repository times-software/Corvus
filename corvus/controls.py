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
        from corvus.abinit import Abinit
        handlers = handlers + [Abinit]
    if config['dmdw']:
        from corvus.dmdw import Dmdw
        handlers = handlers + [Dmdw]
    if config['feff']:
        from corvus.feff import Feff
        handlers = handlers + [Feff]
    if config['vasp_gam'] and config['vasp_std']:
        from corvus.vasp import Vasp
        handlers = handlers + [Vasp]
    if config['nwchem']:
        from corvus.nwchem import Nwchem
        handlers = handlers + [Nwchem]
    if config['ocean']:
        from corvus.ocean import Ocean
        handlers = handlers + [Ocean]
    if config['orca']: 
        # JK - Note that there is another command named orca on many linux systems.
        # Need to figure out how to check that this is indeed orca electronic structure code.
        from corvus.orca import Orca
        handlers = handlers + [Orca]
    if config['siesta']:
        from corvus.siesta import Siesta
        handlers = handlers + [Siesta]
    if config['phsf']:
        from corvus.phsf import phsf
        handlers = handlers + [phsf]
    if config['cif2cell']:
        from corvus.Cif2Cell import Cif2Cell
        handlers = handlers + [Cif2Cell]

    from corvus.mbconv import mbconv
    handlers = handlers + [mbconv]
    from corvus.filereader import filereader
    handlers = handlers + [filereader]
    from corvus.helper import helper
    handlers = handlers + [helper]
    # The following are pure python handlers. Use only if module import throws no error
    # import only if module lmfit exists (fit dependency). Should probably
    # do this with numpy and scipy as well. 
    try:
        import lmfit
        from corvus.fit import fit
        handlers = handlers + [fit]
    except ImportError:
        print("Warning: lmfit not found. fit handler will be disabled.")
        pass
    
    try:
        import pymatgen
        from corvus.PyMatGen import PyMatGen
        handlers = handlers + [PyMatGen]
    except ImportError:
        print("Warning: pymatgen not found. PyMatGen handler will be disabled.")
        pass



    return handlers

# Commenting out Vasp for now since we have no content in the manual for it for
# now. I will readd once I include Scott's Handler.
#   return [Feff, FeffRixs, Dmdw, Abinit, Vasp, Nwchem, Orca]

def configure(config):
    from configparser import RawConfigParser
    from platform import system
    from pathlib import Path
    
    # Store path to corvus module
    # J. Kas - Now we are going to store things in home directories under the .Corvus folder.
    #config['bin'] = os.path.dirname(os.path.abspath(__file__))
    config['bin'] = os.path.join(str(Path.home()),'.Corvus')

    # Store path to current working directory (where user is running Corvus)
    config['cwd'] = os.getcwd()

    # Load Corvus defaults
    rcp = RawConfigParser(allow_no_value=True)
    #configFile = os.path.join(os.path.dirname(config['bin']), 'corvus','config')
    configFile = os.path.join(config['bin'], 'corvus.conf')


    # JJK: Make configuration file if it isn't already there. To start we will just 
    # search for feff10. Move definition of execs up here so that we can search for 
    # all of them.
    Execs_Dict = { 'feff'  :['atomic','compton','crpa','dmdw','eels',
                             'ff2x','fms','genfmt','ldos',
                             'mkgtr','opconsat','path','pot','rdinp','rhorrp',
                             'rixs','screen','sfconv','xsph'],
                   'dmdw'  :['dmdw'],
                   'orca'  :['orca'],
                   'abinit':['abinit','anaddb','mrgddb','mrggkk'],
                   'vasp'  :['vasp_gam','vasp_std'],
                   'nwchem':['nwchem'], 
                   'siesta':['siesta'],
                   'phsf'  :['phsf'],
                   'ocean' :['ocean.pl'],
                   'cif2cell':['cif2cell']}

    # If corvus.config doen't exist, try to find executables and write
    # config file.
    if not Path(configFile).is_file():
        print('\n\n\n')
        print('  Configuration file not found. Will search for feff10. ')
        # Try to find feff10
        if system() == 'Linux':
            search_paths=[Path.home() / Path('JFEFF_FEFF10/feff10/linux/')]
            ext=''
        elif system() == 'Darwin':
            search_paths=[Path.home() / Path('jfeff10.app/Contents/Resources/JFEFF/feff10/mac/')]
            ext=''
        elif system() == 'Windows':
            search_paths=[Path.home() / Path('JFEFF_FEFF10/feff10/bin/'),Path('C:/Program Files (x86)/JFEFF_FEFF10/feff10/bin/'),Path('C:/Program Files/JFEFF_FEFF10/feff10/bin/')]
            ext='.exe'
        else:
            print('')
            print('')
            print('')
            print('')
            print('')
            print('')
            print('############################################################################')
            print('############################################################################')
            print('##              Unknown operating system: Exiting                         ##')
            print('############################################################################')
            print('############################################################################')
            sys.exit()
            
        print('  System:', system())
        # Check if rdinp exists in any of the search paths
        feff_found = False

        for pth in search_paths:
            print('  Checking path:', pth)
            for exe in Execs_Dict['feff']:
                checkFile = pth / Path(exe + ext)
                if checkFile.is_file():
                    allExecs = True
                else:
                    allExecs = False
                    break
            if allExecs:
                feff_found = True
                feff_path = pth
                break

        if not feff_found:
            print('Writing empty corvus.conf')
            f = open(configFile, "w")
            
            f.write("[Executables]\n")
            f.write("dmdw     :\n")
            f.write("feff     :path_to_feff\n")
            f.write("abinit   :\n")
            f.write("nwchem   :\n")
            f.write("orca     :\n")
            f.write("gaussian :\n")
            f.write("vasp     :\n")
            f.write("siesta   :\n")
            f.write("ocean    :\n")
            f.write("cif2cell :\n")
            f.write("phsf     :\n")
            f.write(" \n")
            f.write("[Defaults]\n")
            f.write("prefix      : Corvus\n")
            f.write("inputsuffix : .inp\n")
            f.write("savesuffix  : .nest\n")
            f.write("checkpoints : off\n")
            f.write("parallelrun :\n")
            
            f.close()
            print('')
            print('')
            print('')
            print('')
            print('')
            print('')
            print('############################################################################')
            print('############################################################################')
            print('##    Feff10 was not found on this system.                                ##')
            print('##    Please register for feff10 and install it                           ##')
            print('##    before installing corvus.                                           ##')
            print('##                                                                        ##')
            print('##    Alternatively, edit the file                                        ##')
            print('##      ', configFile, '                                             ##')
            print('##    and insert the correct path to feff10.                              ##')
            print('############################################################################')
            print('############################################################################')
            sys.exit()

        else:
            # Feff was found. Write the config file
            print('Writing corvus.conf')
            f = open(configFile, "w")
            
            f.write("[Executables]\n")
            f.write("dmdw     :\n")
            f.write("feff     : " + str(feff_path) + "\n")
            f.write("abinit   :\n")
            f.write("nwchem   :\n")
            f.write("orca     :\n")
            f.write("gaussian :\n")
            f.write("vasp     :\n")
            f.write("siesta   :\n")
            f.write("ocean    :\n")
            f.write("cif2cell :\n")
            f.write("phsf     :\n")
            f.write(" \n")
            f.write("[Defaults]\n")
            f.write("prefix      : Corvus\n")
            f.write("inputsuffix : .inp\n")
            f.write("savesuffix  : .nest\n")
            f.write("checkpoints : off\n")
            f.write("parallelrun :\n")
            
            f.close()


    if configFile not in rcp.read(configFile):
        printAndExit('Error reading corvus.conf')
    config['pathprefix'] = rcp.get('Defaults', 'prefix')
    config['inputsuffix'] = rcp.get('Defaults', 'inputsuffix')
    config['savesuffix'] = rcp.get('Defaults', 'savesuffix')
    config['inputFile'] = config['pathprefix'] + config['inputsuffix']
    config['saveFile'] = config['pathprefix'] + config['savesuffix']
    config['checkpoints'] = rcp.getboolean('Defaults', 'checkpoints') 
    config['parallelRun'] = rcp.get('Defaults', 'parallelrun')
    utilpath = os.path.join(config['bin'], 'corvutils')
# Modified by FDV
    config['parsnipConf'] = os.path.join(utilpath, 'parsnip.corvus.config')
    config['parsnipForm'] = os.path.join(utilpath, 'parsnip.corvus.formats')
    
# Here we add a check to see if we are running under Cygwin. In Cygwin
# executables have a ".exe" extension. We also define all the different
# executables in one spot to make things more clear, and generalize the code.

    #Execs_Dict = { 'feff'  :['atomic','compton','crpa','dmdw','eels',
    #                         'ff2x','fms','genfmt','ldos',
    #                         'mkgtr','opconsat','path','pot','rdinp','rhorrp',
    #                         'rixs','screen','sfconv','xsph'],
    #               'dmdw'  :['dmdw'],
    #               'orca'  :['orca'],
    #               'abinit':['abinit','anaddb','mrgddb','mrggkk'],
    #               'vasp'  :['vasp_gam','vasp_std'],
    #               'nwchem':['nwchem'], 
    #               'siesta':['siesta'],
    #               'phsf'  :['phsf'],
    #               'ocean' :['ocean.pl'],
    #               'cif2cell':['cif2cell']}


    exe_ext = ''
    if system()[0:6]=='CYGWIN':
      exe_ext = '.exe'

    for code in list(Execs_Dict.keys()):
      path2code = rcp.get('Executables', code)
      for cmd in Execs_Dict[code]:
        config[cmd] = which(cmd+exe_ext, path2code)

# NOTE FDV: Here I try to add code to make JJK's version work. We need the path
# to feff to be avilable in config.
    path2feff = rcp.get('Executables', 'feff')
    config['feff'] = path2feff

    path2siesta = rcp.get('Executables', 'siesta') 
    path2cif2cell = rcp.get('Executables', 'cif2cell') 
    path2ocean = rcp.get('Executables', 'ocean') 
    config['siesta'] = path2siesta
    config['cif2cell'] = path2cif2cell
    config['ocean'] = path2ocean
    path2phsf = rcp.get('Executables', 'phsf')
    config['phsf'] = path2phsf
# Initialize system with user input
def initializeSystem(config, system, doc):
    import corvutils.parsnip

    conf    = checkFile(config['parsnipConf'])
    inp     = config['inputFile']

# Debug: FDV
#print parsnip.parse(conf, inp, mode=['read','UseDefaults'])
#   sys.exit()

    if os.path.exists(inp):
        system.update(corvutils.parsnip.parse(conf, inp, mode=['read','UseDefaults']))


    with open(conf, 'r',encoding="utf8") as conf_file:
        doc2=corvutils.parsnip.readConfig_for_Doc(conf_file).asDict()
    
    doc2 =  {k.lower(): v for k, v in doc2.items()}
    doc.update(doc2)
    

# Debug: FDV
#   pp_debug.pprint(system)
#   sys.exit()

# Basic Workflow Generator
# JK - added system as input to generateWorkflow to accomodate handlers that call other 
# handlers or workflows such as "fit" or "average".
# JK - added config as input to allow use of paths to check for necessary executables
# before adding handlers to available handler list.
def generateWorkflow(target, handlers, system, config, desc=''):
    from corvus.structures import Workflow

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
      if handler_name in list(availableHandlers_map.keys()):
        useHandlers.append(availableHandlers_map[handler_name])
        #print useHandlers
      else:
        print(("Handler %s not in list of available handlers:" % (handler_name)))
        for s in list(availableHandlers_map.keys()):
          print(("%s" % s))
        #sys.exit()
        exitOneshot()
        return

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
    # JK - For now, this assumes that the list of targets supplied by the user is
    # an ordered list that should be run sequentially. We should really have some way to input
    # a list of lists of targets (we can actually do that easily) so that targets in a
    # given list are to be run in parallel, while separate lists will be run sequentially.
    # This gives the user easy access to simple workflows, while still allowing the power 
    # of the underlying corvus machinery. 

    # JK - Changing the operation of user defined handler list. Now we will prefer user defined handlers, but not
    # require them. If a target cannot be produced by any of the user defined handlers, search for one that can produce
    # that target. 
    # First check what handlers can produce each target. There might be more than one for a 
    # given target.
    for t in reversed(mainTargets):
       targets = set([t])
       while len(targets) > 0:
           noMatch = True

           # Loop through all handlers that can produce target t and find all handlers that can produce the given target
           #possibleHandlers = [ h for h
           possibleHandlers = []
           for target in targets:
               found = False
               # Loop through user supplied handlers first to see if we can find one that produces this target.
               for h in useHandlers:
                   if h.canProduce(target):
                       possibleHandlers = possibleHandlers + [h]
                       found = True
                       break # For now just use the first handler available
               if not found:
                   # loop through all available handlers now since we didn't find one in the user
                   # supplied list. For now, use the first handler available. In future, we might
                   # want to have some other method of chosing which handler to use based on input
                   # or cost etc. 
                   for h in availHandlers:
                        if h.canProduce(target):
                            possibleHandlers = possibleHandlers + [h]
                            break # For now just use the first handler available
                   
           for h in possibleHandlers:
               htargets = [target for target in targets if h.canProduce(target)]
               if htargets:
                   noMatch = False
                   # JK - The following loops over all possible subsets of the list of targets
                   # from largest subset to smallest, to try to find the smallest number
                   # of handlers to handle all of the target properties. Note: CanProduce only
                   # checks if a property is implemented, even though other properties are specified
                   # as output of an implemented property. 
                   for subset in reversed(sorted(subs(htargets), key=len)):
                       l = list(subset)
                       if h.canProduce(l):
                           workflow.addExchangeAt(0, h.sequenceFor(l,system))
                           targets.difference_update(subset)  
                           targets.update(set(workflow.getRequiredInput()))
                           #allTargets=allTargets | targets 
                           break
               #print(targets)
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
    elif isinstance(data, list) and len(data) > 1:
        # We have columns of data
        if isinstance(data[0], list) and len(data[0]) == len(data[1]):
            for row in zip(*data):
                f.write('    '.join(map(str,row)) + '\n')
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

def exitOneshot():
    if exitOneshot.sys_exit:
        sys.exit()
    else:
        return

exitOneshot.sys_exit = True

def setExit(sys_exit=True):
    exitOneshot.sys_exit=sys_exit

def printAndExit(msg):
    sys.stderr.write(msg + '\n')
    sys.stderr.flush()
    exitOneshot()
    return

def generateAndRunWorkflow(config, system, targetList):
    import copy
    #if 'target_list' in system.keys():
    #  targetList = system['target_list']
    #else:
    #  print 'Provide target properties or Workflow'
    #  sys.exit()

    if 'usehandlers' in list(system.keys()):
      handlerList = system['usehandlers'][0]
    else:
      #print('Provide the handler list to be used in this calculation')
      #exitOneshot()
      #return
      handlerList = []

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

# JK change to use of argparse in place of passing argv
def oneshot(): 
#def oneshot(argv,sys_exit=True):
    import getopt, pickle, argparse

    # JJK - argv used to be passed to oneshot. Now we get the arguements
    #       directly here. Should change so that oneshot can be called
    #       as a function as well, but later. 
    argv = sys.argv
    setExit(True) # JJK - No arguments should be passed to oneshot,
                  # so setting sys_exit to True here.

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
    shortopts = 'vcrt:i:w:s:j:h:'
    longopts = ['target=','input=','workflow=','checkpoints','resume','save=',
               'jump=','prefix=','parallelrun=','help','version']
    try:
        opts, args = getopt.getopt(argv[1:], shortopts, longopts)
    except getopt.GetoptError as err:
        print((str(err)))
        usage()
        exitOneshot()
        return

    ### Read user options ###
    resume = checkpoints = False
    jumpPoint = 0
# Modified by FDV:
# Removing the target option from the cli. From now on we do it through the
# input.
#   inputFile = saveFile = workflowFile = targetList = pathPrefix = parallelRun = None
    helpOnly = False
    inputFile = saveFile = workflowFile = pathPrefix = parallelRun = None
    for opt, arg in opts:
        #if opt in ('-t', '--target'):
        #    targetList = arg.split(',')
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
        elif opt in ('-v', '--version'):
            print_version()
            sys.exit()
        elif opt in ('-h', '--help'):
            # Print workflow and requirements for this target and exit.
            helpOnly = True
            targetList = [arg.split(',')]
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
    doc = {}
    workflowStart = 0

    # Read in saved state if restarting calculation midway
    if resume and os.path.exists(config['saveFile']):
        saveState = pickle.load(open(config['saveFile'], 'r'))
        system.update(saveState['system'])
        workflow = saveState['workflow']
        workflowStart = saveState['index']

    # Update System with any user input
    initializeSystem(config, system, doc)
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
    if 'target_list' in list(system.keys()):
      targetList = system['target_list']
    elif not helpOnly:
      print('Provide target properties or Workflow')
      exitOneshot()
      return
   

# Added by FDV
# At this point we set the handler list to make the workflow generator work more
# generally.
    if 'usehandlers' in list(system.keys()):
      handlerList = system['usehandlers'][0]
    else:
      #print('Provide the handler list to be used in this calculation')
      #exitOneshot()
      #return
      handlerList = []

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
        if len(workflow.getRequiredInput()) > 1:
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

    # JJK add help by keyword. Print help for keyword and also for all requirements.
    if helpOnly:
        # print doc for each missing required input
        required = workflow.getRequiredInput()
        if len(required) > 1:
            for target in targetList[0]:
               if target.lower() not in doc:
                  print('Help for requested propery ' +  target + '.')
                  print('Required input is as follows:')
            for key in required:
               req = key.lower()
               if req in doc:
                  print('   ' + req + ':')
                  for i,line in enumerate(doc[req]):
                      if line == '%':
                          print('      ' + doc[req][i+1])
        elif len(required) == 1:
            if required[0].lower() == 'all':
                for key in doc:
                    print('\n\n   ' + key + ':')
                    for i,line in enumerate(doc[key]):
                        if line == '%':
                            print('      ' + doc[key][i+1])
            elif required[0].lower() in doc:
                print('\n\n   ' + required[0] + ':')
                for i,line in enumerate(doc[required[0].lower()]):
                    if line == '%':
                            print('      ' + doc[required[0].lower()][i+1])
            else:
                noHelp = True
                for key in doc:
                    if key.startswith(required[0].lower()):
                        noHelp = False
                        print('\n\n   ' + key + ':')
                        for i,line in enumerate(doc[key]):
                            if line == '%':
                                print('      ' + doc[key][i+1])

                if noHelp: 
                    print('No help found for keywords: ' + targetList[0][0]) 

        sys.exit()
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
    if jumpPoint != 0:
        workflowStart = jumpPoint - 1

    # Save intial state
    saveState = {'system':system,'workflow':workflow}
    # Run Workflow
    i = max(0, workflowStart)
    #print workflow.sequence
    #print('Length of sequence', len(workflow.sequence))
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
    with open(config['saveFile'], 'wb') as saveFile:
        pickle.dump(saveState, saveFile, pickle.HIGHEST_PROTOCOL)

    # Print output
    for token in workflow.target:
        if system[token] is not None:
            prettyprint(token, config, system)
        else: 
            printAndExit('Error: target [' + token + '] not produced')

def print_version():
    print('#####################################################')
    print('#                                                   #')
    print('#          Corvus version 1.0.9')
    print('#                                                   #')
    print('#####################################################')
