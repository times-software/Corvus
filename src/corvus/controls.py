import sys, os

# Debug: FDV
import pprint
pp_debug = pprint.PrettyPrinter(indent=4)

# Define a few things to make the output a bit cleaner
# NOTE FDV: Later I will put this in a module that can be loaded where needed
#           so we can keep it consisten over the whole code.
Blank_line = '\n'
Separator1_char = '='
Separator2_char = '-'
Page_width = 80
Separator1 = Page_width*Separator1_char
Separator2 = Page_width*Separator2_char

# Define the available handlers by hand here
# Note FDV: There should be a way to do this automatically, but can't think of
# one right now
def availableHandlers():
    from abinit import Abinit
    from dmdw import Dmdw
    from feff import Feff
    from feffrixs import FeffRixs
    from vasp import Vasp
    from nwchem import Nwchem
    from orca import Orca

# Commenting out Vasp for now since we have no content in the manual for it for
# now. I will readd once I include Scott's Handler.
#   return [Feff, FeffRixs, Dmdw, Abinit, Vasp, Nwchem, Orca]
    return [Feff, FeffRixs, Dmdw, Abinit, Nwchem, Orca]

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
    config['debug']   = int(rcp.get('Defaults','debug'))
    config['verbose'] = int(rcp.get('Defaults','verbose'))
    
# Here we add a check to see if we are running under Cygwin. In Cygwin
# executables have a ".exe" extension. We also define all the different
# executables in one spot to make things more clear, and generalize the code.

    Execs_Dict = { 'feff'  :['atomic','band','compton','crpa','dmdw','eels',
                             'ff2x','fms','fullspectrum','genfmt','ldos',
                             'mkgtr','opconsat','path','pot','rdinp','rhorrp',
                             'rixs','screen','sfconv','xsph'],
                   'dmdw'  :['dmdw_standalone'],
                   'orca'  :['orca'],
                   'abinit':['abinit','anaddb','mrgddb','mrggkk'],
                   'vasp'  :['vasp_gam','vasp_std'],
                   'nwchem':['nwchem'] }


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
def generateWorkflow(config,target, handlers, desc=''):
    from structures import Workflow

# NOTE FDV:
# At some point we will have to put some printouts during the generation of the
# workflow, to make it more inforamtive and easier to troubleshoot. For now I
# just simply added the config variable to the interface of the function so we
# can control output later.

# From the handlers list we generate a mapping dictionary to identify the 
# handlers requested by the user
    availableHandlers_map = { h.__name__:h for h in availableHandlers() }
# Debug: FDV
#   print 'availableHandlers_map = '
#   pp_debug.pprint(availableHandlers_map)
#   sys.exit()

# Now we create the useHandlers list from the input and the available ones
    useHandlers = []
# Debug: FDV
#   print 'handlers = ', handlers
#   print 'availableHandlers_map.keys = ', availableHandlers_map.keys()
#   sys.exit()
    for handler_name in handlers:
      if handler_name in availableHandlers_map.keys():
        useHandlers.append(availableHandlers_map[handler_name])
# Debug: FDV
#       print useHandlers
      else:
        print("Handler %s not in list of available handlers:" % (handler_name))
        for s in availableHandlers_map.keys():
          print("%s" % s)
        sys.exit()
# Debug: FDV
#   sys.exit()

# Debuf: FDV
#   pp_debug.pprint(availableHandlers)
#   pp_debug.pprint(useHandlers)
#   sys.exit()
    workflow = Workflow(target, desc=desc)
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
## Debug: FDV
##   berp = {}
##workflow print statements uncommented, Krsna
#    print type(workflow.sequence[0])
#    print workflow.sequence[0].handler
##   workflow.sequence[0]({1:'a'},{1:'a'},berp)
##   sys.exit()
#-----------------------------------------------------------------------------

    maintargets = workflow.getRequiredInput()
    #print "Done printing ri" # Debug JJK
    for t in maintargets:
       targets = set([t])
       while len(targets) > 0:
           noMatch = True
# Modified by FDV
# NOTE: Can't believe I missed this. 
#          for h in availableHandlers():
           for h in useHandlers:
               #print "Checking Handler: " # Debug JJK
               #print h # Debug JJK
               htargets = [t for t in targets if h.canProduce(t)]
               #print htargets # Debug JJK
               #print "Done checking handler." # Debug JJK
               if htargets:
                   noMatch = False
                   for subset in reversed(sorted(subs(htargets), key=len)):
                       l = list(subset)
                       if h.canProduce(l):
                           workflow.addExchangeAt(0, h.sequenceFor(l))
                           # h.setDefaults(input,h)
                           targets.difference_update(subset)
                           targets.update(set(workflow.getRequiredInput()))
                           break

           if noMatch:
               break

    return workflow

# Formats output strings
def prettyprint(token, config, system):
    import pprint
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
    print "Usage: corvus.py [options]"
    print "       Possible options:"
#uncommented print line below, Krsna
# Commented again by FDV: We do not acquire targets though the CLI anymore.
#   print "         -t, --target      [comma-separated list]"
    print "         -h, --help"
    print "         -a, --avail"
    print "         -v, --verbose     [verbose level: 0, 1, ...]"
    print "         -d, --debug       [debug level: 0, 1, ...]"
    print "         -w, --workflow    [filename]"
    print "         -i, --input       [filename]"  
    print "         -c, --checkpoints"  
    print "         -r, --resume"  
    print "         -s, --save        [filename]"  
    print "             --prefix      [alphanumeric string (-_. okay)]"
    print "             --parallelrun [command string]"

def checkFile(filename):
    if not os.path.exists(filename):
        printAndExit('File does not exist: ' + filename)
    else:
        return filename

def printAndExit(msg):
    sys.stderr.write(msg + '\n')
    sys.stderr.flush()
    sys.exit()

def Available_Handlers_and_Targets():

  import pprint
  pp_debug = pprint.PrettyPrinter(indent=4)

  for h in availableHandlers():
    print 'Handler: {0}'.format(h.__name__)
    Implemented = h.Produces()
#   key_len = max([ len(key) for key in Implemented.keys() ])
#   fmt = '  {:' + str(key_len) + 's} <-'
    fmt = '  {:s} <-'
    for prop in sorted(Implemented.keys()):
      strng = '  ' + prop + ' <-'
      print strng,
      first = 0
      for req in Implemented[prop]['req']:
        print first*(len(strng)+1)*' ' + req
        first = 1
      print ''
    print 64*'-' 

# Handles command line arguments and runs through workflow
def oneshot(argv):
    import getopt, pickle

# Modified by FDV:
# Removing the target option from the cli. From now on we do it through the
# input.
#   shortopts = 'crt:i:w:s:j:'
#   longopts = ['target=','input=','workflow=','checkpoints','resume','save=',
#               'jump=','prefix=','parallelrun=']
    shortopts = 'v:d:ahcr:i:w:s:j:'
    longopts = ['verbose','debug=','avail','help','input=','workflow=','checkpoints','resume','save=',
                'jump=','prefix=','parallelrun=']
    try:
        opts, args = getopt.getopt(argv[1:], shortopts, longopts)
    except getopt.GetoptError as err:
        print str(err)
        usage()
        sys.exit(2)

# Debug
#   print opts
#   sys.exit()

    ### Read user options ###
    Verbose_Opt = {'opt':False,'val':0}
    Debug_Opt   = {'opt':False,'val':0}
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
        elif opt in ('-a', '--avail'):
            Available_Handlers_and_Targets()
            sys.exit()
        elif opt in ('-v', '--verbose'):
            if not arg.isdigit():
                printAndExit('Verbose argument should be a positive integer.')
            Verbose_Opt = {'opt':True,'val':int(arg)}
        elif opt in ('-d', '--debug'):
            if not arg.isdigit():
                printAndExit('Debug argument should be a positive integer.')
            Debug_Opt = {'opt':True,'val':int(arg)}
        elif opt in ('-h', '--help'):
            usage()
            sys.exit()
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
    if Verbose_Opt['opt']:
        config['verbose'] = Verbose_Opt['val']
    if Debug_Opt['opt']:
        config['debug'] = Debug_Opt['val'] 
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

# Added by FDV:
# Improving the output a bit
    if config['verbose'] > 0:
      print_header()

# Print the configuration state to be used in this run
    if config['debug'] > 0:
      print ''
      print 'Current config dictionary:'
      pp_debug.pprint(config)
      print ''

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
#DASb
    if config['debug'] == 10:
      print 'system = ', system
#DASe

# Added by FDV
# At this point we set the target list based on the content of the input file,
# rather than the command line.
    if 'target_list' in system.keys():
      targetList = system['target_list']
    else:
      print 'Provide target properties or Workflow'
      sys.exit()

# Added by FDV
# At this point we set the handler list to make the workflow generator work more
# generally.
    if 'usehandlers' in system.keys():
      handlerList = system['usehandlers'][0]
    else:
      print 'Provide the handler list to be used in this calculation'
      sys.exit()

# Here we should probably add a check to see if the target list if empty, but
# I don't think it can happen. Leaving for future

    if config['verbose'] > 0:
      print 'List of requested targets:'
      print targetList
      print ''
      print 'List of requested handlers:'
      print handlerList
      print ''

    # Set up Workflow or read from file
    if workflowFile is not None:
        # Read Workflow from file
        workflow = pickle.load(open(workflowFile, 'r'))
    elif not resume:
        # Create Workflow
        autodesc = 'Calculate ' + ', '.join(targetList[0])
        workflow = generateWorkflow(config,targetList, handlerList, desc=autodesc)

    if config['verbose'] > 0:
      print 'Workflow for this calculation:'
      print workflow
      print ''

    # Check for any missing user input
    required = set(workflow.getRequiredInput())

# NOTE FDV: This is already printed at the end of workflow above.
#   if config['debug'] > 0:
#     print 'List of required user inputs:'
#     print required
#     print ''

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

    if config['verbose'] > 0:
      print "{0}".format(Separator1)
      print "{0}".format('Starting Workflow Sequence'.center(Page_width))
      print "{0}".format(Separator2)
#     print "{0}".format(1*Blank_line),

    # Run Workflow
    i = max(0, workflowStart)
    while i < len(workflow.sequence):
        if config['checkpoints']:
            saveState['system'] = system
            saveState['index'] = i
            with open(config['saveFile'], 'w') as saveFile:
                pickle.dump(saveState, saveFile, pickle.HIGHEST_PROTOCOL)
        
        config['xcIndex'] = i + 1
        if config['verbose'] > 0:
          print 'Workflow Step [{0}]:'.format(i+1)
        workflow.sequence[i].go(config, system)
        if config['verbose'] > 0:
          if i == len(workflow.sequence)-1:
            print "{0}".format(Separator1)
          else:
            print "{0}".format(Separator2)
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

def print_header():

  from pkg_resources import get_distribution
  from corvutils import parsnip as pnip
  from datetime import datetime

# For now, we get a minimal amount of information from the Corvus package setup
  corvus_pkg = get_distribution('corvus')
  corvus_metadata = corvus_pkg.get_metadata('PKG-INFO').encode('ascii','ignore')
  corvus_info = pnip.parseMetaData(corvus_metadata)

# Print the cover page:
  print "{0}".format(Separator1) 
  print "{0}".format(1*Blank_line),
  print "{0}".format(corvus_info['Name'].capitalize().center(Page_width))
  print "{0}".format(1*Blank_line),
  print "{0}".format(corvus_info['Summary'].center(Page_width))
  print "{0}".format(1*Blank_line),
  print "{0}".format(('Version: '+corvus_info['Version']).capitalize().center(Page_width))
  print "{0}".format(1*Blank_line),
  print "{0}".format(('Web: '+corvus_info['Home-page']).capitalize().center(Page_width))
  print "{0}".format(('E-Mail: '+corvus_info['Author-email']).capitalize().center(Page_width))
  print "{0}".format(1*Blank_line),
  print "{0}".format(corvus_info['Author'].center(Page_width))
  print "{0}".format(1*Blank_line),
  print "{0}".format(Separator1)
  print "{0}".format(1*Blank_line),

