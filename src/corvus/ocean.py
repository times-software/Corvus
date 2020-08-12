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

implemented['OceanXANES'] = {'type':'Exchange','out':['OceanXANES'],'cost':3,
                        'req':['ocean.edges', 'ocean.ecut', 
                            'ocean.pp_list', 'ocean.diemac', 'ocean.xred', 'ocean.typat', 'ocean.natom', 
                            'ocean.znucl', 'ocean.ntypat', 'ocean.rprim', 'ocean.acell', 'ocean.photon.operator',
                            'ocean.photon.polarization'],'desc':'Calculate XANES using ocean.'}



class Ocean(Handler):
    def __str__(self):
        return 'Ocean Handler'

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
            unresolved = {o for o in output if not Ocean.canProduce(o)}
            canProduce = (o for o in output if Ocean.canProduce(o))
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
            return Exchange(Ocean, f('req'), f('out'), cost=f('cost'), desc=f('desc'))

    @staticmethod
    def prep(config):
        subdir = config['pathprefix'] + str(config['xcIndex']) + '_OCEAN'
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
        
        # set atoms and potentials

        # Set directory to ocean executables.
        # Debug: FDV
        #       pp_debug.pprint(config)
        # Debug: FDV
        #       sys.exit()
        dir = config['xcDir']

        # Copy ocean related input to oceanInput here. Later we will be overriding some settings,
        # so we want to keep the original input intact.
        oceanInput = {key:input[key] for key in input if (key.startswith('ocean.') and 'photon.' not in key)}
        photonInput = {key:input[key] for key in input if (key.startswith('ocean.') and 'photon.' in key)}
        photonfile = os.path.join(dir, 'photon1')
        # Write photon1 input file
        lines = []
        lines.append(str(photonInput['ocean.photon.operator'][0][0]))
        lines.append('cartesian ' + ' '.join([str(value) for value in photonInput['ocean.photon.polarization'][0]]))
        lines.append('end')
        if 'ocean.photon.qhat' in photonInput:
            lines.append('cartesian ' + ' '.join([str(value) for value in photonInput['ocean.photon.qhat'][0]]))
        lines.append('end')
        lines.append(str(photonInput.get('ocean.photon.energy',0.0)[0][0]))
        writeList(lines,photonfile)

        #writePhotonInput(photonInput,dir)

        # Generate any data that is needed from generic input and populate oceanInput with
        # global data (needed for all feff runs.)

        # Set directory for this exchange
        oceandir = config['ocean']

        
        # Set input file
        inpf = os.path.join(dir, 'ocean.in')

        # Loop over targets in output. Not sure if there will ever be more than one output target here.
        for target in output:
            if (target == 'OceanXANES'):               
                # Set output and error files
                with open(os.path.join(dir, 'corvus.OCEAN.stdout'), 'w') as out, open(os.path.join(dir, 'corvus.OCEAN.stderr'), 'w') as err:
                    # Get pseudopotentials
                     
                    # Write input file for Ocean.
                    writeXANESInput(oceanInput,inpf)
                    
                    # Copy necessary files to dir
                    # If program is set to hamann, or not set, we are using
                    # John Vinson's version of oncvpsp with q-e
                    # Don't need fhi files, need to append UPF to
                    # files listed in pp_list, and also copy oncvpsp
                    # input file for absorber. In this case we don't 
                    # need .opts or .fill files.
                    if oceanInput.get('ocean.opf.program',[['hamann']])[0][0] == 'hamann':
                        # By default, use q-e? However, this should depend on q-e being available?
                        for file in oceanInput['ocean.pp_list']:
                            fileUPF = file[0] + '.UPF'
                            fileONCV = file[0] + '.in' # Should this always be named .in?
                            # UPF file should exist, .in file should only exist for absorber. For now, check that .in
                            # file exists for at least one of pp_list, and error otherwise.
                            shutil.copy(fileUPF,dir)
                            if os.path.exists(fileONCV):
                                shutil.copy(fileONCV,dir)
                    else:
                        # This should work for abinit or q-e.
                        for file in os.listdir("."):
                            if file.endswith(".UPF"):
                                shutil.copy(file,dir)

                        for file in oceanInput['ocean.pp_list']:
                            shutil.copy(file[0],dir)

                        if 'ocean.opf.fill' in oceanInput:
                            for fill in oceanInput['ocean.opf.fill']:
                                shutil.copy(fill[1],dir)

                        if 'ocean.opf.opts' in oceanInput:
                            for opt in oceanInput['ocean.opf.opts']:
                                shutil.copy(opt[1],dir)


                    # Loop over executables: This is specific to feff. Other codes
                    # will more likely have only one executable.
                    executables = ['ocean.pl']
                    args=['ocean.in']
                    iExec = 0
                    for executable in executables:
                        runExecutable(oceandir,dir,executable,args,out,err)

                
                # For now, I am only passing the directory.
                print('Setting output')
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
        lines = lines + getInpLines(input,key) 
        
    
    # Print ocean input file
    writeList(lines, inpfile)
    
def getInpLines(input,token):
    lines=[]
    block=False
    endblock=' '
    key = token[len('ocean.'):]
    if token in input:
        for element in input[token]: # Takes care of single and multi-line input.
            lines.append(' '.join([str(value) for value in element]))

        if len(input[token]) > 1:
            lines.insert(0,key.upper() + '{')
            endblock = '}'
        else:                                     # Most have arguments on the same line as keyword.
            lines[0] = key.upper() + '{ ' + lines[0] + ' }'


    # Add a blank line after each line
    lines.append(endblock)
    lines.append('')

    return lines

def writeList(lines, filename):
    with open(filename, 'w') as f:
        f.write('\n'.join(lines))

def runExecutable(execDir,workDir,executable, args,out,err):
    # Runs executable located in execDir from working directory workDir.
    # Tees stdout to file out in real-time, and stderr to file err.
    print('Running exectuable: ' + executable)
    # Modified by FDV:
    # Adding the / to make the config more generic
    # Modified by JJK to use os.path.join (even safer than above).
    execList = [os.path.join(execDir,executable)] + args
    print execList
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
            print('Error in executable: ' + executable)
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

def writeXANESInput(input, oceaninp='ocean.in'):
    
    lines = []
    #setInput(input,'feff.print',[[5,0,0,0,0,0]],Force=True)

    writeInput(input,oceaninp)


