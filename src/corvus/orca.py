from structures import Handler, Exchange, Loop, Update
import corvutils.pyparsing as pp
import os, sys, subprocess

# Debug: FDV
import pprint
pp_debug = pprint.PrettyPrinter(indent=4)

# Define dictionary of implemented calculations
implemented = {}
strlistkey = lambda L:','.join(sorted(L))
implemented['OptimizedStructure'] = {'type':'Exchange','out':['OptimizedStructure','cluster'],'cost':2,
                        'req':['cluster'],'desc':'Optimize structure using ORCA.'}
#implemented['cluster'] = {'type':'Exchange','out':['OptimizedStructure','cluster'],'cost':2,
#                        'req':['cluster'],'desc':'Optimize structure using ORCA.'}

# Added by FDV:
# Define the default name for the orca input file
Orca_In_File_Name = 'orca.in'

class Orca(Handler):
    def __str__(self):
        return 'ORCA Handler'

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
            unresolved = {o for o in output if not Orca.canProduce(o)}
            canProduce = (o for o in output if Orca.canProduce(o))
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
            raise LookupError('Corvus cannot currently produce ' + key + ' using ORCA')
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
            raise LookupError('Corvus cannot currently produce ' + key + ' using ORCA')
        f = lambda subkey : implemented[key][subkey]
        if f('type') is 'Exchange':
            return Exchange(Orca, f('req'), f('out'), cost=f('cost'), desc=f('desc'))

    @staticmethod
    def prep(config):
        subdir = config['pathprefix'] + str(config['xcIndex']) + '_ORCA'
        xcDir = os.path.join(config['cwd'], subdir)
        # Make new output directory if if doesn't exist
        if not os.path.exists(xcDir):
            os.mkdir(xcDir)
        # Store current Exchange directory in configuration
        config['xcDir'] = xcDir



    @staticmethod
    def run(config, input, output):
        dir = config['xcDir']
        out = open(os.path.join(dir, 'corvus.ORCA.stdout'), 'w')
        err = open(os.path.join(dir, 'corvus.ORCA.stderr'), 'w')

        if set(output) >= {'OptimizedStructure'}:
            print 'Optimizing strucure using ORCA.' # Debug JJK
            inpf = os.path.join(dir, Orca_In_File_Name)
            writeOptimizeStructureInp(input, orcainp=inpf)
            files = inpf
                 
# Modified by FDV:
# Mimiking what Josh was doing with a script
#           out = open(os.path.join(dir, 'corvus.ORCA.stdout'), 'w')
            out = open(os.path.join(dir, 'orca.out'), 'w')
            err = open(os.path.join(dir, 'corvus.ORCA.stderr'), 'w')
# Debug: FDV
#           print 'config'
#           pp_debug.pprint(config)
#           sys.exit()
# Modified by FDV:
# Add the input file to the subprocess call
            executable = [ config['orca'], Orca_In_File_Name ]
#           p = subprocess.Popen([config['orca']], cwd=dir)
            p = subprocess.Popen(executable, cwd=dir, stdin=None, stdout=out, stderr=err)

            p.wait()
            out.close()
            err.close()

        translateOutput(config, input, output) 
    @staticmethod
    def cleanup(config):
        pass

##### Generic Helper Methods ##########

def translateOutput(config, input, output):
    dir = config['xcDir']
    # Run before loop to prevent possible duplicated call to function
    for target in output:
        if (target == 'OptimizedStructure'):
            output[target] = 'OptimizedStructure'
        elif (target == 'cluster'):
            filename = os.path.join(dir, 'orca.xyz')
            output[target] = readXYZ(filename)

def check(input, token, default=None):
    if token in input:
        return input[token]
    else:
        return default

def writeList(lines, filename):
    with open(filename, 'w') as f:
        f.write('\n'.join(lines))

def readXYZ(filename):
    with open(filename, 'r') as file:
        cleanStr = file.read()
    
    xyzList = [ s.split() for s in cleanStr.split("\n") if s.strip()][2:]
    xyzOut = []
    for xyzLine in xyzList:
       xyzOutLine = [ xyzLine[0], float(xyzLine[1]), float(xyzLine[2]), float(xyzLine[3]) ]
       xyzOut.append(xyzOutLine)
    return xyzOut

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

# Periodic Table information from ORCA 
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
   158.92534, 162.5,     164.93032, 167.26,    168.93421, 173.04,  174.967,   178.49,
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
assert len(atomicMasses) == len(atomicSymbols), "ORCA Handler: Mismatch in periodic table!"
nElem = len(atomicSymbols)
for i in xrange(nElem):
    num = i + 1
    sym = atomicSymbols[i]
    mass = atomicMasses[i]
    ptable[num] = {'mass':mass, 'symbol':sym, 'number':num}
    ptable[sym] = {'mass':mass, 'symbol':sym, 'number':num}

def writeOptimizeStructureInp(input, orcainp=Orca_In_File_Name):
    lines = []
# Modified by FDV:
# This is a temporary fix since I am not running orca in parallel yet
#   lines.append('! BP86 def2-TZVP TightSCF SlowConv Grid3 Opt PAL8')
    lines.append('! BP86 def2-TZVP TightSCF SlowConv Grid3 Opt')
    lines.append('')
# Modified by FDV:
# This is a temporary fix since I am not running orca in parallel yet
# Commenting out the next two lines
#   lines.append('%pal nprocs 16')
#   lines.append('end')
    lines.append('')

    # Get minimum multiplicity
    nElectron = -input['charge'][0][0]
    for atom in input['cluster']:
       nElectron = (nElectron + atomicSymbols.index(atom[0]) + 1)
    if nElectron % 2 == 0:
       multiplicity = ' 1'
    else:
       multiplicity = ' 2'
    
    lines.append('* xyz ' + str(input['charge'][0][0]) + multiplicity)

    for atom in input['cluster']:
        lines.append(''.join(''.join([str(e),' ']) for e in atom)) 
    lines.append('*')
    # Print orca input file
    writeList(lines, orcainp)
