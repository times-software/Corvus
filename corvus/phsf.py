from corvus.structures import Handler, Exchange, Loop, Update
import corvutils.pyparsing as pp
import os, sys, subprocess, shutil #, resource
import re
import numpy as np
import matplotlib.pyplot as plt
import pprint
# Debug: FDVi
import pprint
pp_debug = pprint.PrettyPrinter(indent=4)

# Define dictionary of implemented calculations
implemented = {}
strlistkey = lambda L:','.join(sorted(L))

implemented['spectralFunction'] = {'type':'Exchange','out':['spectralFunction'],'cost':1,
			'req':['siestaBetaofOmega'],'desc':'Calculates the spectral function'}

class phsf(Handler):
    def __str__(self):
        return 'PHSF Handler'

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
            unresolved = {o for o in output if not phsf.canProduce(o)}
            canProduce = (o for o in output if phsf.canProduce(o))
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
        if f('type') == 'Exchange':
            return Exchange(phsf, f('req'), f('out'), cost=f('cost'), desc=f('desc'))

    @staticmethod
    def prep(config):
        subdir = config['pathprefix'] + str(config['xcIndex']) + '_PHSF'
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

        # Set directory to phsf executables.
        # Debug: FDV
        #       pp_debug.pprint(config)
        phsfdir = config['phsf']
        # Debug: FDV
        #       sys.exit()

        # Copy phsf related input to phsfInput here. Later we will be overriding some settings,
        # so we want to keep the original input intact.
        phsfInput = {key:input[key] for key in input if key.startswith('phsf.')}

        # Generate any data that is needed from generic input and populate phsfInput with
        # global data (needed for all feff runs.)

        # Set directory for this exchange
        dir = config['xcDir']

        
        # Set input file
        inpf = os.path.join(dir, 'input.fdf')

        # Loop over targets in output. Not sure if there will ever be more than one output target here.
        for target in output:
            if (target == 'spectralFunction'):               
                print("reached spectralFucntion")
		        # Set output and error files
                with open(os.path.join(dir, 'corvus.PHSF.stdout'), 'w') as out, open(os.path.join(dir, 'corvus.PHSF.stderr'), 'w') as err:
                    
                    # Pull out Beta(omega) from the system dictionary
                    Beta = input['siestaBetaofOmega']
                    Beta = np.array(Beta).T
                    Betafile = os.path.join(dir, "Beta.dat")
                    np.savetxt(Betafile, Beta)
                    
                    # Pull out inputs for phsf code
                    ekeqp = str(input['phsf.ekeqp'][0][0]) + " " + str(input['phsf.ekeqp'][0][1])
                    numtimepoints = str(input['phsf.numtimepoints'][0][0])
                    broadening = str(input['phsf.broadening'][0][0])
                    broadtype = '0'

                    # Create input file phsf.in for phsf code
                    lines = ["Beta.dat", ekeqp, numtimepoints, broadening, broadtype]
                    phsffile = os.path.join(dir, "phsf.in")
                    writeList(lines, phsffile)
                    
                    # Run phsf code using input file that was made in previous step
                    executables = [os.path.join(phsfdir,'phsf')]                        
                    for executable in executables:
                        runExecutable('',dir,executable,[' '],out,err)

                
                # For now, I am only passing the directory.
                print('Setting output')
                outFile=os.path.join(dir,'jkexp.dat')
                output[target] = np.loadtxt(outFile).T.tolist()
                # print(output[target])
                # output[target] = dir
                


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
        
    
    # Print input file
    writeList(lines, inpfile)
    
def getInpLines(input,token):
    lines=[]
    block=False
    endblock=' '
    key = token[len('phsf.'):]
    if key.startswith('block.'):
        block=True
        key = key[len('block.'):]

    if token in input:
        # If the first element is not a boolean, this contains values
        # to be stored after keyword.
        for element in input[token]: # Takes care of single and multi-line input.
            lines.append(' '.join([str(value) for value in element])) 
          	
        if block:
            lines.insert(0,'%block ' + key.upper())           # next line.
            endblock = '%endblock ' + key.upper()
        else:                                     # Most have arguments on the same line as keyword.
            lines[0] = key.upper() + ' ' + lines[0]
 

    # Add a blank line after each line
    lines.append(endblock)
    lines.append('')

    return lines

def writeList(lines, filename):
    with open(filename, 'w') as f:
        f.write('\n'.join(lines))
        f.write('\n')

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
    print(('Running exectuable: ' + executable))
    # Modified by FDV:
    # Adding the / to make the config more generic
    # Modified by JJK to use os.path.join (even safer than above).
    execList = [os.path.join(execDir,executable)] + args
    inFile = open(os.path.join(workDir,'phsf.in'))
    #print execList
    p = subprocess.Popen(execList, bufsize=0, cwd=workDir, stdin=inFile, stdout=out, stderr=err, encoding='utf8')
    out_r = open(out.name,'r')
    err_r = open(err.name,'r')
    while True:
        try:
            output = out_r.readline()
            error = err_r.readline()
        except:
            pass
        if output == '' and p.poll() is not None:
            break
        if output:
            print("\t" + output.strip())
        if error:
            print("\t" + error.strip())
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
        assert isinstance(input[key], str) 
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
    shift = lambda x: map(lambda a,b: a-b, zip(x, coords[0]))
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
    shift = lambda x: map(lambda a,b: a-b, zip(x,dym['atCoords'][iCenter]))
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

def writeXPSInput(input, phsfinp='input.fdf'):
    
    lines = []
    #setInput(input,'feff.print',[[5,0,0,0,0,0]],Force=True)

    writeInput(input,phsfinp)



