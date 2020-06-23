from structures import Handler, Exchange, Loop, Update
import corvutils.pyparsing as pp
import os, sys, subprocess

dmdwInputFile = 'dmdw.inp' # This is hardcoded in DMDW

# Define dictionary of implemented calculations
implemented = {}
strlistkey = lambda L:','.join(sorted(L))
subs = lambda L:[{L[j] for j in range(len(L)) if 1<<j&k} for k in range(1,1<<len(L))]
for s in subs(['selfen','specfn']):
    key = strlistkey(s)
    autodesc = 'Get ' + ', '.join(s) + ' using DMDW'
    input = ['dynmat','pdos','a2f']
    cost = 5
    implemented[key] = {'type':'Exchange','out':list(s),'req':input,
                        'desc':autodesc,'cost':cost}
implemented['s2'] = {'type':'Exchange','out':['s2'],'req':['dynmat','dmdw.paths'],
                     'desc':'Get XAFS DWF using DMDW','cost':1}
implemented['u2'] = {'type':'Exchange','out':['u2'],'req':['dynmat','dmdw.paths'],
                     'desc':'Get crystallographic DWF using DMDW','cost':1}
implemented['vfe'] = {'type':'Exchange','out':['vfe'],'req':['dynmat','pdos','a2f'],
                          'desc':'Get vibrational free energy using DMDW','cost':2}

class Dmdw(Handler):
    def __str__(self):
        return 'DMDW Handler'

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
            unresolved = {o for o in output if not Dmdw.canProduce(o)}
            canProduce = (o for o in output if Dmdw.canProduce(o))
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
            raise LookupError('Corvus cannot currently produce ' + key + ' using DMDW')
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
            raise LookupError('Corvus cannot currently produce ' + key + ' using DMDW')
        f = lambda subkey : implemented[key][subkey]
        if f('type') is 'Exchange':
            return Exchange(Dmdw, f('req'), f('out'), cost=f('cost'), desc=f('desc'))

    @staticmethod
    def prep(config):
        subdir = config['pathprefix'] + str(config['xcIndex']) + '_DMDW'
        xcDir = os.path.join(config['cwd'], subdir)
        # Make new output directory if if doesn't exist
        if not os.path.exists(xcDir):
            os.mkdir(xcDir)
        # Store current Exchange directory in configuration
        config['xcDir'] = xcDir

    @staticmethod
# Modified by FDV
# Change to move the methods to the internal par of run
#   def run(config, files):
    def run(config, input, output):
        files = generateInput(config, input, output)
        dir = config['xcDir']
        out = open(os.path.join(dir, 'corvus.DMDW.stdout'), 'w')
        err = open(os.path.join(dir, 'corvus.DMDW.stderr'), 'w')
        p = subprocess.Popen([config['dmdw']], cwd=dir, stdout=out, stderr=err)
        p.wait()
        out.close()
        err.close()

# Modified by FDV:
# Bringing translateOutput into run
        translateOutput(config, input, output)

    @staticmethod
    def cleanup(config):
        pass

##### Generic Helper Methods ##########
def check(input, token, default=None):
    if token in input:
        return input[token]
    else:
        return default

def writeList(list, filename):
    with open(filename, 'w') as file:
        for line in list:
            file.write(line + '\n') 

def writeXY(data, filename):
    lines = []
    for i in range(len(data[0])):
        x = data[0][i]
        y = data[1][i]
        lines.append(str(x) + ' '*4 + str(y))
    writeList(lines, filename)

#### Specific Helper Methods
def writeDym(dym, dymFilename, shiftCoords=True):

# Debug FDV:
# Don't shift if dym is of type 1, just for testing
    if dym['dymType'] == 1:
      shiftCoords = False
    # If provided, use specified ordering
    if 'printOrder' in dym:
        atRange = dym['printOrder']
    else:
        atRange = range(dym['nAt'])
    # By default, shift atomic coordinates so central atom is at origin
    if shiftCoords:
        shift = lambda x: map(lambda (a,b): a-b, zip(x,dym['atCoords'][atRange[0]]))
    else:
        shift = lambda x: x
    lines = []
    lines.append('{:5d}'.format(dym['dymType']))
    lines.append('{:5d}'.format(dym['nAt'])) # number of atoms in cluster
    for iAt in atRange:
        lines.append('{:5d}'.format(dym['atNums'][iAt]))
    for iAt in atRange:
        lines.append('{:12.6f}'.format(dym['atMasses'][iAt]))
    for iAt in atRange:
        lines.append(('{:14.8f}'*3).format(*shift(dym['atCoords'][iAt])))
    for iShift, iAt in enumerate(atRange):
        for jShift, jAt in enumerate(atRange):
            lines.append(('{:5d}'*2).format(iShift + 1, jShift + 1))
            lines.append('\n'.join(('{:14.6E}'*3).format(*row) for row in dym['dm'][iAt][jAt]))
    if dym['dymType'] in [2]:
        nUnique = dym['nTypeAt'] # number of unique atom types
        nCell = dym['nCellAt']   # number of atoms in cell 
        lines.append(str(nUnique) + ' '*4 + str(nCell))
        for z, zAtomList in dym['centerAt'].iteritems():
            lines.append(str(z) + ' '*4 + str(len(zAtomList)))
            for zAtom in zAtomList:
                clustIndexStr = str(zAtom['clustIndex'] + 1)
                cellIndexStr = str(zAtom['cellIndex'] + 1)
                coordStr = '  '.join(str(x) for x in zAtom['coord'])
                lines.append('    '.join([clustIndexStr, cellIndexStr, coordStr]))
    writeList(lines, dymFilename)

def readDym(dymFilename):
    integer = pp.Word(pp.nums).setParseAction(lambda t: int(t[0]))
    floating = pp.Word(pp.nums + ".+-E").setParseAction(lambda t: float(t[0]))
    vector = pp.Group(floating * 3)
    matrix = pp.Group(vector * 3)

    dymtype = integer
    ncluster = integer.copy() # Uses unique ParseAction for Forward
    atNums = pp.Forward()     # Placeholders for grabbing ncluster entries
    atMasses = pp.Forward()
    atCoords = pp.Forward()
    def countedParseAction(toks):
        n = toks[0]
        atNums << pp.Group(integer * n)
        atMasses << pp.Group(floating * n)
        atCoords << pp.Group(vector * n)
        return None
    ncluster.addParseAction(countedParseAction)
    dm = pp.Group(pp.OneOrMore(pp.Group(integer*2 + matrix)))
    specifics = pp.SkipTo(pp.stringEnd)

    text = dymtype + ncluster + atNums + atMasses + atCoords + dm + specifics
    with open(dymFilename, 'r') as f:
        try:
            data = text.parseString(f.read()).asList()
        except pp.ParseException as pe:
            print('Parsing Error using pyparsing: invalid input:', pe)
            sys.exit()

    # Construct dym object
    dym = {'dymType':data[0]}
    dym['nAt'] = data[1]
    dym['atNums'] = data[2]
    dym['atMasses'] = data[3]
    dym['atCoords'] = data[4]
    dym['dm'] = [[None for j in xrange(dym['nAt'])] for i in xrange(dym['nAt'])]
    for i,j in ((i,j) for i in xrange(dym['nAt']) for j in xrange(dym['nAt'])):
        index = i * dym['nAt'] + j
        if data[5][index][0] == i+1 and data[5][index][1] == j+1:
            dym['dm'][i][j] = data[5][index][2] 

    # Parse dymtype-specific items
    if dym['dymType'] == 2:
        nUnique = integer
        nCell = integer
        centralZAtom = pp.Group(integer + integer + vector)
        zAtomList = pp.Group(integer + pp.countedArray(centralZAtom))
        typeSpecific = nUnique + nCell + pp.Group(pp.OneOrMore(zAtomList))
    
    # Parse
    try:
        data = typeSpecific.parseString(data[6]).asList()
    except pp.ParseException as pe:
        print('Parsing Error using pyparsing: invalid input:', pe)
        sys.exit()
    
    # Store type-specific information
    if dym['dymType'] == 2:
        dym['nTypeAt'] = data[0]
        dym['nCellAt'] = data[1]
        dym['centerAt'] = {}
        for zAtomList in data[2]:
            z = zAtomList[0]
            dym['centerAt'][z] = []
            for centralZAtom in zAtomList[1]:
                atomInfo = {}
                atomInfo['clustIndex'] = centralZAtom[0]
                atomInfo['cellIndex'] = centralZAtom[1]
                atomInfo['coord'] = centralZAtom[2]
            dym['centerAt'][z].append(atomInfo)

    return dym

def readSpecFn(filename):
    header = pp.OneOrMore(pp.pythonStyleComment)
    pt = pp.Group(pp.Word(pp.nums + ".+-E") * 5)
    text = header.suppress() + pp.OneOrMore(pt)
    f = open(filename, 'r')
    try:
        data = text.parseString(f.read()).asList()
    except pp.ParseException as pe:
        print('Parsing Error using pyparsing: invalid input:', pe)
        sys.exit()
    f.close()
    cols = map(list, zip(*data))
    # w [meV], mag, ph, re, im
    return [cols[0], cols[3], cols[4]]

def combineReImSE(reFilename, imFilename):
    header = pp.OneOrMore(pp.pythonStyleComment)
    pt = pp.Group(pp.Word(pp.nums + ".+-E") * 2)
    text = header.suppress() + pp.OneOrMore(pt)
    reFile = open(reFilename, 'r')
    imFile = open(imFilename, 'r')
    try:
        rePts = text.parseString(reFile.read()).asList()
        imPts = text.parseString(imFile.read()).asList()
    except pp.ParseException as pe:
        print('Parsing Error using pyparsing: invalid input:', pe)
        sys.exit()
    reFile.close()
    imFile.close()
    enCol = map(list, zip(*rePts))[0] # in meV?
    reCol = map(list, zip(*rePts))[1] # in meV?
    imCol = map(list, zip(*imPts))[1]
    se = [enCol,[]]
    for i in range(len(enCol)):
        se[1].append(complex(float(reCol[i]), float(imCol[i])))
    return se 

def read_s2(filename, tempgrid):

    with open(filename, 'r') as f:
        cleanStr = f.read()
    # Define Grammars
    header = pp.OneOrMore(pp.pythonStyleComment).setParseAction(pp.replaceWith(''))
    divider = (pp.Literal('---') + pp.restOfLine).setParseAction(pp.replaceWith(''))
    index = pp.Word(pp.nums)
    path = pp.Literal('Path Indices:').suppress() + pp.Group(pp.OneOrMore(index))
    colon = pp.Literal(':')
    flt = pp.Word(pp.nums + '.+-')
    pathlength = pp.Suppress(pp.Literal('Path Len') + pp.SkipTo(colon) + colon) + flt
# Debug: FDV
#   print tempgrid
#   sys.exit()
#   ntemp = int(tempgrid.split()[0])
    ntemp = int(tempgrid[0])
    if ntemp is 1:
        datum = pp.Group(path + pathlength + flt)
    else:
        colHeader = pp.Suppress(pp.Literal('Temp (K)   s^2 (1e-3 Ang^2)'))
        tempTable = pp.Group(pp.OneOrMore(pp.Group(flt * 2)))
        datum = pp.Group(path + pathlength + colHeader + tempTable)
    data = pp.OneOrMore(datum)
    # Parse String
    try:
        for grammar in [header, divider]:
            cleanStr = grammar.transformString(cleanStr)
        ppresults = data.parseString(cleanStr).asList()
    except pp.ParseException as pe:
        print('Parsing Error using pyparsing: invalid input:', pe)
        sys.exit()
    # Format Results
    results = {}
    for ppr in ppresults:
        key = ','.join(ppr[0])
        if ntemp is 1:
            temps = float(tempgrid.split()[1])
            dwfs = float(ppr[2])
        else:
            temps = map(float, zip(*ppr[2])[0])
            dwfs = map(float, zip(*ppr[2])[1])
        results[key] = {'Pathlength (Ang)':float(ppr[1]), 'Temp (K)':temps, 's^2':dwfs}
    return results

# Added by FDV
# The dmdw output format for the u2 has changed, so we just make a different
# parser setup. There is proably a more pythonic way to set this up, but I don't
# konw enough python at this point to do it better.
def read_u2(filename, tempgrid):

    with open(filename, 'r') as f:
        cleanStr = f.read()
    # Define Grammars
    header = pp.OneOrMore(pp.pythonStyleComment).setParseAction(pp.replaceWith(''))
    divider = (pp.Literal('===') + pp.restOfLine).setParseAction(pp.replaceWith(''))
    index = pp.Word(pp.nums)
    at_index = pp.Literal('Atom Index:').suppress() + pp.Group(pp.OneOrMore(index))
    dir   = pp.Word(pp.alphas,max=1)
    coord_dir = pp.Group(pp.OneOrMore('-')).suppress() + pp.Literal('Direction').suppress() + dir + pp.Group(pp.OneOrMore('-')).suppress()
    colon = pp.Literal(':')
    flt = pp.Word(pp.nums + '.+-')
# Modified by FDV
#   ntemp = int(tempgrid.split()[0])
    ntemp = int(tempgrid[0])
    if ntemp is 1:
        colHeader = pp.Suppress(pp.Literal('u^2 (1e-3 Ang^2):'))
        datum = pp.Group(at_index + coord_dir + colHeader + flt)
    else:
        colHeader = pp.Suppress(pp.Literal('Temp (K)   u^2 (1e-3 Ang^2)'))
        tempTable = pp.Group(pp.OneOrMore(pp.Group(flt * 2)))
        datum = pp.Group(at_index + coord_dir + colHeader + tempTable)
    data = pp.OneOrMore(datum)
    # Parse String
    try:
        for grammar in [header, divider]:
            cleanStr = grammar.transformString(cleanStr)
        ppresults = data.parseString(cleanStr).asList()
    except pp.ParseException as pe:
        print('Parsing Error using pyparsing: invalid input:', pe)
        sys.exit()
    # Format Results
    results = {}
    for ppr in ppresults:
        key = ','.join([ppr[0][0],ppr[1]])
        if ntemp is 1:
            temps = float(tempgrid.split()[1])
            dwfs = float(ppr[2])
        else:
            temps = map(float, zip(*ppr[2])[0])
            dwfs = map(float, zip(*ppr[2])[1])
        results[key] = {'Temp (K)':temps, 'u^2':dwfs}
    return results

def readVFE(filename):
    with open(filename, 'r') as f:
        cleanStr = f.read()
    # Define Grammars
    header = pp.OneOrMore(pp.pythonStyleComment).setParseAction(pp.replaceWith(''))
    divider = (pp.Literal('---') + pp.restOfLine).setParseAction(pp.replaceWith(''))
    flt = pp.Word(pp.nums + '.+-')
    colHeader = pp.Suppress(pp.Literal('Temp (K)      vfe (J/mol-c)'))  
    tempTable = pp.Group(pp.OneOrMore(pp.Group(flt * 2)))
    datum = pp.Group(colHeader + tempTable)
    data = pp.OneOrMore(datum)
    # Parse String
    try:
        for grammar in [header, divider]:
            cleanStr = grammar.transformString(cleanStr)
        ppresults = data.parseString(cleanStr).asList()
    except pp.ParseException as pe:
        print('Parsing Error using pyparsing: invalid input:', pe)
        sys.exit()
    # Format Results
    results = {}
    for ppr in ppresults:
        temps = map(float, zip(*ppr[0])[0])
        vfes = map(float, zip(*ppr[0])[1])
        key = 'Temp (K)'
        if key not in results:
            results[key] = temps
        key = 'VFE (J/mol-c)'
        if key not in results:
            results[key] = vfes
        else:
            results[key] = [sum(x) for x in zip(results[key], vfes)]
    return results

# Modified by FDV:
# Getting the generateInput and translateOutput methods out of the class.
# This makes them visible to all the other methods and usable within the handler

def generateInput(config, input, output):
    dir = config['xcDir']

    # Set type of DMDW calculation based on target output
    if set(output) <= {'s2'}:
        propflag = 0
    elif set(output) <= {'vfe'}:
        propflag = 1
    elif set(output) <= {'selfen','specfn'}:
        propflag = 2
    elif set(output) <= {'u2'}:
        propflag = 3
    else:
        return "Error: unrecognized DMDW target(s): " + output

    # Write dym file
    dymFilename = 'corvus.dym'
    writeDym(input['dynmat'], os.path.join(dir, dymFilename))

    # Single temp --> tempgrid
    if 'dmdw.tempgrid' not in input and 'temp' in input:
        input['dmdw.tempgrid'] = '1 ' + str(input['temp'])

    lines = []
    if propflag in [1,2]: # output = spectral function
        # Replace text options with corresponding DMDW int flag
        dispopt = check(input, 'dispopt', default=0)
        if isinstance(dispopt, int) and dispopt in range(4):
            dispflag = dispopt 
        elif isinstance(dispopt, basestring):
            if dispopt.isdigit():
                dispflag = int(dispopt)
            else:
                flagmap = {'x':1, 'y':2, 'z':3, 'full':0}
                if dispopt in flagmap:
                    dispflag = flagmap[dispopt]
                else:
                    return "DMDW Input Error: unrecognized dispopt: " + dispopt
        else:
            return "DMDW Input Error: unrecognized dispopt: " + dispopt
            

        # Write PDOS and electron-phonon couplings
        pdosFilename = 'corvus.pdos'
        a2fFilename = 'corvus.a2f'
        #a2Filename = 'corvus.a2'
        writeXY(input['pdos'], os.path.join(dir, pdosFilename))
        writeXY(input['a2f' ], os.path.join(dir, a2fFilename))
        #writeXY(input[ 'a2' ], os.path.join(dir, a2Filename))
        lines.append(check(input, 'dmdw.ioflag', default='0'))
        lines.append(check(input, 'dmdw.nlanc', default='16'))
        lines.append(check(input, 'dmdw.tempgrid', default='3 100.0 300.0')) 
        lines.append(str(propflag))
        lines.append(str(dispflag))
        ekopt = check(input, 'ekopt', default=1) # meV default
        ekval = check(input, 'ek', default=1.0)
        lines.append(str(ekopt) + ' ' + str(ekval))
        lines.append(dymFilename)
        lines.append(pdosFilename)
        lines.append(a2fFilename)
        #lines.append(a2Filename)

    elif propflag in [0,3]:
        lines.append(check(input, 'dmdw.ioflag', default='0'))
        lines.append(check(input, 'dmdw.nlanc', default='16'))
        lines.append(check(input, 'dmdw.tempgrid', default='3 100.0 300.0')) 
        lines.append(str(propflag))
        lines.append(dymFilename)
        if propflag is 0:
            defaultpaths = '1\n2 1 0 30.0'
        if propflag is 3:
            defaultpaths = '1\n1 0 30.0'
        lines.append(check(input, 'dmdw.paths', default=defaultpaths))

    else:
        return "DMDW Input Error: unrecognized property flag: " + propflag

    print('lines')
    print(lines)
    writeList(lines, os.path.join(dir, dmdwInputFile))
    return [dmdwInputFile]

def translateOutput(config, input, output):
    dir = config['xcDir']
    for target in output:
        if target == 'specfn':
            output[target] = readSpecFn(os.path.join(dir, 'dmdw_Akw.dat'))
        elif target == 'selfen':
            reFilename = os.path.join(dir, 'dmdw_reSE_a2F.dat')
            imFilename = os.path.join(dir, 'dmdw_imSE_a2F.dat')
            output[target] = combineReImSE(reFilename, imFilename)
        elif target == 's2':
            filename = os.path.join(dir, 'dmdw.out')
            output[target] = read_s2(filename, input['dmdw.tempgrid'][0])
        elif target == 'u2':
            filename = os.path.join(dir, 'dmdw.out')
            output[target] = read_u2(filename, input['dmdw.tempgrid'][0])
        elif target == 'vfe':
            filename = os.path.join(dir, 'dmdw.out')
            output[target] = readVFE(filename)

