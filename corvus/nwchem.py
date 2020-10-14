from corvus.structures import Handler, Exchange, Loop, Update
import corvutils.pyparsing as pp
import os, sys, subprocess, math

# Debug: FDV
import pprint
pp_debug = pprint.PrettyPrinter(indent=4)

# Naming grammar for *.files: [xcNum].[tool].[description].files
filesGrammar = pp.delimitedList(pp.Word(pp.alphanums), delim='.')

# Define dictionary of implemented calculations
implemented = {}

# The following are the basic requirement for even the most minimal VASP
# calculation.
basics = ['cell_struc_xyz_red']
strlistkey = lambda L:','.join(sorted(L))
subs = lambda L:[{L[j] for j in range(len(L)) if 1<<j&k} for k in range(1,1<<len(L))]
for s in subs(['ene_int','aopt','qmd','xas','dynmat','opt_dynmat']):
    key = strlistkey(s)
    autodesc = 'Get ' + ', '.join(s) + ' using Nwchem'
    input = basics
    cost =  2
    implemented[key] = {'type':'Exchange','out':list(s),'req':input,
                        'desc':autodesc,'cost':cost}

implemented['xas_avg'] = {'type':'Exchange','out':['xas_avg'],'req':basics+['qmd'],
                          'desc':'Get xas_avg using nwchem','cost':10}
implemented['opt_qmd'] = {'type':'Exchange','out':['opt_qmd'],'req':basics+['aopt'],
                          'desc':'Get opt_qmd using nwchem','cost':2}
#implemented['opt_qmd'] = {'type':'Exchange','out':['opt_qmd'],'req':basics+['aopt']+['qmd'],
                          #'desc':'Get opt_qmd using nwchem','cost':2}
implemented['opt_xas'] = {'type':'Exchange','out':['opt_xas'],'req':basics+['aopt'],
                          'desc':'Get opt_xas using nwchem','cost':2}
implemented['opt_xasavg'] = {'type':'Exchange','out':['opt_xasavg'],'req':basics+['opt_qmd'],
                             'desc':'Get optimized xas_avg using nwchem','cost':10}
implemented['user_xasavg'] = {'type':'Exchange','out':['user_xasavg'],'req':basics,
                          'desc':'Get xas_avg from user qmd snapshots using nwchem','cost':10}

# Debug: FDV
#pp_debug.pprint(implemented)
#sys.exit()

class Nwchem(Handler):
    def __str__(self):
        return 'Nwchem Handler'

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
            unresolved = {o for o in output if not Nwchem.canProduce(o)}
            canProduce = (o for o in output if Nwchem.canProduce(o))
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
    def costOf(output):
        if isinstance(output, list) and output and isinstance(output[0], str):
            key = strlistkey(output)
        elif isinstance(output, str):
            key = output
        else:
            raise TypeError('Output should be token or list of tokens')
        if key not in implemented:
            raise LookupError('Corvus cannot currently produce ' + key + ' using Nwchem')
        return implemented[key]['cost']

    @staticmethod
    def sequenceFor(output):
        if isinstance(output, list) and output and isinstance(output[0], str):
            key = strlistkey(output)
        elif isinstance(output, str):
            key = output
        else:
            raise TypeError('Output should be token or list of tokens')
        if key not in implemented:
            raise LookupError('Corvus cannot currently produce ' + key + ' using Nwchem')
        f = lambda subkey:implemented[key][subkey]
        if f('type') is 'Exchange':
            return Exchange(Nwchem, f('req'), f('out'), cost=f('cost'), desc=f('desc'))
        
    @staticmethod
    def prep(config):
        subdir = config['pathprefix'] + str(config['xcIndex']) + '_Nwchem'
        xcDir = os.path.join(config['cwd'], subdir)
        # Make new output directory if it doesn't exist
        if not os.path.exists(xcDir):
            os.mkdir(xcDir)
        # Remove any old output files
        else :
            for file in os.listdir(xcDir):
                if file.endswith('.out'):
                    os.remove(os.path.join(xcDir, file))
        # Store current Exchange directory in configuration
        config['xcDir'] = xcDir
    
    @staticmethod
# Modified by FDV:
# Change to move the methods to the internal par of run
#   def run(config, files):
    def run(config, input, output):

# Modified by FDV:
# Bringing generateInput into run
        files = generateInput(config, input, output)
#       print(files)
#       sys.exit()
        # Run by looping through input files
        dir = config['xcDir']
        for f in files:
            print(f)
            print('hellofile')
            (num, tool, desc, suffix) = filesGrammar.parseString(f)
            if config['parallelRun']:
                if tool == 'nwchem':
                    prefix = config['parallelRun']
                else:
                    prefix = config['parallelRun'].split()[0] + ' -n 1' 
            else:
                prefix = ''
            executable = (prefix + ' ' + config[tool] + ' ' + f).split()
                
            inp = open(os.path.join(dir, f), 'r')
            log = open(os.path.join(dir, '.'.join([num,tool,desc,'out'])), 'w')
            print((executable,dir,log))
            p = subprocess.run(executable, cwd=dir, stdin=None, stdout=log, stderr=log)
            print(files)
            inp.close()
            log.close()

# Modified by FDV:
# Bringing translateOutput into run
        translateOutput(config, input, output)
	
    @staticmethod
    def cleanup(config):    
        dir = config['xcDir']
        for file in os.listdir(dir):
            if 'CLEANUP' in file:
                os.remove(os.path.join(dir, file))

##### Generic Helper Methods ##########
def check(input, token, default=None):
    if token in input:
        return input[token]
    else:
        return default

def writeDict(inDict, filename):
    with open(filename, 'w') as file:
        for token, value in iter(sorted(inDict.items())):
            file.write(str(value) + '\n')

def writeList(list, filename):
    with open(filename, 'w') as file:
        for line in list:
            file.write(line + '\n') 

def expandedList(string, length=3):
    assert isinstance(string, str)
    value = pp.Word(pp.nums + ".+-E")
    prefix = pp.Combine(pp.Optional(pp.Word(pp.nums)) + pp.Literal('*'))
    numList = pp.OneOrMore(pp.Optional(prefix) + value)
    expanded = []
    num = 1
    for x in numList.parseString(string).asList():
        if num is not 1:
            expanded.extend([x]*num)
            num = 1
        elif '*' in x:
            if x[:-1] is '':
                num = length
            else:
                num = int(x[:-1])
        else:
            expanded.append(x)
    return expanded
    
#### Specific Input Methods #########
def nwchemFiles(nwPrefix, pspList):
    files = []
    files.append(nwPrefix + '.in')
    files.append(nwPrefix + '.out')
    files.append(nwPrefix + '.i')
    files.append(nwPrefix + '.o')
    files.append(nwPrefix + '.CLEANUP')
    # Get absolute file location(s)
    for psp in pspList:
        files.append(os.path.abspath(psp))
    return files

def ddbFiles(mddbPrefix, datasetPrefix, qptRange):
    files = []
    files.append(mddbPrefix + '.MDDB')
    files.append('Merged Derivative Database')
    files.append(str(len(qptRange)))
    for i in qptRange:
        files.append(datasetPrefix + str(i) + '_DDB')
    return files

def gkkFiles(mgkkPrefix, datasetPrefix, qptRange, atdirRange):
    files = []
    files.append(mgkkPrefix + '.MGKK_bin')
    files.append('0  # binary output')
    files.append(datasetPrefix + '2_WFK')
    numfiles = len(qptRange)*len(atdirRange)
    files.append(' '.join(['0', str(numfiles), str(numfiles)]))
    for i in qptRange:
        for j in atdirRange:
            files.append(datasetPrefix + str(i) + '_GKK' + str(j))
    return files

def anaddbFiles(anaPrefix, mddbFile, mgkkFile):
    files = []
    files.append(anaPrefix + '.in')
    files.append(anaPrefix + '.out')
    files.append(mddbFile)
    files.append(anaPrefix + '.band2eps')
    files.append(mgkkFile)
    files.append(anaPrefix + '.dummyLine1')
    files.append(anaPrefix + '.dummyLine2')
    return files

def baseInput(input):

    import corvutils.corvus_misc_funcs as cu

    dict = {}
    temp = {}
    s = ' '
    seq0 = ('geometry','units','angstroms','print','xyz','noautosym')
    seq0_ = str(s.join(seq0))
    xc = check(input,'nwchem.xc', default = 'b3lyp')
    xcstr = ''.join(str(r) for v in xc for r in v)
    seq1 = ('xc', xcstr)
    seq1_ = str(s.join(seq1))
    multstr = spinmult(input)
    #multstr = ''.join(str(r) for v in mult for r in v)
    seq2 = ('mult', multstr)
    seq2_ = str(s.join(seq2))
# NOTE FDV:
# Hack to get better DMs
#   seq3 = ('grid','nodisk')
    seq3 = ('grid','nodisk','xfine')
    seq3_ = str(s.join(seq3))
    seq4 = ('convergence','nolevelshifting')
    seq4_ = str(s.join(seq4))
    temp['method1'] = 'dft'
    temp['dft_xc'] = seq1_
    temp['dft_mult'] = seq2_
    temp['dft_disc'] = seq3_
    temp['end'] = 'end'
    ite = check(input,'nwchem.iter', default = '200')
    iterstr = ''.join(str(r) for v in ite for r in v)
    iterseq = ('iterations', iterstr)
    iteration = str(s.join(iterseq))
    prnt = 'print "final vectors analysis"'
    coordstring = [[str(j) for j in i] for i in input['cell_struc_xyz_red']]
    b = [s.join(x) for x in coordstring]
    c = '\n'.join(b)
    targ = check(input,'target_list', default = None)
    target = ''.join(str(r) for v in targ for r in v)
    if target == 'qmd':
        basisH = 'H library "def2-svp"'
        basiselse = '* library "def2-tzvp"'
        basisset = basisH + '\n' + basiselse
    elif target == 'aopt':
        basisH = 'H library "def2-svp"'
        basiselse = '* library "def2-tzvp"'
        basisset = basisH + '\n' + basiselse
    elif target == 'dynmat':
        basisH = 'H library "def2-svp"'
        basiselse = '* library "def2-tzvp"'
        basisset = basisH + '\n' + basiselse
# NOTE FDV: This seems a bit wrong. s2 is not really a target for nwchem
# so setting variables based on "target" is not exactly the right way to
# proceed. Hacking it for now.
    elif target == 's2':
        basisH = 'H library "def2-svp"'
        basiselse = '* library "def2-tzvp"'
        basisset = basisH + '\n' + basiselse
    elif target == 'ene_int':
        basisH = 'H library "def2-svp"'
        basiselse = '* library "def2-tzvp"'
        basisset = basisH + '\n' + basiselse
    elif target == 'opt_qmd':
        basisH = 'H library "def2-svp"'
        basiselse = '* library "def2-tzvp"'
        basisset = basisH + '\n' + basiselse
    elif target == 'xas':
        basisH = 'H library "def2-svp"'
        basiselse = '* library "def2-svpd"'
        X = check(input, 'nwchem.xaselem', default = None)
        if X is not None:
            Xstr = ''.join(str(r) for v in X for r in v)
            basisX = Xstr + ' library "def2-tzvpd"'
            basisset = basisH + '\n' + basiselse + '\n' + basisX
        elif X is None:
            basisset = basisH + '\n' + basiselse
            print('Must specify target element, if no energy window, for x-ray spectroscopy')
    elif target == 'xas_avg':
        basisH = 'H library "def2-svp"'
        basiselse = '* library "def2-svpd"'
        X = check(input, 'nwchem.xaselem', default = None)
        if X is not None:
            Xstr = ''.join(str(r) for v in X for r in v)
            basisX = Xstr + ' library "def2-tzvpd"'
            basisset = basisH + '\n' + basiselse + '\n' + basisX
        elif X is None:
            basisset = basisH + '\n' + basiselse
            print('Must specify target element, if no energy window, for x-ray spectroscopy')
    elif target == 'opt_xas':
        basisH = 'H library "def2-svp"'
        basiselse = '* library "def2-svpd"'
        X = check(input, 'nwchem.xaselem', default = None)
        if X is not None:
            Xstr = ''.join(str(r) for v in X for r in v)
            basisX = Xstr + ' library "def2-tzvpd"'
            basisset = basisH + '\n' + basiselse + '\n' + basisX
        elif X is None:
            basisset = basisH + '\n' + basiselse
            print('Must specify target element, if no energy window, for x-ray spectroscopy')
    elif target == 'opt_xasavg':
        basisH = 'H library "def2-svp"'
        basiselse = '* library "def2-svpd"'
        X = check(input, 'nwchem.xaselem', default = None)
        if X is not None:
            Xstr = ''.join(str(r) for v in X for r in v)
            basisX = Xstr + ' library "def2-tzvpd"'
            basisset = basisH + '\n' + basiselse + '\n' + basisX
        elif X is None:
            basisset = basisH + '\n' + basiselse
            print('Must specify target element, if no energy window, for x-ray spectroscopy')
    elif target == 'user_xasavg':
        basisH = 'H library "def2-svp"'
        basiselse = '* library "def2-svpd"'
        X = check(input, 'nwchem.xaselem', default = None)
        if X is not None:
            Xstr = ''.join(str(r) for v in X for r in v)
            basisX = Xstr + ' library "def2-tzvpd"'
            basisset = basisH + '\n' + basiselse + '\n' + basisX
        elif X is None:
            basisset = basisH + '\n' + basiselse
            print('Must specify target element, if no energy window, for x-ray spectroscopy') 

# NOTE FDV: Here we have a small problem now, because the format of the data
# that is in the input (system dict) is different than the one the default
# (basisset) has. We can try to get the input into the default format or the
# other way around. I have decided for now to modify the input to match the
# requirements here. This is simpler than the opposite, but is not as pleasing
# from a coding point of view.

#   basis = check(input,'nwchem.basis', default= basisset)
    if 'nwchem.basis' in list(input.keys()):
# Here we convert the input list into a single string with the right format
      basis = cu.StrList2Str(input['nwchem.basis'])
    else:
      basis = basisset

# Debug: FDV
#   print(basis)
#   sys.exit()
    basisstr = ''.join(str(r) for v in basis for r in v)
#   basisstr = ' '.join(str(r) for v in basis for r in v)
    chrg = check(input,'nwchem.charge', default = '0')
    chrgstr = ''.join(str(r) for v in chrg for r in v)
    chrgseq = ('charge', chrgstr)
    charge = str(s.join(chrgseq))
    temp['echo'] = 'echo'
    temp['charge'] = charge
    temp['mem'] = 'memory 4096 mb'
    dict['coords'] = seq0_ + '\n' + c + '\n' + temp['end'] + '\n'
    dict['basis'] = temp['echo'] + '\n' + '\n' + temp['charge'] + '\n'+ '\n'+ temp['mem'] + '\n' + '\n' + 'basis' + '\n' + basisstr + '\n' + temp['end'] + '\n'
    dict['method'] = temp['method1'] + '\n' + temp['dft_xc'] + '\n' + temp['dft_mult'] + '\n' + temp['dft_disc'] + '\n' + seq4_ + '\n' + iteration + '\n' + 'direct' + '\n' + prnt + '\n' + temp['end'] + '\n'

# NOTE FDV:
# Hack to get better DMs
    dict['addendum'] = 'scf'+'\n'+'thresh 1e-10'+'\n'+'end'+'\n'+'driver'+'\n'+'tight'+'\n'+'maxiter 80'+'\n'+'end'

    print(dict)
    return dict

def optInput(input):
    dict = baseInput(input)
    sp = ' '
    seq3 = ('task','dft','optimize')
    seq3_ = str(sp.join(seq3))
    dict['task'] = seq3_
    return dict

def dynmatInput(input):
    dict = baseInput(input)
    sp = ' '
    seq3 = ('task','dft')
    seq3_ = str(sp.join(seq3))
    dict['task'] = 'task dft freq\n'
    return dict

def opt_dynmatInput(input):
    dict = baseInput(input)
    sp = ' '
    seq3 = ('task','dft','optimize')
    seq3_ = str(sp.join(seq3))
    dict['task'] = 'task dft optimize\ntask dft freq\n'
    return dict

def qgridInput(input):
    dict = baseInput(input)
    dict['kptopt'] = 1
    dict['ngkpt'] = input['ngqpt']
    dict['nshiftk'] = 1
    dict['shiftk'] = '0 0 0'
    dict['prtvol'] = -1
    dict['ecut'] = 10.0
    return dict
    
def qmdInput(input):
    dict = baseInput(input)
    a = ' '
    emp = {}
    emp['method'] = 'qmd'
    nstep = check(input, 'nwchem.qmd.nstep_nucl', default = '2000')
    nstepstr = ''.join(str(r) for v in nstep for r in v)
    seq_step = ('nstep_nucl', nstepstr)
    seq_stepj = str(a.join(seq_step))
    dt = check(input, 'nwchem.qmd.dt_nucl', default = '10')
    dtstr = ''.join(str(r) for v in dt for r in v)
    seq_dt = ('dt_nucl', dtstr)
    seq_dtj = str(a.join(seq_dt))
    temp = check(input, 'nwchem.qmd.targ_temp', default = '200')
    tempstr = ''.join(str(r) for v in temp for r in v)
    seq_temp = ('targ_temp', tempstr)
    seq_tempj = str(a.join(seq_temp))
    therm = check(input, 'nwchem.qmd.thermostat', default = 'svr 100')
    thermstr = ''.join(str(r) for v in therm for r in v)
    seq_therm = ('thermostat', thermstr)
    seq_thermj = str(a.join(seq_therm))
    xyz = check(input, 'nwchem.qmd.print_xyz', default = '5')
    xyzstr = ''.join(str(r) for v in xyz for r in v)
    seq_xyz = ('print_xyz',xyzstr)
    seq_xyzj = str(a.join(seq_xyz))
    emp['end'] = 'end'
    dict['method2'] = emp['method']+'\n'+seq_stepj+'\n'+seq_dtj+'\n'+seq_tempj+'\n'+seq_thermj+'\n'+seq_xyzj+'\n'+emp['end']+'\n'
    seq4 = ('task','dft','qmd')
    seq4_ = str(a.join(seq4))
    dict['task2'] = seq4_
    return dict

def UInput(input):
    dict = baseInput(input)
    #dict['toldfe'] = check(input, 'toldfe', default='1.0d-12')
    #dict['nstep'] = check(input, 'nstep', default=50)
    sp = ' '
    seq3 = ('task','dft','energy')
    seq3_ = str(sp.join(seq3))
    dict['task'] = seq3_
    return dict
    
def xasInput(input):
    dict = baseInput(input)
    b = ' '
    xas = {}
    xas['method'] = 'tddft'
    alpha = check(input, 'nwchem.xas.alpha', default = None)
    ewin = check(input, 'nwchem.xas.xrayenergywin', default = None)
    beta = ''
    if alpha is not None:
        alphastr = ''.join(str(r) for v in alpha for r in v)
        seq_alpha = ('alpha', alphastr)
        seq_alphaj = str(b.join(seq_alpha))
        alphaflag = '\n' + seq_alphaj
        seq_beta = ('beta', alphastr)
        beta = str(b.join(seq_beta))
    elif alpha is None:
        alphaflag = ''
#       seq_beta = ('beta', '1 1')
#       beta = str(b.join(seq_beta))
        print('Need either alpha flag, energy window, or spectroscopy element for x-ray spectroscopy')
    if ewin is not None:
        ewinstr = '  '.join(str(r) for v in ewin for r in v)
        seq_ewin = ('ewin', ewinstr)
        seq_ewinj = str(b.join(seq_ewin))
        ewinflag = '\n' + seq_ewinj
#       seq_beta = ('beta', '1 1')
#        beta = str(b.join(seq_beta))
    elif ewin is None:
        ewinflag = ''
#       seq_beta = ('beta', '1 1')
#        beta = str(b.join(seq_beta))
        print('Need either alpha flag, energy window, or spectroscopy element for x-ray spectroscopy')
    xas['cis'] = 'cis'
    nroots = check(input, 'nwchem.xas.nroots', default = '10')
    nrootsstr = ''.join(str(r) for v in nroots for r in v)
    seq_nroots = ('nroots', nrootsstr)
    seq_nrootsj = str(b.join(seq_nroots))
    xas['trip'] = 'notriplet'
    xas['end'] = 'end'
    if int(nrootsstr)*15 > 100:
        vec = str(int(nrootsstr)*15)
    elif int(nrootsstr)*15 < 100:
        vec = '100'
    maxvec = check(input,'nwchem.xas.vec', default = vec)
    maxvecstr = ''.join(str(r) for v in maxvec for r in v)
    seq_vec = ('maxvecs', maxvecstr)
    maxvecs = str(b.join(seq_vec))
    if spinmult(input) == '2':
        xas['method1'] = xas['method'] + alphaflag + ewinflag + '\n' + beta + '\n'
    elif spinmult(input) == '1':
        xas['method1'] = xas['method'] + alphaflag + ewinflag + '\n'
    dict['method2'] = xas['method1'] + xas['cis'] + '\n' + seq_nrootsj + '\n' + xas['trip'] + '\n' + maxvecs + '\n' + xas['end'] + '\n'
    seqtask = ('task', 'tddft')
    dict['task'] = str(b.join(seqtask))
    print(dict)
    return dict

def elphonInput(input):
    dict = baseInput(input)
    qpts = input['qptlist']
    dict['ndtset'] = 2 + len(qpts)
    # Density calculation
    dict['tolvrs1'] = check(input, 'tolvrs1', default='1.0d-18')
    dict[ 'nline1'] = check(input, 'nline1', default=8)
    dict['rfphon1'] = 0
    dict[  'nqpt1'] = 0
    dict['getwfk1'] = 0
    dict['prtden1'] = 1
    dict['kptopt1'] = 1
    # Wavefunction calculation
    dict['tolwfr2'] = check(input, 'tolwfr2', default='1.0d-22')
    dict[  'iscf2'] = check(input, 'iscf2', default=-3)
    dict['rfphon2'] = 0
    dict[  'nqpt2'] = 0
    dict['getwfk2'] = 0
    dict['getden2'] = 1
    dict[ 'prtwf2'] = 1
    # Calculation for each qpt
    for i, qpt in enumerate(qpts):
        dict['qpt' + str(i+3)] = qpt
    dict[    'ixc'] = check(input, 'ixc', default=7)
    dict[ 'diemac'] = check(input, 'diemac', default='1.0E+06')
    dict[ 'tolvrs'] = check(input, 'tolvrs', default='1.0e-10')
    dict[ 'istwfk'] = '*1'
    dict[ 'rfphon'] = 1
    dict[   'nqpt'] = 1
    dict[  'rfdir'] = '1 1 1'
    dict['rfatpol'] = '1 ' + str(input['natom'])
    dict[  'prtwf'] = 0
    dict[ 'getwfk'] = 2
    dict['prepgkk'] = 1
    dict[ 'prtgkk'] = 1
    dict[ 'kptopt'] = 3
    dict[  'ngkpt'] = check(input, 'ngkpt', default='4 4 4')
    dict['nshiftk'] = 1
    dict[ 'shiftk'] = '0.0 0.0 0.0'
    dict[   'ecut'] = check(input, 'ecut', default='50.0')
    if 'ecutsm' in input:
        dict['ecutsm'] = input['ecutsm']
    dict[  'nstep'] = check(input, 'nstep', default=800)
    if 'verbatim' in input:
        dict[''] = input['verbatim']
    return dict

def ifcInput(input):
    dict = {}
    dict[ 'elphflag'] = 1
    dict[   'prtdos'] = 1
    dict[  'prt_ifc'] = 1
    dict[  'ifcflag'] = 1
    dict[   'ifcana'] = 1
    dict[   'ifcout'] = check(input, 'ifcout', default=200) #TODO: create shell/nifc function
    dict[   'natifc'] = input['natom']
    dict[    'atifc'] = ' '.join(str(i) for i in range(1, 1 + int(input['natom'])))
    dict[    'ngqpt'] = check(input, 'ngqpt', default='4 4 4')
    dict[   'ng2qpt'] = check(input, 'ng2qpt', default='16 16 16')
    dict[      'asr'] = check(input, 'asr', default=2)
    dict['symdynmat'] = check(input, 'symdynmat', default=0)
    dict[    'nph1l'] = check(input, 'nph1l', default=1)
    dict[    'qph1l'] = check(input, 'qph1l', default='0.0 0.0 0.0 1')
    dict[   'nqpath'] = check(input, 'nqpath', default=5)
    qpaths = '0.0 0.0 0.0\n0.0 0.0 1/8\n0.0 0.0 1/4\n 0.0 0.0 1/2\n0.0 0.0 1\n'
    dict[    'qpath'] = check(input, 'qpath', default=qpaths)
    dict[   'mustar'] = check(input, 'mustar', default=0.1)
    return dict

#### Specific Translation Methods #########
def getacell(file):
    converge = pp.Literal('Optimization converged')
    geometry = pp.Literal('Output coordinates in angstroms')
    dash = pp.Word('-')+pp.Word('-')+pp.Word('-')+pp.Word('-')+pp.Word('-')+pp.Word('-')
    end = pp.Literal('Atomic Mass')
    vector = pp.Group(pp.Word(pp.nums + ".+-E")*3)
    #align = pp.Group(pp.Word(pp.nums) + pp.Word(pp.alphas) + pp.OneOrMore(pp.Word(pp.nums + ".+-E")))
    #end = pp.Literal('END DATASET(S)')
    acell = (pp.SkipTo(converge) + pp.SkipTo(geometry) + pp.SkipTo(dash, include = True)).suppress()
    bcell = acell + pp.SkipTo(end, include = False)
    f = open(file, 'r')
    try:
        a = bcell.parseString(f.read()).asList()[0]
    except pp.ParseException as pe:
        print(('Parsing Error using pyparsing: invalid input:', pe))
        sys.exit()
    f.close()
    #return ''.join(a)
    b = a.split("\n")
    e = []
    for line in b:
        c = line.split()
        if len(c)>5:
            d = [c[1],c[3],c[4],c[5]]
            e.append(d)
    e = list(map(list, list(zip(*e))))     
    return e
    
def getdynmat(input,file):

  from . import conversions as cn
  import corvutils.corvus_misc_funcs as cu

  Buffer = []
  for lineraw in open(file,'r'):
    Buffer.append(float(lineraw.rstrip().lstrip().replace('D','e')))

  nFC = len(Buffer)
  nCoord = int((math.sqrt(8*nFC+1)-1)/2);
  nAtoms = nCoord/3;

# Create the empty FC matrix
  FC = [[None for i in range(nCoord)] for j in range(nCoord)]

# Transform the cartesian force constants vector buffer into a symmetric array
  iE = 0
  for i in range(nCoord):
    for j in range(i+1):
      FC[i][j] = Buffer[iE]
      FC[j][i] = FC[i][j]
      iE += 1

# Convert the symmetric array of FCs into blocked format
  FCb = [[ [[None for p in range(3)] for q in range(3)] \
           for iAt in range(nAtoms)] for jAt in range(nAtoms)]
  for iAtom in range(nAtoms):
    for jAtom in range(nAtoms):
      for p in range(3):
        for q in range(3):
          FCb[iAtom][jAtom][p][q] = FC[3*iAtom+p][3*jAtom+q]

# Prepare the dym output
  dym = {'dymType':1}
  dym['nAt'] = nAtoms

# Debug
  pp_debug.pprint(input['cell_struc_xyz_red'])
# sys.exit()

  dym['atNums']   = [ cu.AtSym2AtNum(atominf[0]) for atominf in input['cell_struc_xyz_red'] ]

# NOTE FDV:
# Setting these to -1 means that dmdw will use its own set of masses
  dym['atMasses'] = [ -1.0 ]*nAtoms

  dym['atCoords'] = [ cn.angstrom2bohr(atominf[1:]) for atominf in input['cell_struc_xyz_red'] ]

  dym['dm'] = FCb

  return dym
    
def optqmdcell(file):
    converge = pp.Literal('Optimization converged')
    geometry = pp.Literal('Output coordinates in angstroms')
    dash = pp.Word('-')+pp.Word('-')+pp.Word('-')+pp.Word('-')+pp.Word('-')+pp.Word('-')
    end = pp.Literal('Atomic Mass')
    vector = pp.Group(pp.Word(pp.nums + ".+-E")*3)
    #align = pp.Group(pp.Word(pp.nums) + pp.Word(pp.alphas) + pp.OneOrMore(pp.Word(pp.nums + ".+-E")))
    #end = pp.Literal('END DATASET(S)')
    acell = (pp.SkipTo(converge) + pp.SkipTo(geometry) + pp.SkipTo(dash, include = True)).suppress()
    bcell = acell + pp.SkipTo(end, include = False)
    f = open(file, 'r')
    try:
        a = bcell.parseString(f.read()).asList()[0]
    except pp.ParseException as pe:
        print(('Parsing Error using pyparsing: invalid input:', pe))
        sys.exit()
    f.close()
    #return ''.join(a)
    b = a.split("\n")
    e = []
    for line in b:
        c = line.split()
        if len(c)>5:
            d = [c[1],c[3],c[4],c[5]]
            e.append(d)	     
# Added by FDV:
# Before we return, we convert the coordinates to floats, as
# they should be
    e = [ [ l[0], float(l[1]), float(l[2]), float(l[3]) ] for l in e ]
#   print e
#   sys.exit()
    return e

def getagrid(acell, n):
    if isinstance(acell, str): 
        a0 = list(map(float, expandedList(acell)))
    start = 0.98
    end = 1.03
    step = (end - start) / (n - 1.0)
    grid = []
    for i in range(n):
        scale = start + i * step
        grid.append('  '.join([str(x * scale) for x in a0]))
    return grid

def getU(file, addDS=''):
    token = pp.Literal('Total DFT energy =' + addDS)
    value = pp.Word(pp.nums + ".+-E")
    #end = pp.Literal('END DATASET(S)')
    etotal = (pp.SkipTo(token) + token).suppress() + value
    f = open(file, 'r')
    try:
        u = etotal.parseString(f.read()).asList()[0]
    except pp.ParseException as pe:
        print(('Parsing Error using pyparsing: invalid input:', pe))
        sys.exit()
    f.close()
    return float(u)
    
def snapshots(file,input):
    nstep = check(input, 'nwchem.qmd.nstep_nucl', default = '2000')
    nstepstr = ''.join(str(r) for v in nstep for r in v)
    nstep1 = int(nstepstr)
    half = int(nstep1/10)
    #print nstep1
    numatoms = int(open(file).readline().rstrip(' '))
    n1 = int(numatoms + 2)
    snapstep = check(input, 'nwchem.qmd.snapstep', default = '5')
    snapstepstr = ''.join(str(r) for v in snapstep for r in v)
    stepint = int(snapstepstr)
    with open(file, 'r') as f:
        data = f.readlines()
        data1 = [x.strip('\n') for x in data]
        result = [data1[i:i+n1] for i in range(0,len(data),n1)]
        #for x in result:
            #del x[0:2]
        #for x in result:
            #x = map(list, zip(*x))
        result = result[-half:]
        res = result[::stepint]
        s = '\n'.join(str(r) for v in res for r in v)
    return s   
        
def getSpect(file):
    eV = []
    osc = []
    with open(file) as f:
        for line in f:
            if 'Root' in line:
                parts = line.split()
                eV.append(parts[-2])
            elif 'Dipole Oscillator Strength' in line:
                osc.append(line)
    dipole = [' '.join(x.split()[3:4]) for x in osc]
    return eV, dipole

def avgSpect(input):
    xasgrid = input['xas-grid']
    xas = [x[1] for x in xasgrid]
    x = [[float(j) for j in i] for i in xas]
    averages = [sum(items)/len(x) for items in zip(*x)]
    avgstr = [str(item) for item in averages]
    ev = [x[0] for x in xasgrid]
    energy = ev[0]
    return energy, avgstr

def kgrid2qgrid(file):
    token = pp.Keyword('kpt')
    qpt = pp.Group(pp.Word(pp.nums + ".+-E")*3)
    kpt = (pp.SkipTo(token) + token).suppress() + pp.OneOrMore(qpt)
    f = open(file, 'r')
    try:
        pts = kpt.parseString(f.read()).asList()
    except pp.ParseException as pe:
        print(("Parsing Error using pyparsing: invalid input:", pe))
        sys.exit()
    f.close()
    qpts = []
    for p in pts:
        qpts.append(' '.join(p))
    return qpts

def anaddb2cols(file):
    header = pp.OneOrMore(pp.pythonStyleComment)
    pt = pp.Group(pp.Word(pp.nums + ".+-E") * 2)
    text = header.suppress() + pp.OneOrMore(pt)
    f = open(file, 'r')
    try:
        data = text.parseString(f.read()).asList()
    except pp.ParseException as pe:
        print(("Parsing Error using pyparsing: invalid input:", pe))
        sys.exit()
    f.close() 
    return list(map(list, list(zip(*data))))

def eli2couplings(filePrefix):
    header = pp.OneOrMore(pp.pythonStyleComment)
    pt = pp.Group(pp.Word(pp.nums + ".+-E") * 2)
    text = header.suppress() + pp.OneOrMore(pt)
    pdsFile = open(filePrefix + 'PDS', 'r')
    a2fFile = open(filePrefix + 'A2F', 'r')
    try:
        pdsRow = text.parseString(pdsFile.read()).asList()
        a2fRow = text.parseString(a2fFile.read()).asList()
    except pp.ParseException as pe:
        print(('Parsing Error using pyparsing: invalid input:', pe))
        sys.exit()
    pdsFile.close()
    a2fFile.close()
    enCol = list(map(list, list(zip(*pdsRow))))[0]
    pdsCol = list(map(list, list(zip(*pdsRow))))[1]
    a2fCol = list(map(list, list(zip(*a2fRow))))[1]
    couplings = [enCol,[]]
    for i, pds in enumerate(pdsCol):
        if float(pds) == 0.0:
            couplings[1].append(0.0)
        else:
            couplings[1].append( float(a2fCol[i]) / float(pds) )
    return couplings

def ifc2dym(file, input):
    from decimal import Decimal, ROUND_UP
    import math
    index = pp.Word(pp.nums)
    vector = pp.Group(pp.Word(pp.nums + ".+-E")*3)
    matrix = pp.Group(vector * 3)
    ifcblock = pp.Group(index*2 + vector + matrix)
    text = pp.OneOrMore(ifcblock)
    f = open(file, 'r')
    try:
        data = text.parseString(f.read()).asList()
    except pp.ParseException as pe:
        print(('Parsing Error using pyparsing: invalid input:', pe))
        sys.exit()
    f.close()
    natom = int(data[-1][0])
    # if natom is not int()  # check natom against System?
    nifc = int(data[-1][1])

    # Create and fill ifcInfo data structure
    ifcInfo = [[{} for b in range(nifc)] for a in range(natom)]
    for info in data:
        atomIndex = int(info[0]) - 1
        ifcIndex = int(info[1]) - 1 
        atomPos = list(map(float, info[2]))
        ifcMatrix= [list(map(float, row)) for row in info[3]]
        ifcInfo[atomIndex][ifcIndex] = {'pos':atomPos,'ifc':ifcMatrix} 

    # Convert VASP input strings into indexed numerical data
    acell = list(map(float, expandedList(input['acell'])))
    try:
        parsedMatrix = matrix.parseString(input['rprim']).asList()[0]
        rprim = [list(map(float, row)) for row in parsedMatrix]
    except pp.ParseException as pe:
        print(('Parsing Error using pyparsing: invalid input:', pe))
        sys.exit()
    znucl = list(map(int, expandedList(input['znucl'], length=natom)))
        
    ncluster = 300.0
    n = int(math.ceil((pow(ncluster/natom, 1.0/3.0) - 1.0) / 2.0))
    copyRange = list(range(-n, n+1))
    copyGrid = ((i,j,k) for i in copyRange for j in copyRange for k in copyRange)

    # Build cluster
    cluster = []
    for ijk in copyGrid:
        dx = [a*sum(i*c for i,c in zip(ijk,row)) for row, a in zip(rprim,acell)]
        for iAt in range(natom):
            clusterAtom = {'cellIndex':iAt}
            cellCoord = list(map(float, ifcInfo[iAt][0]['pos']))
            clusterAtom['coord'] = list(map(sum, list(zip(cellCoord, dx)))) 
            clusterAtom['dist'] = math.sqrt(sum(a*a for a in clusterAtom['coord']))
            clusterAtom['Z'] = znucl[iAt]
            clusterAtom['mass'] = amu[znucl[iAt]]
            if all(i is 0 for i in ijk):
                clusterAtom['central'] = True
            cluster.append(clusterAtom)
    ru8 = lambda x: Decimal(x).quantize(Decimal('0.12345678'),rounding=ROUND_UP)
    cluster.sort(key=lambda atom:(ru8(atom['dist']), 
                                  atom['cellIndex'],
                                  ru8(atom['coord'][0]),
                                  ru8(atom['coord'][1]),
                                  ru8(atom['coord'][2])
                                 ))
    ncluster = len(cluster)

    # Get distance/displacement for cluster pairs used for geometric matching
    clustPair = [[{} for j in range(ncluster)] for i in range(ncluster)]
    for i,j in ((i,j) for i in range(ncluster) for j in range(i,ncluster)):
        if i is j:
            clustPair[i][j]['dist'] = 0.0
            clustPair[i][j]['xi-xj'] = [0.0, 0.0, 0.0]
        else:
            displ = list(map(lambda x,y:x-y, cluster[i]['coord'], cluster[j]['coord']))
            clustPair[i][j]['dist'] = math.sqrt(sum(x*x for x in displ))
            clustPair[j][i]['dist'] = clustPair[i][j]['dist']
            clustPair[i][j]['xi-xj'] = displ
            clustPair[j][i]['xi-xj'] = [-x for x in displ]

    # Get distance/displacement for reference pairs used for geometric matching
    refPair = [[{} for b in range(nifc)] for a in range(natom)]
    for a,b in ((a,b) for a in range(natom) for b in range(nifc)):
        displ = list(map(lambda x,y:x-y, ifcInfo[a][0]['pos'], ifcInfo[a][b]['pos']))
        refPair[a][b]['dist'] = math.sqrt(sum(x*x for x in displ))
        refPair[a][b]['xa0-xab'] = displ

    # Start building a type-2 DYM dict
    dym = {'dymType':2} 
    
    # Section 1: Number of atoms in cluster
    dym['nAt'] = ncluster
    
    # Section 2: List of atomic numbers for each atom in cluster
    dym['atNums'] = [atom['Z'] for atom in cluster]

    # Section 3: List of atomic masses for each atom in cluster
    dym['atMasses'] = [atom['mass'] for atom in cluster]

    # Section 4: List of positions for each atom in cluster
    dym['atCoords'] = [atom['coord'] for atom in cluster]

    # Section 5: 3x3 block of dynamical matrix for each pair of atoms in cluster
    maxDiff = 0.001
    cutoffDist = 10.0
    areClose = lambda x,y: abs(x-y) <= maxDiff
    dym['dm'] = [[None for j in range(ncluster)] for i in range(ncluster)]
    for i,j in ((i,j) for i in range(ncluster) for j in range(ncluster)):
        dist = clustPair[i][j]['dist']
        displ = clustPair[i][j]['xi-xj']
        if dist > cutoffDist:
            dym['dm'][i][j] = [[0.0]*3]*3
        else:
            found = False
            # Loop over reference pairs of same original cell index
            atomCellIndex = cluster[i]['cellIndex']
            for b,refPairB in enumerate(refPair[atomCellIndex]):
                refDist = refPairB['dist']
                refDispl = refPairB['xa0-xab']
                # Found reference pair with same geometry
                if areClose(dist, refDist) and all(map(areClose, displ, refDispl)):
                    clustPair[i][j]['refifc'] = b
                    dym['dm'][i][j] = ifcInfo[atomCellIndex][b]['ifc']
                    found = True
                # At end of reference list with no match found
            if not found:
                raise Exception('Error: cannnot find reference for pair ' + str((i+1,j+1)))
    
    # Ordering of cluster atoms for printing purposes
    dym['printOrder'] = list(range(ncluster))

    # Section Type-2: Info for central cell atoms
    typat = list(map(int, expandedList(input['typat'], length=natom)))
    dym['nTypeAt'] = len(set(typat))
    dym['nCellAt'] = natom
    dym['centerAt'] = {z:[] for z in znucl}
    for i,atom in ((i,a) for i,a in enumerate(cluster) if 'central' in a):
        atomInfo = {'clustIndex':i}
        atomInfo.update({k:atom[k] for k in ('cellIndex','coord')})
        dym['centerAt'][atom['Z']].append(atomInfo)

    return dym

def nifc(acell, rprim, cutoff=10.0):
    import math
    dxBasis = [[r*a for r in col] for col, a in zip(list(zip(*rprim)),acell)]
    minDist = math.sqrt(min(sum(x*x for x in dx) for dx in dxBasis))
    nord = int(math.ceil(cutoff/minDist))
    iRange = list(range(-nord, nord+1))
    grid = ((i,j,k) for i in iRange for j in iRange for k in iRange)
    count = 0
    for ijk in grid:
        dx = [a*sum(i*c for i,c in zip(ijk,row)) for row, a in zip(rprim,acell)]
        dist = math.sqrt(sum(a*a for a in dx))
        if dist < cutoff:
            count += 1
    return count
    
def xasSnap(file):
    snapread = []
    with open(file) as f:
        s = f.readline()
        n = int(s)
        for line in f:
            snapread.append(line.split())
    a = [item for item in snapread if len(item) > 5]
    coordlist = [a[x:x+n] for x in range(0,len(a),n)]
    return coordlist
    
def mySequenceFor(output):
        if isinstance(output, list) and output and isinstance(output[0], str):
            key = strlistkey(output)
        elif isinstance(output, str):
            key = output
        else:
            raise TypeError('Output should be token or list of tokens')
        if key not in implemented:
            raise LookupError('Corvus cannot currently produce ' + key + ' using Nwchem')
        f = lambda subkey:implemented[key][subkey]
        if f('type') is 'Exchange':
            return Exchange(Nwchem, f('req'), f('out'), cost=f('cost'), desc=f('desc'))
    
#def loop()
    #exchange = self.Exchange(output)
    #xas_loop = self.Loop(exchange)
    #return xas_loop
    
def spinmult(input):

    import corvutils.corvus_misc_funcs as cu

#   Z = {}
#   Z['H'] = 1
#   Z['He'] = 2
#   Z['Li'] = 3
#   Z['Be'] = 4
#   Z['B'] = 5
#   Z['C'] = 6
#   Z['N'] = 7
#   Z['O'] = 8
#   Z['F'] = 9
#   Z['Ne'] = 10
#   Z['Na'] = 11
#   Z['Mg'] = 12
#   Z['Al'] = 13
#   Z['Si'] = 14
#   Z['P'] = 15
#   Z['S'] = 16
#   Z['Cl'] = 17
#   Z['Ar'] = 18
#   Z['K'] = 19
#   Z['Ca'] = 20
    coord = input['cell_struc_xyz_red']
    atoms = [item[0] for item in coord]
#   n = {k:Z[k] for k in atoms if k in Z}
    n = {k:cu.AtSym2AtNum(k) for k in atoms}
# Debug
#   print atoms
#   print n
#   print np
#   sys.exit()
    cnt = dict((x,atoms.count(x)) for x in set(atoms))
    totelect = sum(n[k]*cnt[k] for k in n)
    if totelect % 2 == 0:
        mult = 1
    elif totelect % 2 == 1:
        mult = 2
    return(str(mult))
  
    
# Atomic mass table from ABINIT
amu = {}
amu.update({  1:1.00794,    2:4.002602,   3:6.941,      4:9.012182,    5:10.811})
amu.update({  6:12.011,     7:14.00674,   8:15.9994,    9:18.9984032, 10:20.1797})
amu.update({ 11:22.989768, 12:24.3050,   13:26.981539, 14:28.0855,    15:30.973762})
amu.update({ 16:32.066,    17:35.4527,   18:39.948,    19:39.0983,    20:40.078})
amu.update({ 21:44.955910, 22:47.88,     23:50.9415,   24:51.9961,    25:54.93805})
amu.update({ 26:55.847,    27:58.93320,  28:58.69,     29:63.546,     30:65.39})
amu.update({ 31:69.723,    32:72.61,     33:74.92159,  34:78.96,      35:79.904})
amu.update({ 36:83.80,     37:85.4678,   38:87.62,     39:88.90585,   40:91.224})
amu.update({ 41:92.90638,  42:95.94,     43:98.9062,   44:101.07,     45:102.9055})
amu.update({ 46:106.42,    47:107.8682,  48:112.411,   49:114.82,     50:118.710})
amu.update({ 51:121.753,   52:127.60,    53:126.90447, 54:131.29,     55:132.90543})
amu.update({ 56:137.327,   57:138.9055,  58:140.115,   59:140.90765,  60:144.24})
amu.update({ 61:147.91,    62:150.36,    63:151.965,   64:157.25,     65:158.92534})
amu.update({ 66:162.50,    67:164.93032, 68:167.26,    69:168.93421,  70:173.04})
amu.update({ 71:174.967,   72:178.49,    73:180.9479,  74:183.85,     75:186.207})
amu.update({ 76:190.2,     77:192.22,    78:195.08,    79:196.96654,  80:200.59})
amu.update({ 81:204.3833,  82:207.2,     83:208.98037, 84:209.0,      85:210.0})
amu.update({ 86:222.0,     87:223.0,     88:226.0254,  89:230.0,      90:232.0381})
amu.update({ 91:231.0359,  92:238.0289,  93:237.0482,  94:242.0,      95:243.0})
amu.update({ 96:247.0,     97:247.0,     98:249.0,     99:254.0,     100:253.0})
amu.update({101:256.0,    102:254.0,    103:257.0,    104:260.0})

# Modified by FDV:
# Getting the generateInput and translateOutput methods out of the class.
# This makes them visible to all the other methods and usable within the handler

def generateInput(config, input, output):

# Debug: FDV
#  pp_debug.pprint(input)
#  pp_debug.pprint(output)
#  sys.exit()

# Start with the preparation of the POTCAR since that is needed for all
# calculations

# In VASP the order of the content of the POTCAR and that of the POSCAR are
# related, so we have to be careful on how we order things and we need to
# remember that order for future calculations.

# First we parse cell_red_xyz to find the atom labels

# Debug: FDV
   print('random stuff')
   print(input)
#  sys.exit()
        
        # Switch based on what output is required
   files = []
   dir = config['xcDir']
   ppFilelist = pp.OneOrMore(pp.Word(pp.printables)) 
#  pseudos = ppFilelist.parseString(input['pspfiles']).asList()
   pseudos = []

# Lattice optimization
   if set(output.keys()) == set(['aopt']):
      filebase = '0.nwchem.opt'
      pathbase = os.path.join(dir, filebase)
      nwchemFileList = nwchemFiles(pathbase, pseudos)
      writeList(nwchemFileList, pathbase + '.files')
      writeDict(optInput(input), pathbase + '.in')
      files.append(filebase + '.in')
      print((pathbase + '.in'))

#  Dynamical matrix
   if set(output.keys()) == set(['dynmat']):
      filebase = '0.nwchem.opt'
      pathbase = os.path.join(dir, filebase)
      nwchemFileList = nwchemFiles(pathbase, pseudos)
      writeList(nwchemFileList, pathbase + '.files')
# Debug
      print('\n optInput(input) \n')
      pp_debug.pprint(optInput(input))
#     sys.exit()
      writeDict(dynmatInput(input), pathbase + '.in')
      files.append(filebase + '.in')
# Debug
#     print pathbase + '.in'
#     sys.exit()

# Structure optimization + dynamical matrix
   if set(output.keys()) == set(['opt_dynmat']):
      filebase = '0.nwchem.opt'
      pathbase = os.path.join(dir, filebase)
      nwchemFileList = nwchemFiles(pathbase, pseudos)
      writeList(nwchemFileList, pathbase + '.files')
# Debug
      print('\n optInput(input) \n')
      pp_debug.pprint(optInput(input))
#     sys.exit()
      writeDict(opt_dynmatInput(input), pathbase + '.in')
      files.append(filebase + '.in')
# Debug
#     print pathbase + '.in'
#     sys.exit()

# Symmetry-reduced q-grid for use with ANADDB 
   elif set(output.keys()) == set(['qptlist']):
      filebase = '0.abinit.qgrid'
      pathbase = os.path.join(dir, filebase)
      abinitFileList = abinitFiles(pathbase, pseudos)
      writeList(abinitFileList, pathbase + '.files')
      writeDict(qgridInput(input), pathbase + '.in')
      files.append(filebase + '.files')

# Internal energy
   elif set(output.keys()) ==set(['ene_int']):
      filebase = '1.nwchem.U'
      pathbase = os.path.join(dir, filebase)
      nwchemFileList = nwchemFiles(pathbase, pseudos)
      writeList(nwchemFileList, pathbase + '.files')
      writeDict(UInput(input), pathbase + '.in') 
      files.append(filebase + '.in')
      print((pathbase + '.in'))

# Electron-phonon calculation using ANADDB
   elif set(output.keys()).issubset(set(['pdos','a2f','a2','eint'])):
      filebase = '1.abinit.elphon'
      pathbase = os.path.join(dir, filebase)
      abFileList = abinitFiles(pathbase, pseudos)
      writeList(abFileList,  pathbase + '.files')
      writeDict(elphonInput(input), pathbase + '.in')
      dsPrefix = abFileList[3] + '_DS'
      files.append(filebase + '.files')
      
      filebase = '2.mrgddb.elphon'
      pathbase = os.path.join(dir, filebase)
      qptRange = list(range(3, 3 + len(input['qptlist'])))
      ddbFileList = ddbFiles(pathbase, dsPrefix, qptRange)
      writeList(ddbFileList, pathbase + '.files')
      files.append(filebase + '.files')

      filebase = '3.mrggkk.elphon'
      pathbase = os.path.join(dir, filebase)
      atdirRange = list(range(1, 1 + 3*int(input['natom'])))
      gkkFileList = gkkFiles(pathbase, dsPrefix, qptRange, atdirRange)
      writeList(gkkFileList, pathbase + '.files')
      files.append(filebase + '.files')
      
      filebase = '4.anaddb.ifc'
      pathbase = os.path.join(dir, filebase)
      anaFileList = anaddbFiles(pathbase, ddbFileList[0], gkkFileList[0]) 
      writeList(anaFileList, pathbase + '.files')
      writeDict(ifcInput(input), pathbase + '.in')
      files.append(filebase + '.files')
    
# Quantum Molecular Dynamics   
   elif set(output.keys()) == set(['qmd']):
      filebase = '1.nwchem.qmd'
      pathbase = os.path.join(dir, filebase)
      nwchemFileList = nwchemFiles(pathbase, pseudos)
      writeList(nwchemFileList, pathbase + '.files')
      writeDict(qmdInput(input), pathbase + '.in') 
      files.append(filebase + '.in')
      print((pathbase + '.in'))
    
# qmd with optimized geometry
   elif set(output.keys()) == set(['opt_qmd']):
      filebase = '1.nwchem.optqmd'
      pathbase = os.path.join(dir, filebase)
      nwchemFileList = nwchemFiles(pathbase, pseudos)
      writeList(nwchemFileList, pathbase + '.files')
      writeDict(qmdInput(input), pathbase + '.in') 
      files.append(filebase + '.in')
      print((pathbase + '.in'))
    
# x-ray spectrum(single shot)
   elif set(output.keys()) == set(['xas']):
      filebase = '1.nwchem.xas'
      pathbase = os.path.join(dir, filebase)
      nwchemFileList = nwchemFiles(pathbase, pseudos)
      writeList(nwchemFileList, pathbase + '.files')
      writeDict(xasInput(input), pathbase + '.in') 
      files.append(filebase + '.in')
      print((pathbase + '.in'))
    
# xas with optimized geometry
   elif set(output.keys()) == set(['opt_xas']):
      filebase = '1.nwchem.optxas'
      pathbase = os.path.join(dir, filebase)
      nwchemFileList = nwchemFiles(pathbase, pseudos)
      writeList(nwchemFileList, pathbase + '.files')
      writeDict(xasInput(input), pathbase + '.in') 
      files.append(filebase + '.in')
      print((pathbase + '.in'  ))

# average x-ray spectrum
   elif set(output.keys()).issubset(set(['xas_avg','user_xasavg'])):
      filebase = 'snapshots.xyz'
      pathbase = os.path.join(config['cwd'], filebase)
      coords = xasSnap(pathbase)
      my_exchange = mySequenceFor('xas')
      my_loop = Loop(my_exchange, 'cell_struc_xyz_red', None, coords)
      my_loop.go(config,input)
#     print my_exchange
#     print config['cwd']
#     print coords
#     print pathbase
#     print 'not implemented'
#     print input['xas-grid']

# optimized average x-ray spectrum   
   elif set(output.keys()) == set(['opt_xasavg']):
      filebase = 'snapshots.xyz'
      pathbase = os.path.join(config['cwd'], filebase)
      coords = xasSnap(pathbase)
      my_exchange = mySequenceFor('xas')
      my_loop = Loop(my_exchange, 'cell_struc_xyz_red', None, coords)
      my_loop.go(config,input)

# Return ordered list of *.files
   return files

def translateOutput(config, input, output): 

   dir = config['xcDir']
   for target in output:
      if target == 'aopt':
         file = os.path.join(dir, '0.nwchem.opt.out')
         output[target] = getacell(file)
         input['cell_struc_xyz_red'] = optqmdcell(file)
      elif target == 'dynmat':
         file_out  = os.path.join(dir, '0.nwchem.opt.out')
         file_hess = os.path.join(dir, '0.nwchem.opt.hess')
# NOTE: This should probably be changed to output?
#        input['cell_struc_xyz_red'] = optqmdcell(file_out)
         output[target] = getdynmat(input,file_hess)
      elif target == 'opt_dynmat':
         file_out  = os.path.join(dir, '0.nwchem.opt.out')
         file_hess = os.path.join(dir, '0.nwchem.opt.hess')
# NOTE: This should probably be changed to output?
         input['cell_struc_xyz_red'] = optqmdcell(file_out)
         output[target] = getdynmat(input,file_hess)
      elif target == 'acell-grid':
         num = int(check(input, 'nagrid', default=9))
         output[target] = getagrid(input['acell'], num)
      elif target == 'qptlist':
         file = os.path.join(dir, '0.abinit.qgrid.out')
         output[target] = kgrid2qgrid(file)
      elif target == 'ene_int':
         if set(output) == {'ene_int'}:
            file = os.path.join(dir, '1.nwchem.U.out')
            output[target] = getU(file)
         elif len(set(output).intersection({'dynmat','pdos','a2','a2f'})) > 0:
            file = os.path.join(dir, '1.abinit.elphon.out')
            output[target] = getU(file, addDS='1')
      elif target == 'pdos':
         file = os.path.join(dir, '4.anaddb.ifc.out_ep_PDS')
         output[target] = anaddb2cols(file)
      elif target == 'a2f':
         file = os.path.join(dir, '4.anaddb.ifc.out_ep_A2F')
         output[target] = anaddb2cols(fisle)
      elif target == 'a2':
         file = os.path.join(dir, '4.anaddb.ifc.out_ep_')
         output[target] = eli2couplings(file) 
      elif target == 'dynmat':
         file = os.path.join(dir, 'ifcinfo.out')
         output[target] = ifc2dym(file, input)
      elif target == 'qmd':
         file = os.path.join(dir, '1.nwchem.qmd.xyz')
         output[target] = snapshots(file,input)
         f = open('snapshots.xyz','w')
         snaps = snapshots(file,input)
         f.write(snaps)
         f.close()
      elif target == 'opt_qmd':
         file = os.path.join(dir, '1.nwchem.optqmd.xyz')
         output[target] = snapshots(file,input)
         f = open('snapshots.xyz','w')
         snaps = snapshots(file,input)
         f.write(snaps)
         f.close()
      elif target == 'xas':
         file = os.path.join(dir, '1.nwchem.xas.out')
         output[target] = getSpect(file)
         eV, dipole = getSpect(file)
         file2 = os.path.join(dir, 'spectrum.out')
         f = open(file2,'w')
         for index in range(len(eV)):
            f.write(eV[index] + ' ' + dipole[index] + '\n')
         f.close()
      elif target == 'opt_xas':
         file = os.path.join(dir, '1.nwchem.optxas.out')
         output[target] = getSpect(file)
         eV, dipole = getSpect(file)
         file2 = os.path.join(dir, 'spectrum.out')
         f = open(file2,'w')
         for index in range(len(eV)):
            f.write(eV[index] + ' ' + dipole[index] + '\n')
         f.close()
      elif target == 'xas_avg':
         file = os.path.join(dir, '1.nwchem.xas_avg.out')
         output[target] = avgSpect(input)
         energy, avgstr = avgSpect(input)
         f = open('avgSpect.out','w')
         for index in range(len(energy)):
            f.write(energy[index] + ' ' + avgstr[index] + '\n')
         f.close()
      elif target == 'user_xasavg':
         file = os.path.join(dir, '1.nwchem.xas_avg.out')
         output[target] = avgSpect(input)
         energy, avgstr = avgSpect(input)
         f = open('avgSpect.out','w')
         for index in range(len(energy)):
            f.write(energy[index] + ' ' + avgstr[index] + '\n')
         f.close() 
      elif target == 'opt_xasavg':
         file = os.path.join(dir, '1.nwchem.opt_xasavg.out')
         output[target] = avgSpect(input)
         energy, avgstr = avgSpect(input)
         f = open('avgSpect.out','w')
         for index in range(len(energy)):
            f.write(energy[index] + ' ' + avgstr[index] + '\n')
         f.close()

