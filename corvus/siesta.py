from corvus.structures import Handler, Exchange, Loop, Update
import corvutils.pyparsing as pp
import os, sys, subprocess, shutil #, resource
import numpy as np
import re
# Debug: FDV
import pprint
pp_debug = pprint.PrettyPrinter(indent=4)

# Define dictionary of implemented calculations
implemented = {}
strlistkey = lambda L:','.join(sorted(L))

implemented['siestaCoreResponse'] = {'type':'Exchange','out':['siestaCoreResponse'],'cost':3,
                        'req':['supercell'],'desc':'Calculate XPS using siesta.'}

implemented['siestaFourierTransform'] = {'type':'Exchange','out':['siestaFourierTransform'],'cost':1,
            'req':['siestaCoreResponse'],'desc':'Calculated Desity of States'}
implemented['siestaBetaofOmega'] = {'type':'Exchange','out':['siestaBetaofOmega'],'cost':1,
			'req':['siestaFourierTransform'],'desc':'Calculate Beta of Omega'}
implemented['siestaDOS'] = {'type':'Exchange','out':['siestaDOS'],'cost':1,
            'req':['siestaCoreResponse'],'desc':'Calculated Desity of States'}

class Siesta(Handler):
    def __str__(self):
        return 'SIESTA Handler'

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
            unresolved = {o for o in output if not Siesta.canProduce(o)}
            canProduce = (o for o in output if Siesta.canProduce(o))
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
            return Exchange(Siesta, f('req'), f('out'), cost=f('cost'), desc=f('desc'))

    @staticmethod
    def prep(config):
        subdir = config['pathprefix'] + str(config['xcIndex']) + '_SIESTA'
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
        

        # Set directory to siesta executables.
        # Debug: FDV
        #       pp_debug.pprint(config)
        siestadir = config['siesta']
        # Debug: FDV
        #       sys.exit()

        # Copy siesta related input to siestaInput here. Later we will be overriding some settings,
        # so we want to keep the original input intact.
        siestaInput = {key:input[key] for key in input if key.startswith('siesta.')}
       
        siestaInput['siesta.td.numberoftimesteps'] = input.get('siesta.td.numberoftimesteps',[[500]])
        siestaInput['siesta.td.timestep'] = input.get('siesta.td.timestep',[[0.5]])
        siestaInput['siesta.td.shapeofefield'] = input.get('siesta.td.Shapeofefield',[['core']])
        siestaInput['siesta.td.coreperturbationcharge'] = input.get('siesta.td.coreperturbationcharge',[[0.02]])
        siestaInput['siesta.td.mxpc'] = input.get('siesta.td.mxpc',[[2]])
        siestaInput['siesta.solutionmethod'] = [[ 'diagon' ]]
        # Now get siesta inputs from general input.
        if 'absorbing_atom' in input:
            siestaInput['siesta.td.coreexcitedAtom'] = input['absorbing_atom']
        if 'siesta.block.atomiccoordinatesandatomicspecies' not in input:
            siestaInput['siesta.atomiccoordinatesformat'] = [['Fractional']]
            siestaInput['siesta.numberofatoms'] = input['number_of_atoms']
            siestaInput['siesta.numberofspecies'] = input['number_of_species']
            species_label = []
            for i,s in enumerate(input['species']):
                species_label = species_label + [[i+1, s[0], s[1]]]


            sites = []
            for j,site in enumerate(input['cell_struc_xyz_red']):
                for i,s in enumerate(input['species']):
                    if site[0] == s[1]:
                        sites = sites + [[site[1],site[2],site[3],i+1]]

            siestaInput['siesta.block.chemicalspecieslabel'] = species_label
            siestaInput['siesta.latticeconstant'] = [[1.0, 'Ang']]
            siestaInput['siesta.block.atomiccoordinatesandatomicspecies'] = sites
            siestaInput['siesta.block.latticeparameters'] = [input['cell_scaling_abc'][0] + input['cell_angles_abg'][0]]
            

        # Generate any data that is needed from generic input and populate siestaInput with
        # global data (needed for all feff runs.)

        # Set directory for this exchange
        dir = config['xcDir']

        
        # Set input file
        inpf = os.path.join(dir, 'input.fdf')
        # Define EIG file
        EIG = os.path.join(dir, 'siesta.EIG')

        # Loop over targets in output. Not sure if there will ever be more than one output target here.
        for target in output:
            if (target == 'siestaFourierTransform'):
                print("Calculating siestaFourierTransform")
                
                # Remove contribution from zero frequency	
                CoreResponse = np.asarray([input['siestaCoreResponse'][0],input['siestaCoreResponse'][1]])
                print(CoreResponse[1])
                tmp_avg = np.average(CoreResponse[1])
                CoreResponse[1] = np.subtract(CoreResponse[1], tmp_avg)	
                
                # Broadening the coreresponse, default 0.5
                broadfactor = input['siesta.coreresponse.broadening'][0][0]
                print('broadfactor=',broadfactor)
                CoreResponse[1] = CoreResponse[1]*np.exp(-(broadfactor/27.21*CoreResponse[0]/0.02418884326505)**2)
                
                # Remove contribution from zero frequency again	
                tmp_avg = np.average(CoreResponse[1])
                CoreResponse[1] = np.subtract(CoreResponse[1], tmp_avg)	
                
                # Switch to Atomic Units
                CoreResponse[1] = CoreResponse[1]/27.21 

                # Calculate frequency array
                timestep = (CoreResponse[0][1]-CoreResponse[0][0])
                print(timestep)
                timestep = timestep/(0.02418884326505) #Switching from fs to au
                freq = np.fft.fftfreq(len(CoreResponse[1]), timestep)
                holder = []
                for i in freq:
                    if (i>=0):
                        holder.append(i) 
                freq = np.array(holder)
                #print(freq)

                # Implementing trapezoidal fourier transform
                f_w = np.zeros(len(freq))
                timegrid = np.asarray(CoreResponse[0])/(0.02418884326505)
                
                def Trapz_Int(function,dx):
                    return 0.5*dx*(2*np.sum(function)-function[0]-function[-1])

                for i in range(len(freq)):
                    Func = np.exp(2*np.pi*1j*freq[i]*timegrid)*CoreResponse[1]
                    f_w[i] = Trapz_Int(Func,(timegrid[1]-timegrid[0]))
                
                # Switching from Hartree Atomic Units to eV, multiply by 2pi for angular
                fourierTrapOutput = [np.multiply(freq,27.21*2*np.pi).tolist(),(f_w/(2*np.pi)).tolist()]
                
                # For graphing: Charles
                # np.savetxt("fourierTrapOutput.dat", np.array(fourierTrapOutput).T)

                outFile=os.path.join(dir,'siestaFourierTransformResult')
                #temp = np.array(fourierTrapOutput)
                #for i in range(len(temp[0])):
                #    temp[1][i] = temp[0][i]*temp[1][i]
                #np.savetxt("BetaRaw", temp.T)
                output[target] = fourierTrapOutput 

            elif (target == 'siestaBetaofOmega'):
                print("Calculating BetaofOmega")
                
                def writeBETA():
                    # Define Beta
                    Beta = np.array(input['siestaFourierTransform'])

                    # multiply by omega for Beta(omega)
                    for i in range(len(Beta[0])):
                        Beta[1][i] = Beta[0][i]*Beta[1][i]
                    
                    #Leave all broadening in time domian in Core Response
                    '''
                    def Broadening(arr, broad):
                        temp = []
                        stepsize = arr[0][1]-arr[0][0]
                        for i in arr[0]:
                            kernel = broad/(np.pi) * 1/((arr[0]-i)**2 + broad**2)
                            temp.append(sum(kernel*arr[1]) * stepsize)
                        arr[1] = np.array(temp)
                        return arr

                    Beta = Broadening(Beta, input['siesta.Beta.Broadening'][0][0])
                    #print(Beta[0][1]-Beta[0][0])
                    # For graphing: Charles
                    # np.savetxt("Beta_before.dat", np.array(Beta).T)
                    '''

                    # Flip across y-axis
                    Beta[0,:] = Beta[0,:]*(-1)
                    
                    # Pad with zeros
                    # Right now it just pads with 50 zeros on each side, should be updated to pad 20 eV
                    timestep = Beta[0][1]-Beta[0][0]
                    print("Beta Timestep = " + str(abs(timestep)))
                    back_hold = (Beta.min(1))[0]
                    front_temp = [[],[]]
                    back_temp = [[],[]]
                    x = 0
                    y = 0
                    z = 0
                    while(x<100):
                        front_temp[0].append(y-timestep)
                        y = y - timestep
                        front_temp[1].append(0)
                        
                        back_temp[0].append(back_hold+z+timestep)
                        z = z + timestep
                        back_temp[1].append(0)
                        
                        x = x+1
                    
                    front_temp = np.array(front_temp)
                    back_temp = np.array(back_temp)
                    
                    Beta = np.concatenate((Beta, front_temp), axis=1)
                    Beta = np.concatenate((back_temp, Beta), axis=1)
                    
                    # Sort Beta
                    Beta = Beta.T
                    Beta = Beta[np.argsort(Beta[:,0])]
                    
                    temp = Beta.T
                    for i in range(len(temp[0])):
                        if temp[1][i] < 0.0:
                            temp[1][i] = 0.0
                    
                    output[target]= temp.tolist()

                writeBETA()		
            
            elif (target == 'siestaDOS'):
                print("hi")
                '''
                print("Calculating Density of States")
               
                fermiE = 0
                def readEIG(filename):
                    doc = open(filename, "r")
                    holder = []
                    for i in doc:
                        for j in i.split():
                            holder.append(float(j))
                    fermiE = holder[0]
                    return np.asarray(holder[5:])-fermiE

                def makedos(arr, spacer):
                    roundnum = len(str(spacer))-2
                    arr = (np.asarray(np.round(arr,roundnum))).tolist()
                    temp = [[],[]]
                    temp[0].append(arr[0])
                    temp[1].append(1)
                    for i in range(len(arr)-1):
                        if(arr[i+1] == arr[i]):
                            temp[1][-1] = temp[1][-1]+1
                        else:
                            temp[0].append(arr[i+1])
                            temp[1].append(1)
                    
                    print(temp)
                    data = np.arange(-30,140,spacer)
                    data = np.asarray([data,np.zeros(np.shape(data))])
                    for i in range(len(temp[0])):
                        x = np.arange(-30,140,spacer)
                        arr = np.asarray([x,gaussian(x,temp[0][i],temp[1][i],0.5)])
                        data[1] = data[1] + arr[1]
                    
                    return data
                
                EIG = readEIG(EIG)
                EIG = makedos(EIG, spacer=0.01)
                EIG[1] = EIG[1]/(sum(EIG[1])*np.abs(EIG[0][0]-EIG[0][1])

                #output[target] = EIG.tolist()
                '''

            elif (target == 'siestaCoreResponse'):               
                # Set output and error files
                with open(os.path.join(dir, 'corvus.SIESTA.stdout'), 'w+') as out, open(os.path.join(dir, 'corvus.SIESTA.stderr'), 'w+') as err:
                    # Get pseudopotentials
                    
                    # Write input file for FEFF.
                    writeXPSInput(siestaInput,inpf)
                    
                    # Copy pseudos to dir
                    for file in os.listdir("."):
                        if file.endswith(".psf"):
                           shutil.copy(file,dir)

                    # Loop over executables: This is specific to feff. Other codes
                    # will more likely have only one executable.
                    executables = siestaInput.get('siesta.mpi.cmd',[['mpirun']])[0]
                    args = siestaInput.get('siesta.mpi.args',[['-n 4']])[0] + [os.path.join(siestadir,'siesta')]
                    iExec = 0
                    for executable in executables:

                        runExecutable('',dir,executable,args,out,err)

                
                # For now, I am only passing the directory.
                print('Setting output')
                outFile=os.path.join(dir,'coreresponse.vs.time')
                output[target] = np.loadtxt(outFile).T.tolist()



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
        
    
    # Print siesta input file
    writeList(lines, inpfile)
    
def getInpLines(input,token):
    lines=[]
    block=False
    endblock=' '
    key = token[len('siesta.'):]
    if key.startswith('block.'):
        block=True
        key = key[len('block.'):]
        lines = lines + ['%block ' + key.upper()]
        endblock = '%endblock ' + key.upper()

    if token in input:
        # If the first element is not a boolean, this contains values
        # to be stored after keyword.
        for element in input[token]: # Takes care of single and multi-line input.
            lines.append(' '.join([str(value) for value in element])) 

        if not block:
            lines = [key.upper() + ' ' + lines[0]]
          	
 

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
    print(('Running exectuable: ' + executable))
    # Modified by FDV:
    # Adding the / to make the config more generic
    # Modified by JJK to use os.path.join (even safer than above).
    execList = [os.path.join(execDir,executable)] + args
    inFile = open(os.path.join(workDir,'input.fdf'))
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
    #rc = p.poll() 
    #out.close()
    #err.close()
    
def readColumns(filename, columns=[1,2]):
    # Read file and clear out comments
    with open(filename, 'r') as file:
        cleanStr = file.read()
    comments = pp.ZeroOrMore(pp.pythonStyleComment).setParseAction(pp.replaceWith(''))
    try:
        cleanStr = comments.transformString(cleanStr)
    except pp.ParseException as pe:
        print(('Parsing Error using pyparsing: invalid input:', pe))
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
        print(('Parsing Error using pyparsing: invalid input:', pe))
        sys.exit()
    cols = list(map(list, list(zip(*data))))
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

def writeXPSInput(input, siestainp='input.fdf'):
    
    lines = []
    #setInput(input,'feff.print',[[5,0,0,0,0,0]],Force=True)

    writeInput(input,siestainp)


