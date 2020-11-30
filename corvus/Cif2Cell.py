from corvus.structures import Handler, Exchange, Loop, Update
import corvutils.pyparsing as pp
import os, sys, subprocess, shutil, resource
import re
# Debug: FDV
import pprint
# The below is from the pycifrw package. It comes with 
# cif2cell, so it should work if cif2cell is installed.
from CifFile import ReadCif
# This one is from the cif2cell package. It allows calculations of 
# cell properties from cif input.
from cif2cell.uctools import *

pp_debug = pprint.PrettyPrinter(indent=4)


# Define dictionary of implemented calculations
implemented = {}
strlistkey = lambda L:','.join(sorted(L))
subs = lambda L:[{L[j] for j in range(len(L)) if 1<<j&k} for k in range(1,1<<len(L))]
for s in subs(['cell_vectors', 'cell_struct_xyz_red', 'cell_scaling_iso', 'cell_scaling_abc', 'number_density']):
    key = strlistkey(s)
    autodesc = 'Get ' + ', '.join(s) + ' using cif2cell'
    cost = 10
    implemented[key] = {'type':'Exchange','out':list(s),'req':['cif_input'],
                        'desc':autodesc,'cost':cost}

implemented['cell_structure'] = {'type':'Exchange','out':['cell_structure'],'cost':0,
                        'req':['cell_vectors','cell_struct_xyz_red','cell_scaling_iso','cell_scaling_abc'],'desc':'Calculate cell structure from cif file using cif2cell.'}


#implemented['cell_structure'] = {'type':'Exchange','out':['cell_structure','cell_struc_xyz','cell_scaling_abc','cell_scaling_iso'],'cost':0,
#                        'req':['cif_input'],'desc':'Calculate cell structure from cif file using cif2cell.'}


#implemented['cluster'] = 


class Cif2Cell(Handler):
    def __str__(self):
        return 'cif2cell Handler'

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
            unresolved = {o for o in output if not Cif2Cell.canProduce(o)}
            canProduce = (o for o in output if Cif2Cell.canProduce(o))
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
        if f('type') is 'Exchange':
            return Exchange(Cif2Cell, f('req'), f('out'), cost=f('cost'), desc=f('desc'))

    @staticmethod
    def prep(config):
        subdir = config['pathprefix'] + str(config['xcIndex']) + '_CIF2CELL'
        xcDir = os.path.join(config['cwd'], subdir)
        # Make new output directory if if doesn't exist
        if not os.path.exists(xcDir):
            os.mkdir(xcDir)
        # Store current Exchange directory in configuration
        config['xcDir'] = xcDir

    #@staticmethod
    #def setDefaults(input,target):

    @staticmethod
    def run(config, input, output):
        
        dir = config['xcDir']

        # Copy cif2cell related input to cif2cellInput here. Later we will be overriding some settings,
        # so we want to keep the original input intact.
        cif2cellInput = {key:input[key] for key in input if key.startswith('cif2cell.')}

        # Generate any data that is needed from generic input and populate cif2cellInput with
        # global data (needed for all feff runs.)
        cif2cellInput['cif2cell.cif_input'] = input['cif_input']

        # Set executable directory for this exchange
        cif2celldir = config['cif2cell']

        # Set input file - cif2cell takes input from command line, so no input file, however, maybe we can take commands from input file?
        #inpf = os.path.join(dir, 'cif2cell.in')

        # Loop over targets in output. Not sure if there will ever be more than one output target here.
        if set(output.keys()).issubset(set(['cell_vectors', 'cell_struct_xyz_red', 'cell_scaling_iso', 'cell_scaling_abc','number_density'])):
            # Set output and error files
            with open(os.path.join(dir, 'corvus.CIF2CELL.stdout'), 'wb') as out, open(os.path.join(dir, 'corvus.CIF2CELL.stderr'), 'wb') as err:
                # Copy necessary files to dir
                cif_file=cif2cellInput['cif2cell.cif_input'][0][0]
                shutil.copy(cif_file,dir)

                # Loop over executables: This is specific to feff. Other codes
                # will more likely have only one executable.
                executables = ['cif2cell']

                args=[cif_file, '-p', 'cif', '-o', 'out.cif']
                iExec = 0
                for executable in executables:
                    runExecutable(cif2celldir,dir,executable,args,out,err)

                
                # For now, I am only passing the directory.
                print('Setting output')
                outCif = os.path.join(dir,'out.cif')
                cifFile=ReadCif(outCif)
                # cifFile now contains an object that acts like a dictionary
                # of dictionaries, with the outer keys the data for each 
                # structure listed in the cif, next the normal cif data, like
                # '_cell_angle_alpha', and '_atom_site_fract_[xyz]'.
                # Note that as indicated above, the cif can hold multiple 
                # structures. For now I will assume that there is only one. 
                cif_dict = cifFile[list(cifFile.keys())[0]]
                                 
                # Now lets make all of the data structures that we want out of
                # this. We are using the CellData structure from cif2cell's
                # uctools package. Right now, I have cif2cell write another
                # cif file first, then we read it in and analyze it. This is
                # because we don't want to reproduce all of the sanity check
                # that cif2cell already has in it. 
                cell_data = CellData()
                # Fill the cell data object from the cif file.
                cell_data.getFromCIF(cif_dict)

                # Calculate cell structure for the primitive cell for now. Supercells are
                # represented as P1, so this will not harm that.
                cell_data.primitive()
                
                # The cell_data structure now has all cell data that we need 
                # in it, or a method to get that data. 
                if 'cell_vectors' in output:
                    output['cell_vectors'] = cell_data.latticevectors

                if 'cell_struct_xyz_red' in output:
                    xred = []
                    for atom in cell_data.atomdata:
                        if atom:
                            xred = xred + [atom[0].position]

                    output['cell_struct_xyz_red'] = xred

                if 'cell_scaling_iso' in output:
                    output['cell_scaling_iso'] = [[cell_data.lengthscale]]
                    
                if 'cell_scaling_abc' in output:
                    # For now I will output cell_data_abc as 1, since cif2cell seems to put things into lengthscale
                    # and a,b,c are for convenience.
                    output['cell_scaling_abc'] = [[1.0,1.0,1.0]]

        elif 'cell_structure' in output:
            # Return path to cif file.
            output['cell_structure'] = [['Working']]



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
        
    
    # Print cif2cell input file
    writeList(lines, inpfile)
    
def getInpLines(input,token):
    lines=[]
    block=False
    endblock=' '
    key = token[len('cif2cell.'):]
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
    print(('Running exectuable: ' + executable))
    # Modified by FDV:
    # Adding the / to make the config more generic
    # Modified by JJK to use os.path.join (even safer than above).
    execList = [os.path.join(execDir,executable)] + args
    print(execList)
    result = subprocess.run(execList, cwd=workDir, stdout=out, stderr=err)
    #print("stdout:", result.stdout)
    #out.write(result.stdout)
    #print("stderr:", result.stderr)
    #out.write(result.stderr)
#    p = subprocess.Popen(execList, bufsize=0, cwd=workDir, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#    while True:
#        pout = p.stdout.readline()
#        if pout == '' and p.poll() is not None:
#            break
#        if pout:
#            print((pout.strip()))
#            out.write(pout)
#
#    while True:
#        perr = p.stderr.readline()
#        if perr == '' and p.poll() is not None:
#            break
#        if perr:
#            print('###################################################')
#            print('###################################################')
#            print((perr.strip()))
#            print('###################################################')
#            print('###################################################')
#            err.write(perr)

            
#    print('waiting')
#    p.wait()
#    print('done waiting')
    
    


    
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

def writeXANESInput(input, cif2cellinp='cif2cell.in'):
    
    lines = []
    #setInput(input,'feff.print',[[5,0,0,0,0,0]],Force=True)

    writeInput(input,cif2cellinp)


