#!/usr/bin/env python

# Debug: FDV
import pprint
pp_debug = pprint.PrettyPrinter(indent=4)

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

# Convert a string with certain content into a python boolean
def boolean(string):

# Remove all leading and trailing blanks and turn into lowercase, to ensure
# things are ok and make them simpler
  string = string.strip().lower()

# Define the possible values of a true boolean
  True_Values  = ['1','t','true','.true.','y']

# Define the possible values of a false boolean
  False_Values = ['0','f','false','.false.','n']

  if string in True_Values and string not in False_Values:
    return True
  if string not in True_Values and string in False_Values:
    return False
  if string not in True_Values and string not in False_Values:
    raise ValueError('Error in corvus_misc_funcs.boolean: Unrecognized value ' + string)

# Try to automatically determine a value's type and convert
def autotype(string):

# Remove all leading and trailing blanks and turn into lowercase, to ensure
# things are ok and make them simpler
  string = string.strip().lower()

# Debug:FDV
  value = string

# There is probably a better way to do this is python, but I hack a quick
# solution here. Will return to this later.

# Check if the string is an integer
  int_chars = set('+-0123456789')
  if set(string).difference(int_chars) == set([]):
    try:
      value = int(string)
    except:
      pass
    return value

  flt_chars = set('+-0123456789.e')
  if set(string).difference(flt_chars) == set([]):
    try:
      value = float(string)
    except:
      pass
    return value

  return value

# Periodic Table functions
def AtLbl2AtSym(Label):

  Symb1 = Label[0:2]
  if Symb1 in atomicSymbols:
    Symb = Symb1
  else:
    Symb2 = Label[0:1]
    if Symb2 in atomicSymbols:
      Symb = Symb2
    else:
      raise ValueError('Error in corvus_misc_funcs.AtLbl2AtSym: Unrecognized element label ' + Label)

  return Symb

def AtSym2AtNum(Symb):

  AtNum = atomicSymbols.index(Symb)+1

  return AtNum

# Convert a list of list of strings into a single string where the internal
# lists are concatenated with spaces and then those are concatenated with
# newlines.
def StrList2Str(Line_List):

  Str = ''

  for Line in Line_List:
    Str = Str + ' '.join(Line) + '\n'

  return Str.rstrip()

# This helper function simplifies the calls to the external subprocesses and
# handles the printing to stdout and stderr depending on the verbosity
# requested
def subp(config,exe,cwdir,out,err):

  import sys, subprocess

  if config['verbose'] > 0:
    p = subprocess.Popen(exe, cwd=cwdir,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for line in p.stdout:
      sys.stdout.write(line)
      out.write(line)
    for line in p.stderr:
      sys.stdout.write(line)
      err.write(line)
  else:
    p = subprocess.Popen(exe, cwd=cwdir, stdout=out, stderr=err)
  p.wait()

# Helper function to create a gnupllot script that plots data
# NOTE FDV: This current implementation is a hack. We need to improve it to
# make it more flexible. For now we assume the data file has comments defined
# by "#" and then a series of lines of data with the same number of columns.
def Make_Gnuplot(filename,Ctrl):

  import sys, os
  import corvutils.pyparsing as pp

# Debug: FDV
# print filename
# pp_debug.pprint(Ctrl)

# NOTE FDV: This is not the best way to parse this data, but I want to get
# something going quickly. Will need to make a bit prettier later once I learn
# more about pyparsing
  Comm_Chars ='#'
  commentGrammar = pp.Word(Comm_Chars) + pp.restOfLine
  num = pp.Word(pp.nums + ".+-ED")
  row = pp.OneOrMore(num)
  input = open(filename, 'r').read()
  input = commentGrammar.suppress().transformString(input)
  data = []
  for line in input.split('\n'):
    if line != '':
      data.append(row.parseString(line).asList())

# Sanity checks
  nCol = [ len(row) for row in data ]
  if min(nCol) != max(nCol):
    print "Error in Make_Gnuplot: Variable number of columns"
    sys.exit()
  nCol = min(nCol)

# NOTE FDV:
# Here I should make sure that I print to the right directory. For now I will
# assume (as it is usually correct) that the cwd is the Corvus running
# directory. This can be checked by uncommenting the following two lines
# cwd = os.getcwd()
# print cwd

  for key, value in Ctrl.iteritems():
    xcol = value[0][0]
    ycol = value[0][1]
    xlab = value[1][0]
    ylab = value[1][1]
    if xcol < 1 or xcol > nCol:
      print 'Error in Make_Gnuplot: Wrong x column index for data ', key
      sys.exit()
    if ycol < 1 or ycol > nCol:
      print 'Error in Make_Gnuplot: Wrong x column index for data ', key
      sys.exit()
    file = open(key+".gplot","w")
#  Write a simple file header with some info.
# NOTE FDV: We should add more info here
    file.write('# Corvus run output: '+key+'\n')
    file.write('set grid\n')
    file.write('set xlabel "'+xlab+'"'+'\n')
    file.write('set ylabel "'+ylab+'"'+'\n')
    file.write('# Data: '+xlab+' '+ylab+'\n')
    file.write('plot "-" with lines title '+'"'+ylab+'"'+'\n')
    for line in data:
      file.write(line[xcol-1]+' '+line[ycol-1]+'\n')
    file.write('e'+'\n')
    file.write('pause -1\n')
    file.close()

# Debug: FDV
  sys.exit()

# Should add testing when this is executed
def main():
  pass
