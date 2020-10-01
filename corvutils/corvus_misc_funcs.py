#!/usr/bin/env python

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

# Should add testing when this is executed
def main():
  pass
