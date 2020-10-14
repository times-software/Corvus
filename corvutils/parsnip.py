#!/usr/bin/env python

# Debug: FDV
import pprint
pp_debug = pprint.PrettyPrinter(indent=4)

import sys, getopt
from .pyparsing import *

#----------------------------------------------------------------------------
# Added by FDV
# See comments in formatInputs for an explanation of what these definitions do
import corvutils.corvus_misc_funcs as cu

# Allowed format type indicators for the parsnip.corvus.formats file
Str_Typ_Indicat = 'S'
Boo_Typ_Indicat = 'L'
Int_Typ_Indicat = 'I'
Flt_Typ_Indicat = 'F'
Aut_Typ_Indicat = 'A'
Pass_Indicat    = 'P'
RepData_Indicat = '...'

# Define the separator between the default and format part
DF_Sep = '|X|'

# Define the documentation marker character
Doc_Chars = '%'

# Define the comment characters
Comm_Chars ='#!'

# Define the conversion functions associated with each type, to avoid having
# to code a bunch of "ifs"
Typ_Conv_Funcs = {Str_Typ_Indicat:str,Boo_Typ_Indicat:cu.boolean,Int_Typ_Indicat:int,Flt_Typ_Indicat:float,Aut_Typ_Indicat:cu.autotype,RepData_Indicat:None,Pass_Indicat:None}

Allowed_Formats = set(Typ_Conv_Funcs.keys())
#----------------------------------------------------------------------------

# Default parsing modes replicate orignal parse script
defaultMode = ['WriteToFile','UseDefaults']

# Define core grammars
# 1. Token: the lookup key/label for the quantity
# 2. Body: text enclosed by curly braces
# 3. Stub: text (w/o whitespace) following curly braces
tokenGrammar = Word(alphanums + '.-_')
lbrace = Literal('{').suppress()
rbrace = Literal('}').suppress()
value = Word(printables.replace('{','').replace('}',''))
keepWS = Combine(value + ZeroOrMore(White() + value), adjacent=False)
bodyGrammar = (lbrace + Optional(keepWS, default=None) + rbrace) | value
stubGrammar = Word(printables)

# Ignored text
# Modified by FDV
# Changed to make the documentation tags be considered comments in the normal
# configuration process
#commentGrammar = Word(Comm_Chars) + restOfLine
commentGrammar = Word(Comm_Chars+Doc_Chars) + restOfLine
leadingWhiteGrammar = LineStart().leaveWhitespace() +  OneOrMore(White(' \t'))
commentGrammar.setParseAction(replaceWith(''))
leadingWhiteGrammar.setParseAction(replaceWith(''))

# Default mode recreates shell script parse
#def parse(configFile, formatFile, inputFile, mode=defaultMode):
def parse(configFile, inputFile, mode=defaultMode):
# Modified by FDV
    # Make sure we can open files
    with open(configFile, 'r') as config:
        configDict = readConfig(config).asDict()
    with open(inputFile, 'r') as input:
        inputDict = readInput(input).asDict()

# Testing: FDV
# The values that come out in the dictionaries have type pyparsing.ParseResults,
# which is probably ok, but can be pretty confusing when you print the data to
# explore its contents and it uses the pyparsing __repr__ method, which does
# funcky things to the data. Here I just test a reformat of the data to turn
# all the dictionary values into basic python lists.
# NOTE: Apparently inputDict already comes out in the right way. Will have to
# check readInput to see what is different from readConfig.
# J.K. In python 3, this seems to already be a list.
#    configDict = { key:value.asList() for key,value in configDict.items() }
#   inputDict = { key:value.asList() for key,value in inputDict.items() }

# Debug: FDV
#   pp_debug.pprint(configDict)
#   pp_debug.pprint(inputDict)
#   sys.exit()

# Now, before we proceed as before, we take the information we just read and
# split it into the dictionaries we originally had. In this way the code changes
# are minimal.
    
# Initialize the defaults and formats dicts
    defaultsDict = {}
    formatsDict  = {}

# Debug: FDV
# Here we split the data in the configDict into the parts that represent
# defaults and formats. This makes things easier later on allowing us to check
# the defaults in the same way we check the normal input.
# NOTE: The code is more complicated than it should be, but I don't want to
#       import the regular expressions module
    for key in configDict:
      val_default = ''
      val_format  = ''
      for line in configDict[key][0].split('\n'):
        line = line.strip()
        DF_Sep_Ind = line.find(DF_Sep)
        val_default = val_default+line[0:DF_Sep_Ind] + '\n'
        val_format  = val_format +line[DF_Sep_Ind+len(DF_Sep):len(line)] + '\n'
      val_default = val_default.strip()
      val_format  = val_format.strip()
# Debug: FDV
#     print [ val_default ]
#     print [ val_format ]
      defaultsDict.update({key:[val_default,configDict[key][1]]})
      formatsDict.update({key:[val_format,configDict[key][1]]})

# Debug: FDV
#   pp_debug.pprint(defaultsDict)
#   pp_debug.pprint(formatsDict)
#   sys.exit()

# At this point we have two dictionaries: one with the default values, and the
# other with the formats that can be used for both the defaults and the input.
# Here we use the formats for the defaults to amke sure the defaults are
# correctly defined and formatted.

# Loop over each key in and split the default values from the formats
#   for key in configDict:

# Ensure each line in the config value has the separator in it
#     if False in [ DF_Sep in line for line in configDict[key][0].split('\n') ]:
#       print 'Parsnip Error: Missing value/format separator for ', key
#       sys.exit()

#     lst = [ line.split(DF_Sep) for line in configDict[key][0].split('\n') ]

# Debug: FDV
#     print lst

# Now, before we put the default value back together, we need to manage an
# unintended consequence of mixing the default configuration input with the
# format input: The mixing makes it so that now there are always "defaults",
# some of which are null. Here we deal with those to eliminate them and pack
# things properly. It is a bit of a kludge, but it is needed to keep a nice
# and compact configuration file.
# The solution is to ensure that:
# 1) If a token has several lines of input (as implied by its associated format)
#    then each line must have non-null value.
# 2) To make a token required get its input from the input file, then each line
#    in its configuration must have a null format.
# 3) We don't accept mixed defaults with values in some lines and nulls in
#    others.

#     value_list = [line[0].strip() for line in lst]

# Check if all lines have a null default
#     if all([ elem is '' for elem in value_list ]):
#       value = ''
#     else:

# Check if all lines have an actual default
#       if all([ elem is not '' for elem in value_list ]):
#         value  = '\n'.join([line[0] for line in lst])
#       else:
# If we get here it means there was a problem
#         print 'Parsnip Error: Mixed null default for', key
#         sys.exit()

#     configDict[key][0]  = value

# Debug: FDV
#     print key, value

# The formats are easy to handle. Just join again and update the empty dict
#     format = ['\n'.join([line[1] for line in lst]),configDict[key][1]]
#     formatsDict.update({key:format})
    
# Debug: FDV
#   sys.exit()

# Debug: FDV
#   pp_debug.pprint(configDict)
#   pp_debug.pprint(defaultsDict)
#   pp_debug.pprint(formatsDict)
#   pp_debug.pprint(inputDict)
#   sys.exit()

# NOTE: FDV
# The check below is understandable, but I am not totally certain it is
# what we want. What is "required" will depend on what we are running, not
# just the difference between what is in the configuration file and the input.
# I am commenting this out.

# Check for any missing required inputs => error + exit
#   requiredTokens = [k for k in configDict.keys() if configDict.get(k)[0] is ''] 
#   missing = list(set(requiredTokens) - set(inputDict.keys()))
#   if missing and 'UseDefaults' in mode:
#       missing.sort()
#       output = 'Parsnip Error: missing required input(s): ' + str(missing)
#       raise LookupError(output)

# I am changing the original below, where unexpected inputs were just ignored
# and a warning issued. It is too risky to do it this way. Maybe people will
# make a typo on the input, Corvus will use the default instead and the user
# might not realize this.

# Debug: FDV
#   print 'configDict.keys()'
#   print configDict.keys()
#   print 'inputDict.keys()'
#   print inputDict.keys()
#   sys.exit()
# Check for any unexpected inputs => just issue a warning [FDV: not anymore]
    extra = list(set(inputDict.keys()) - set(configDict.keys()))
    if extra: 
#       print 'Parsnip Warning: unsupported option(s):', extra
        print(('Parsnip Error: unsupported option(s):', extra))
        sys.exit()

# NOTE FDV: Here we make a change and only initialize with the default if
# the default value isn't a ''. This means that if something is required and
# doesn't have a value the check has to be done later in the workflow.
# This is probably better since there are tokens that are not required by
# certain types of calculations but we can't really provide a default for.
# A typical example is coordinates, which can be reduced PBC, or just simply
# xyz, depending on the type of calculation.are tokens that are not required by
# certain types of calculations but we can't really provide a default for.

# Write out user input (or defaults) to associated files
    parsedInputs = {}
    for token in list(configDict.keys()):
        default,stub = defaultsDict[token]
        if inputDict.get(token) is None:
            if 'UseDefaults' in mode and default:

                parsedInputs[token] = default
                if 'WriteToFile' in mode:
                    with open(stub, 'w') as file:
                        file.write(default + '\n')
        else:
            parsedInputs[token] = inputDict[token]
            if 'WriteToFile' in mode:
                with open(stub, 'w'):
                    file.write(inputDict[token] + '\n')

# Debug:FDV
#   pp_debug.pprint(parsedInputs)
#   sys.exit()

# Added by FDV
    formatInputs(parsedInputs,formatsDict)

    return parsedInputs

# Added by FDV
def formatInputs(parsedInputs,formatsDict):

# Before returning the parsed input, we process it into a more useful format.
# In this way we don't have to be parsing and reparsing later on when we need
# to use the information.
# The format of each token is defined in the parsnip.corvus.formats file. For
# now we have the following simple definitions:
# 1:
#  T1
#    Token can only have a single value of type T1
# 2:
#  T1 T2
#    Token can have two values of types T1 and T2 (in a single line)
# 3:
#  T1 T2
#  ...
#    Token can have two values of type string per line, on one of more lines
# 4:
#  T1
#  T2 T3
#  ...
#    Token has value of type T1, followed by one or more lines of type T1 and T2
# etc.

# The general idea should be obvious. This format covers the most common
# formats for input tokens. We can try to make even more flexible later one.
# For instance, we could try to add more flexibility within a single input line.

# There are two special format indicators:
#  P, or pass:
#    This means to just leave the value as is, a string
#  A, or autotype:
#    This type indicates to try and find the type automatically, based on either
#    the default or the value entered by the user.

# The remaning values will be converted to only four types:
#  I type, or integer: 
#                      Example: 42
#  F type, or float:
#                     Example: -1.0e+10
#  L type, or boolean:
#                     Example: 0, 1, True, False, .true., .false., T, F, etc
#  S type, or string: everything else

# Make sure that formatsDict has the proper type indicators and format
  for (key,value) in list(formatsDict.items()):
    format = value[0].replace('\n',' ')

# Debug: FDV
#   print format

# Make sure that the first format indicator is not 'repeat data'
    if format.split()[0] == RepData_Indicat:
      print(('Parsnip Error: First format indicator in ' + \
        key + ' can not be \'' + RepData_Indicat + '\''))
      sys.exit() 

# Make sure that the format is allowed
    if not set(format.split()).issubset(Allowed_Formats):
      print(('Parsnip Error: Unrecognized format type ' + \
         ' for token ' + key))
      sys.exit() 

# Loop over the parsedInputs dictionary to process each key:value pair
# For each pair we format the value according to the format we defined in
# parsnip.corvus.formats. If we don't have a format defined for key, then
# we leave it as is but issue a warning.

  for (key,value_block) in list(parsedInputs.items()):

# Before we do anything, make sure this key is present in the formats. If it
# is not, then we print a warning and just keep the original format.
    if key not in list(formatsDict.keys()):
      print(('Parsnip Warning: no format found for ' + key + ', keeping as is.'))
      continue

# We also check to see if the format associated with this key is 'P', or pass.
# This means we again leave the key value unformatted.
# NOTE: We check the first element of the format list ONLY, ignoring the rest.
    if formatsDict[key][0].split('\n')[0].split()[0] == Pass_Indicat:
      continue

# Initialize the new value block
    new_value_block = []

# Split the input value block into a list of lines
    value_list = value_block.split('\n')

# Get the key-associated format in list form too
    format_list = [ l.split() for l in formatsDict[key][0].split('\n') ]

# Debug: FDV
#   print value_list
#   print format_list

# Check if any lines in in the format are empty, just in case
    if not all(format_list):
      print(('Parsnip Error: Missing format type' + \
         ' for token ' + key))
      sys.exit() 

# Check if we have RepData in the last line.
# If we do we generate a format list that has the last format line repeated
# sufficient times that all the data values are covered.
    if format_list[-1][0] == RepData_Indicat:

# Check we have the minimum number of values to put into the listed formats
      if len(value_list) < len(format_list)-1:
        print(('Parsnip Error: Not enough values for format of token ' + key))
        sys.exit() 
      
# Get the number of times we have to repeat the last format line
      nrep = len(value_list) - (len(format_list)-2)

# Redefine the format list now with the right number of repetitions
      format_list = format_list[0:-2] + [format_list[-2]]*nrep

# Debug: FDV
#     print format_list

# Check if the values and formats lists are consistent
# This is obvious for the ones we created for a repeated format, but might not
# be the case for defined length formats
    if ( len(format_list) != len(value_list) ):
      print(('Parsnip Error: Value for ' + key + ' has wrong format.'))
      sys.exit()

# After this we are ready to process the values line by line
    for (value_line,format_line) in zip(value_list,format_list):
      conv_line = formatLine(key,value_line,format_line)
      new_value_block.append(conv_line)

# Finally, we update the original value block with the newly formatted one
    parsedInputs.update({key:new_value_block})

  return parsedInputs

def formatLine(key,value_line,format_list):

# Convert split the value line, which is currently a string
  value_list = value_line.split()

# Debug: FDV
# print value_list, format_list

# Initialize the formatted value line
  conv_line = []

# Check if the current format line has repeated data
  if format_list[-1] == RepData_Indicat:

# Check we have the minimum number of values to put into the listed formats
    if len(value_list) < len(format_list)-1:
      print(('Parsnip Error: Not enough values for format of token ' + key))
      sys.exit() 
      
# Get the number of times we have to repeat the last format
    nrep = len(value_list) - (len(format_list)-2)

# Redefine the format list now with the right number of repetitions
    format_list = format_list[0:-2] + [format_list[-2]]*nrep

# Debug: FDV
#   print format_list

# Check if the current value and format lists have the same number of elements
  if len(value_list) != len(format_list):
    print(('Parsnip Error: Value for ' + key + ' has wrong format.'))
    sys.exit()

# After this we are ready to process the values element by element
  for (value,format) in zip(value_list,format_list):
    try:
      new_value = Typ_Conv_Funcs[format](value)
      conv_line.append(new_value)
    except ValueError:
      print(('Parsnip Error: Value ' + value + ' for ' + key + ' has wrong format.'))
      sys.exit()

  return conv_line

def readConfig(file):
    # Define config grammar
    config = Group(tokenGrammar + bodyGrammar + stubGrammar)
    configFile = Dict(OneOrMore(config))

    try:
        cleanConfig = commentGrammar.transformString(file.read())
        return configFile.parseString(cleanConfig)
    except ParseException as pe:
        print(("Parsnip Error: invalid input:", pe))

def readConfig_for_Doc(file):

# Redefine the comment grammar to keep the documentation lines that start
# with "%"
    commentGrammar = Word(Comm_Chars) + restOfLine
    commentGrammar.setParseAction(replaceWith(''))

# Redefine the components of the config grammar to keep the documentation
# information
    tokenGrammar = Word(alphanums + '.-_')
    docGrammar = Word('%') + restOfLine
    lbrace = Literal('{').suppress()
    rbrace = Literal('}').suppress()
    value = Word(printables.replace('{','').replace('}',''))
    keepWS = Combine(value + ZeroOrMore(White() + value), adjacent=False)
    bodyGrammar = (lbrace + Optional(keepWS, default=None) + rbrace) | value
    stubGrammar = Word(printables)

    # Define config grammar
    config = Group(tokenGrammar + ZeroOrMore(docGrammar) + bodyGrammar + stubGrammar)
    configFile = Dict(OneOrMore(config))

    try:
        cleanConfig = commentGrammar.transformString(file.read())
        return configFile.parseString(cleanConfig)
    except ParseException as pe:
        print(("Parsnip Error: invalid input:", pe))

def parseMetaData(String):

    testdata = """
Metadata-Version: 1.0
Name: corvus
Version: 0.9.0
Summary: Property-driven Scientific Workflow Manager
Home-page: http://feffproject.org
Author: S. Story, F. D. Vila, J. J. Kas, S. D. Pemmaraju
Author-email: feff@uw.edu
License: UNKNOWN
Description: UNKNOWN
Platform: UNKNOWN
"""

# Define the components of the config grammar to keep the documentation
# information
    key    = Word(alphanums + '-')
    colon  = Suppress(':')
    blanks = ZeroOrMore(' ').suppress()
    value  = restOfLine('value')
    metadataGrammar = Dict(OneOrMore(Group(key + colon + blanks + value)))

#   print testdata
#   print metadataGrammar.parseString(testdata).asDict()

    try:
        return metadataGrammar.parseString(String).asDict()
    except ParseException as pe:
        print(("Parsnip Error: invalid input in parseMetaData:\n", pe))

def readInput(file):
    # Here we do some preprocessing of the input file
    cleanInput = file.read()
    for g in [commentGrammar, leadingWhiteGrammar]:
        cleanInput = g.transformString(cleanInput)
    
    # Define input grammar
    input = Group(tokenGrammar + bodyGrammar)
    inputFile = Dict(OneOrMore(input))

    # This is used to pull out a list of just the tokens to check for duplicates
    tokenList = OneOrMore(tokenGrammar + Suppress(bodyGrammar))

    # Parse input file
    try: 
        # Check for duplicate tokens => error + exit
        tokens = tokenList.parseString(cleanInput).asList()
        duplicates = [t for t in tokens if tokens.count(t) > 1]
        if duplicates:
            output = 'Parsnip Error: duplicate input entries detected...'
            for t in set(duplicates):
                ouptut.append(str(duplicates.count(t)) + ' x ' + t)
            raise LookupError(output)
        # Return parsed dictionary            
        return inputFile.parseString(cleanInput)
    except ParseException as pe:
        print(("Parsnip Error: invalid input:", pe ))

def usage():
    print("Usage: python parsnip.py -c <configfile> -i <inputfile>")

def main(argv):
    # Parse command line arguments
    try: 
        shortopts = 'c:i:m:'
        longopts = ['config=', 'input=', 'mode=']
        opts, args = getopt.getopt(argv[1:], shortopts, longopts)
    except getopt.GetoptError as err:
        print((str(err)))
        usage()
        sys.exit(2)
    if not opts:
        usage()
        sys.exit()

    configFile = ''
    inputFile = ''
    mode = defaultMode
    for opt, arg in opts:
        if opt in ('-c', '--config'):
            configFile = arg
        elif opt in ('-i', '--input'):
            inputFile = arg
        elif opt in ('-m', '--mode'):
            mode = arg
   
    parse(configFile, inputFile, mode)

if __name__=='__main__': main(sys.argv)
