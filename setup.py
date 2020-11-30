import setuptools 
#from os import walk
from shutil import copy
import sys
from pathlib import Path

# Generate the examples data by scavenging the examples directory
#Examples_Dir = 'examples'

# J. Kas - copy corvus.conf to corvus/config for backward compatibility.
# Make directory ~/.Corvus if it doesn't exist
config_pth = Path.home() / ".Corvus"
config_pth.mkdir(exist_ok=True)
util_pth = config_pth / "corvutils" 
util_pth.mkdir(exist_ok=True)
try:
    # J. Kas - now going to copy all data files to ~/.Corvus
    copy('corvus.conf',str(config_pth))
    copy(str(Path('./corvutils/parsnip.corvus.config')),str(util_pth))
except:
    print('')
    print('')
    print('')
    print('')
    print('')
    print('')
    print('############################################################################')
    print('############################################################################')
    print('#        Corvus has not been configured. Please make corvus.conf file.     #')
    print('#        You can start from corvus.conf.template for reference.            #')
    print('############################################################################')
    print('############################################################################')
    sys.exit('')
# J. Kas - Don't copy examples. We don't need them in the python path.
#Examples_Data = []
#for root, dirs, files in walk(Examples_Dir):
#  for file in files:
#    Examples_Data.append((root,[root+'/'+file]))

# Debug: FDV
#print Examples_Data
#print [('doc',['doc/Corvus_Manual.txt'])]+Examples_Data
#sys.exit()

setuptools.setup(name='corvus',
      version='1.0',
      description='Property-driven Scientific Workflow Manager',
      author='S. Story, F. D. Vila, J. J. Kas, S. D. Pemmaraju, J. J. Rehr',
      author_email='feff@uw.edu',
      maintainer='F. D. Vila, J. J. Kas, S. D. Pemmaraju, J. J. Rehr',
      maintainer_email='feff@uw.edu',
      url='http://feffproject.org',
      scripts=['bin/run-corvus'],
      packages=setuptools.find_packages(),
      # J Kas - Moved corvus.conf to corvus/config since pip/setuptools don't like names that start with the module name?
      package_data={'corvutils':['parsnip.corvus.config','parsnip.corvus.formats'],'corvus':['config']},
      install_requires=['cif2cell','lmfit'],
      )
