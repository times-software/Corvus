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
# Still copy corvus.conf if one exists. Just don't do anything
# if it doesn't exist. We will search for programs on first
# run if the program doesn't find $HOME/.Corvus/corvus.conf.
if Path('corvus.conf').is_file():
    # J. Kas - now going to copy all data files to ~/.Corvus
    copy('corvus.conf',str(config_pth))
# always copy parsnip.corvus.config
copy(str(Path('./corvutils/parsnip.corvus.config')),str(util_pth))
#else:
    #print('Writing corvus.conf')
    #f = open("corvus.conf", "w")
    #print(os.getcwd())
    #f.write("[Executables]\n")
    #f.write("dmdw     :\n")
    #f.write("feff     : " + str(feff_path) + "\n") 
    #f.write("abinit   :\n")
    #f.write("nwchem   :\n")
    #f.write("orca     :\n")
    #f.write("gaussian :\n")
    #f.write("vasp     :\n")
    #f.write("siesta   :\n")
    #f.write("ocean    :\n")
    #f.write("cif2cell :\n")
    #f.write("phsf     :\n")
    #f.write(" \n")
    #f.write("[Defaults]\n")
    #f.write("prefix      : Corvus\n")
    #f.write("inputsuffix : .inp\n")
    #f.write("savesuffix  : .nest\n")
    #f.write("checkpoints : off\n")
    #f.write("parallelrun :\n")
#
    #f.close()
    #copy('corvus.conf',str(config_pth))
    #copy(str(Path('./corvutils/parsnip.corvus.config')),str(util_pth))
#
setuptools.setup(name='corvus',
      version='1.0.9',
      description='Property-driven Scientific Workflow Manager',
      author='S. Story, F. D. Vila, J. J. Kas, S. D. Pemmaraju, J. J. Rehr',
      author_email='feff@uw.edu',
      maintainer='F. D. Vila, J. J. Kas, S. D. Pemmaraju, J. J. Rehr',
      maintainer_email='feff@uw.edu',
      url='http://feffproject.org',
      #scripts=['bin/run-corvus'],
      packages=setuptools.find_packages(),
      # J Kas - Moved corvus.conf to corvus/config since pip/setuptools don't like names that start with the module name?
      package_data={'corvutils':['parsnip.corvus.config','parsnip.corvus.formats'],'corvus':['config']},
      install_requires=['lmfit','pymatgen'],
      )
