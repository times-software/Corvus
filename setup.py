import setuptools 
#from os import walk
import sys

# Generate the examples data by scavenging the examples directory
#Examples_Dir = 'examples'

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
      package_data={'':['corvus.conf'],'corvutils':['parsnip.corvus.config','parsnip.corvus.formats']},
      data_files=[('.',['corvus.conf']),('doc',['doc/Corvus_Manual.txt'])], #+Examples_Data
      install_requires=['numpy','lmfit','cif2cell'],
      python_requires='>=2.7,<3'
      )
