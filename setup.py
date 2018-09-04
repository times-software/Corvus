from distutils.core import setup
from os import walk
import sys

# Generate the examples data by scavenging the examples directory
Examples_Dir = 'examples'

Examples_Data = []
for root, dirs, files in  walk(Examples_Dir):
  for file in files:
    Examples_Data.append((root,[root+'/'+file]))

# Debug: FDV
#print Examples_Data
#print [('doc',['doc/Corvus_Manual.txt'])]+Examples_Data
#sys.exit()

setup(name='corvus',
      version='0.9.0',
      description='Property-driven Scientific Workflow Manager',
      author='S. Story, F. D. Vila, J. J. Kas, S. D. Pemmaraju',
      author_email='feff@uw.edu',
      maintainer='S. Story, F. D. Vila, J. J. Kas, S. D. Pemmaraju',
      maintainer_email='feff@uw.edu',
      url='http://feffproject.org',
      scripts=['bin/run-corvus'],
      packages=['corvus', 'corvutils', ''],
      package_dir={'corvus':'src/corvus', 'corvutils':'src/util'},
      package_data={'':['corvus.conf'],'corvutils':['parsnip.corvus.config','parsnip.corvus.formats']},
      data_files=[('doc',['doc/Corvus_Manual.txt'])]+Examples_Data
      )
