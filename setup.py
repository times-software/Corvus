from setuptools import setup, find_packages
import os
import re
#from os import walk
from shutil import copy
import sys
from pathlib import Path

def get_version_info_from_metadata(default_version="0.0.0"):
    """
    Reads the archived metadata file to extract the last active tag.
    Falls back to a default version if the file doesn't exist.
    """
    # Adjust this if your primary branch is 'main' instead of 'master'
    meta_branch = "master_shit"
    setup_dir = os.path.dirname(os.path.abspath(__file__))
    meta_path = os.path.join(setup_dir, "archived", "meta-master.txt")
    meta_version=default_version
    meta_commit_hash="fuck me"
    
    if not os.path.exists(meta_path):
        return default_version

    try: 
        with open(meta_path, "r") as f:
            for line in f:
                if line.startswith("Last Active Tag:"):
                    # Extract the tag name (e.g., 'v1.2.3' or '1.2.3')
                    tag = line.split(":", 1)[1].strip()
                    
                    if tag == "no-tag-found" or not tag:
                        meta_version = default_version
                    else:
                        meta_version = re.sub(r'^v', '', tag)
                    
                    # Clean the tag for PEP 440 compliance (e.g., stripping a leading 'v')
                    # Setuptools prefers strict semantic versioning formats
                elif line.startswith("Commit Hash:"):
                    meta_commit_hash = line.strip()
                    # Checks if the line contains nothing but "Commit Hash:", spaces, or newlines
                    if meta_commit_hash.replace("Commit Hash:", "").strip() == "":
                        meta_commit_hash = "commit: No hash found"
                        
                elif line.startswith("Branch:"):
                    meta_branch = line.strip()
                    # Checks if the line contains nothing but "Branch:", spaces, or newlines
                    if meta_branch.replace("Branch:", "").strip() == "":
                        meta_branch = "No branch found"

#                elif line.startswith("Commit Hash:"):
#                    commit_hash = line
#                    if not commit_hash:
#                        commit_hash = "commit: No hash found"
#                elif line.startswith("Branch:"):
#                    branch = line
#                    if not branch:
#                        branch = "No branch found"

    except Exception:
        pass
        
    return (meta_version, meta_branch, meta_commit_hash)
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
version_info = get_version_info_from_metadata()
setup(name='corvus',
      version=version_info[0],
      python_requires=">=3.12, <3.14",
      description='Property-driven Scientific Workflow Manager. ' + version_info[1] + ', ' + version_info[2],
      author='S. Story, F. D. Vila, J. J. Kas, S. D. Pemmaraju, J. J. Rehr',
      author_email='feff@uw.edu',
      maintainer='F. D. Vila, J. J. Kas, S. D. Pemmaraju, J. J. Rehr',
      maintainer_email='feff@uw.edu',
      url='http://feffproject.org',
      #scripts=['bin/run-corvus'],
      packages=find_packages(),
      # J Kas - Moved corvus.conf to corvus/config since pip/setuptools don't like names that start with the module name?
      package_data={'corvutils':['parsnip.corvus.config'],'corvus':['config']},
      #install_requires=['more_itertools','h5py==3.15.1','lmfit','mp_api','pymatgen','orjson']
      install_requires=['more_itertools','h5py','lmfit','mp_api','pymatgen','orjson']
      )
