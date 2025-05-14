import subprocess, sys
from pathlib import Path

def check_abort(p=None,cf=None):
    return False
    print('Checking abort signal.', cf)
    path = Path('abort_corvus.txt')
    if path.exists():
        print('File exists.')
        f = open('abort_corvus.txt','r')
        line = f.readlines()[0].strip()
        f.close()
        if line == 'True':
            print('Aborting ...')
            if p is None:
                print('no process. Exiting system.')
                sys.exit()
                return True
            else:
                print('Killing process.')
                p.kill()
                sys.exit()
                return True
        else:
            print('Not aborting.')
            return False
    else:
        print('Not aborting. Abort file not found.')
        pass
        return False
        
