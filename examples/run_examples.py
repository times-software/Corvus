import re
import sys, os
from contextlib import contextmanager
from corvus.controls import oneshot
from concurrent.futures import ThreadPoolExecutor, as_completed
# Run corvus with --version to load the libraries. This will make the first run load faster.
infiles=[
"./Feff_XANES/KSpace/Graphite.in",
"./Feff_XANES/GeCl4/GeCl4.in",
"./Feff_RIXS/LMnACACNBPh4/LMnACACNBPh4.inp",
"./Doping/SrTixSn1-xO3/SrTiSnO3.in",
"./fit/GeCl4_Full/GeCl4.in",
"./fit/GeCl4_Fast/GeCl4.in",
"./opcons/Diamond/Diamond.in",
"./opcons/Corundum/Corundum.in",
"./cfavg/CaCoO2.in",
"./Feff_XES/GeCl4/GeCl4.in",
"./loop/corvus.in"
]

print('\n\n#####################################################')
print('   Running examples. This will take time.')
print('#####################################################\n\n')
@contextmanager
def redirect_fd(filename):
    with open(filename, "w") as f:
        old_stdout = os.dup(1)
        old_stderr = os.dup(2)
        os.dup2(f.fileno(), 1)
        os.dup2(f.fileno(), 2)
    try:
        yield
    finally:
        os.dup2(old_stdout, 1)
        os.dup2(old_stderr, 2)
        os.close(old_stdout)
        os.close(old_stderr)




def main(files,max_threads):
    start_dir = os.getcwd()

    with ThreadPoolExecutor(max_workers=max_threads) as executor:
        futures = [
            executor.submit(process_file, f, start_dir)
            for f in files
        ]

        for future in as_completed(futures):
            try:
                completed_file = future.result()
                print(f"Completed: {completed_file}")
            except Exception as e:
                print(f"Failed: {e}")

def process_file(file,start_dir):
    os.chdir(start_dir)
    directory = os.path.dirname(os.path.abspath(file))
    print(directory)
    filename = os.path.basename(file)
    print("    File: ", file)
    try:
        os.chdir(directory)
        with redirect_fd(filename + '.out'):
            sys.argv = ['run-corvus','-i',filename]
            sys.argv[0] = re.sub(r'(-script\.pyw|\.exe)?$', '', sys.argv[0])
            try:
                oneshot()
            except SystemExit:
                print('Corvus failed for file: ', file)

    except Exception as e:
        print(f"Error processing {file}: {e}")

if __name__ == "__main__":
    # Same default logic used by ThreadPoolExecutor
    max_threads = min(len(infiles), os.cpu_count() - 4)
    print('Running with ', max_threads, ' processes.')
    start_dir = os.getcwd()
    main(infiles,max_threads)
