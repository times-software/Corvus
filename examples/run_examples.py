#!/usr/bin/env python3

import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

infiles = [
    "./Feff_XANES/GeCl4/GeCl4.in",
    "./Feff_RIXS/LMnACACNBPh4/LMnACACNBPh4.inp",
    "./Doping/SrTixSn1-xO3/SrTiSnO3.in",
    "./fit/GeCl4_Fast/GeCl4.in",
    "./opcons/Diamond/Diamond.in",
    "./cfavg/CaCoO2.in",
    "./Feff_XES/GeCl4/GeCl4.in",
    "./loop/corvus.in",
]
outfiles = [
    "./Feff_XANES/GeCl4/Corvus.xanes.out",
    "./Feff_RIXS/LMnACACNBPh4/Corvus1_FEFF/rixsET-sat.dat",
    "./Doping/SrTixSn1-xO3/Corvus.cfavg_xanes.out",
    "./fit/GeCl4_Fast/fitconvergence.dat",
    "./opcons/Diamond/opconsKK.dat",
    "./cfavg/Corvus_cfavg.xanes.out",
    "./Feff_XES/GeCl4/Corvus.xes.out",
    "./loop/Corvus.loop.out",
]

assert len(infiles) == len(outfiles)
#infiles = ["./Doping/SrTixSn1-xO3/SrTiSnO3.in"]


def process_file(rel_path, expected_output, start_dir):
    print("   RUNNING: ", rel_path,flush=True)
    abs_path = os.path.abspath(
        os.path.join(start_dir, rel_path)
    )

    directory = os.path.dirname(abs_path)
    filename = os.path.basename(abs_path)
    outfile = filename + ".out"

    expected_path = os.path.join(
        directory,
        expected_output
    )
    # Remove stale output from previous runs
    if os.path.exists(expected_path):
        os.remove(expected_path)

    with open(os.path.join(directory, outfile), "w") as out:

        result = subprocess.run(
            ["run-corvus", "-i", filename],
            cwd=directory,
            stdout=out,
            stderr=out,
            text=True,
        )
    success = os.path.isfile(expected_path)
    return rel_path, success


def main(infiles,outfiles):

    start_dir = os.getcwd()

    cpu_count = os.cpu_count() or 1
    max_workers = min(len(infiles), max(1, cpu_count - 2))

    print(
        f"Running {len(infiles)} jobs using "
        f"{max_workers} worker processes."
    )

    with ProcessPoolExecutor(
        max_workers=max_workers
    ) as executor:


        futures = {}

        for infile, outfile in zip(infiles, outfiles):
            print(f"QUEUE  {infile}")

            future = executor.submit(
                process_file,
                infile,
                outfile,
                start_dir
            )

            futures[future] = infile


        for future in as_completed(futures):

            try:
                file, rc = future.result()

                if rc == 0:
                    print(f"PASS  {file}")
                else:
                    print(
                        f"FAIL  {file} "
                        f"(return code {rc})"
                    )

            except Exception as e:
                print(
                    f"CRASH {futures[future]}: {e}"
                )


if __name__ == "__main__":
    main(infiles,outfiles)
