#!/bin/tcsh
#SBATCH --job-name=Crv_Tst
#SBATCH --partition=debug
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --error=%j.serr
#SBATCH --output=%j.sout
##SBATCH --signal=B:INT@180
##SBATCH -C knl,quad,cache
#SBATCH -C haswell

module load python nwchem

which python

setenv PYTHONPATH /global/u2/f/fvila/bin/Corvus-Merge/lib/python2.7/site-packages

/global/u2/f/fvila/bin/Corvus-Merge/bin/run-corvus -i GeCl4.in -c --parallelrun "srun --exclusive -n 4"

