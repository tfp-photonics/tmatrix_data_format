#!/bin/sh
#SBATCH --account=hk-project-bit2blobs
#SBATCH --job-name=param_sweep
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=24:00:00
#SBATCH --partition=cpuonly
#SBATCH --output=%j_%x.log

CONDAENV=$TMP/conda

source $(conda info --base)/etc/profile.d/conda.sh
srun --exclusive -n 1 -c 16 --mem=125400mb \
    conda env create -f environment.yml --prefix $CONDAENV

conda activate $CONDAENV

mpirun -n 64 --bind-to core --map-by core python main.py
