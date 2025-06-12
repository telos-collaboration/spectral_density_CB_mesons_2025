#!/bin/bash

#SBATCH --time=1-23:59:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --job-name chimera_lsd
#SBATCH --account=scw1019
#SBATCH --partition=compute
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err


module load anaconda/2024.06
source activate

conda activate snakemake_x86
snakemake --cores --use-conda
