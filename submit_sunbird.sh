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

conda activate my-new-env2

bash reproduce_everything.sh

#cd plateaus

#bash run_plateaus.sh

#conda activate analysis-env4

#cd lsd_out


#bash run_spectral_densities.sh
#bash reproduce_everything.sh
#python3 analyse_data.py
#python3 print_samples.py
