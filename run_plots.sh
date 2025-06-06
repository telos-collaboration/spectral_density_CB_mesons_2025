#!/bin/bash

set -eu

bash check_latex.sh

#cd just_plotting/code
#cd topologies
#bash main.sh


cd just_plotting/code
cd final_spectrum
python spectrum_MN.py
python spectrum_ensembles_showing.py
cd ../final_matrixel
python matrix_MN.py
cd ../stability_plot
python stability.py
cd ../sfs
python sfs.py
cd ../../../lsd_out
python elaborate.py
cd ../CB_autocorrelation_decay_constant
conda activate wall_decay_constant
bash main.sh
cp -r assets/wall_comparison.pdf ../plots
cp -r assets/local_smeared_decay_constants.tex ../tables
cp -r assets/ensembles.tex ../tables
