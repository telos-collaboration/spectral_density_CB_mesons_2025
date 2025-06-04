#!/bin/bash

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
