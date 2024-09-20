#!/bin/bash

bash check_latex.sh

cd plateaus

python analysis_GEVP_F.py
python analysis_GEVP_AS.py
python analysis_GEVP_MIX_F.py
python analysis_GEVP_MIX_AS.py

