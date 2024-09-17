#!/bin/bash

bash check_latex.sh

cd plateaus
#python record_meson_F.py
#python analysis_F.py

#python record_meson_AS.py
#python analysis_AS.py

#python record_GEVP_F.py
python analysis_GEVP_F.py

#python record_GEVP_AS.py
python analysis_GEVP_AS.py

#python record_GEVP_mixF.py
python analysis_GEVP_MIX_F.py

#python record_GEVP_mixAS.py
python analysis_GEVP_MIX_AS.py

#python plot_GEVP_F.py
#python plot_GEVP_AS.py
