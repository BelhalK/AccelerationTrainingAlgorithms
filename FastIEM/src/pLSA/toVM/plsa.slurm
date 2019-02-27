#!/bin/sh
#SBATCH --job-name=belhal_plsa1
#SBATCH --account=small 
#SBATCH --mail-user=belhal.karim@polytechnique.edu --mail-type=ALL

#SBATCH --ntasks=4
#SBATCH --partition=interactive
python plsaforVM_EM.py
python plsaforVM_IEM.py
python plsaforVM_OEM.py
python plsaforVM_OEMVR.py

