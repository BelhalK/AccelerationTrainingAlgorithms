#!/bin/sh
#SBATCH --job-name=belhal_plsa1
#SBATCH --account=small 
#SBATCH --mail-user=belhal.karim@polytechnique.edu --mail-type=ALL

#SBATCH --ntasks=1
#SBATCH --partition=interactive
python plsaforVM_OEM.py -d dataset1k.txt -s 20 -e 10