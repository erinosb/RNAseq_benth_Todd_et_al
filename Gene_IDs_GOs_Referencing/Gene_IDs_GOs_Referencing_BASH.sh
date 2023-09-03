#!/bin/bash
#SBATCH --nodes=1                       # Number of requested nodes
#SBATCH --ntasks=4                      # Number of requested cores
#SBATCH --time=0:10:00                  # Max walltime
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --constraint=ib
#SBATCH --output=./output/python_%j.out         # Output file name

# purge all existing modules
module purge

# Load the python module
module load python
module load intel impi
module load anaconda
conda activate mycustomenv

# Run Python Script
mpirun -np 4 python Gene_IDs_GOs_Referencing.py
