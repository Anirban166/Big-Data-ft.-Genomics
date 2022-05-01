#!/bin/bash
#SBATCH --job-name=suffixTreeB         
#SBATCH --output=/scratch/ac4743/suffixTreeB.txt
#SBATCH --time=00:00:30 
#SBATCH --mem=10000

srun ./suffixTreeConstructs Problem1B SARSCoV2.txt