#!/bin/bash
#SBATCH --job-name=prefixTreeB         
#SBATCH --output=/scratch/ac4743/Problem1B.txt
#SBATCH --time=00:01:45 
#SBATCH --mem=10000

srun ./prefixTreeConstructs Problem1B SARSCoV2.txt