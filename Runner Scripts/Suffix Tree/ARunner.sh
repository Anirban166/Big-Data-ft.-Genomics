#!/bin/bash
#SBATCH --job-name=suffixTreeA         
#SBATCH --output=/scratch/ac4743/suffixTreeA.txt
#SBATCH --time=00:02:45 
#SBATCH --mem=10000

srun ./suffixTreeConstructs Problem1A SARSCoV2.txt