#!/bin/bash
#SBATCH --job-name=prefixTreeA         
#SBATCH --output=/scratch/ac4743/Problem1A.txt
#SBATCH --time=00:01:45 
#SBATCH --mem=10000

srun ./prefixTreeConstructs Problem1A SARSCoV2.txt