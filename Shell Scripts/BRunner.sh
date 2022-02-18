#!/bin/bash
#SBATCH --job-name=bruteforceSearchForDSOneInTwo         
#SBATCH --output=/scratch/ac4743/Problem1B.txt
#SBATCH --time=96:00:00 
#SBATCH --mem=5000

srun ./matrix Problem1B sampleDataset.fa