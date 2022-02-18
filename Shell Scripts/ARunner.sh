#!/bin/bash
#SBATCH --job-name=extractSequenceFragmentCounts             
#SBATCH --output=/scratch/ac4743/Problem1A.txt
#SBATCH --time=00:30:00 
#SBATCH --mem=1000

srun ./matrix Problem1A sampleDataset.fa