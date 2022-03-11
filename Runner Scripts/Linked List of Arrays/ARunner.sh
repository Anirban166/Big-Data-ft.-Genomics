#!/bin/bash
#SBATCH --job-name=fragmentSearch         
#SBATCH --output=/scratch/ac4743/Problem1A.txt
#SBATCH --time=00:30:00 
#SBATCH --mem=10000

srun ./linkedListConstructs Problem1A sampleDataset.fa