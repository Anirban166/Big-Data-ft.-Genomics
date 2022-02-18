#!/bin/bash
#SBATCH --job-name=sortDSTwoAndBinarySearchForDSOneInTwo    
#SBATCH --output=/scratch/ac4743/Problem1C.txt
#SBATCH --time=96:00:00 
#SBATCH --mem=5000

srun ./arrayConstructs Problem1C sampleDataset.fa