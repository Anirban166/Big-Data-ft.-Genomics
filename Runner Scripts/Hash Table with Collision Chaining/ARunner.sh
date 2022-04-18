#!/bin/bash
#SBATCH --job-name=hashTableConstructsA         
#SBATCH --output=/scratch/ac4743/Problem1A.txt
#SBATCH --time=00:10:00 
#SBATCH --mem=10000

srun ./hashTableConstructs Problem1A sampleGenome.fasta