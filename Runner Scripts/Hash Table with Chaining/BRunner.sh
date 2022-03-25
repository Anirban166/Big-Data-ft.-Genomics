#!/bin/bash
#SBATCH --job-name=hashTableConstructsB         
#SBATCH --output=/scratch/ac4743/Problem1B.txt
#SBATCH --time=00:10:00 
#SBATCH --mem=10000

srun ./hashTableConstructs Problem1B sampleGenome.fasta