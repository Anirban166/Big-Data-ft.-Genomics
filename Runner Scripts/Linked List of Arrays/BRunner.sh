#!/bin/bash
#SBATCH --job-name=queryReadsetAgainstGenome       
#SBATCH --output=/scratch/ac4743/Problem1B.txt
#SBATCH --time=69:00:00 
#SBATCH --mem=10000

srun ./linkedListConstructs Problem1B sampleDataset.fa sampleGenome.fasta