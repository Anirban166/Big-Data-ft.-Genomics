#!/bin/bash
#SBATCH --job-name=BLASTAlignment        
#SBATCH --output=/scratch/ac4743/BLASTAlignment.txt
#SBATCH --time=12:00:00 
#SBATCH --mem=10000

srun ./alignmentConstructs Problem1A sampleDataset.fa sampleGenome.fasta