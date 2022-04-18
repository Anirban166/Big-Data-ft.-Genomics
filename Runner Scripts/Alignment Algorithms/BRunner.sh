#!/bin/bash
#SBATCH --job-name=BLASTAlignmentPerfectHits         
#SBATCH --output=/scratch/ac4743/BLASTAlignmentPerfectHits.txt
#SBATCH --time=20:00:00 
#SBATCH --mem=10000

srun ./alignmentConstructs Problem1B sampleDataset.fa sampleGenome.fasta