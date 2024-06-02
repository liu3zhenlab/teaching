#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8g
#SBATCH --time=1-00:00:00
module load STAR
STAR --runThreadN 4 \
  --runMode genomeGenerate \
  --genomeDir "." \
  --genomeFastaFiles chr10.fasta \
  --sjdbGTFfile chr10.gtf \
  --sjdbOverhang 100

