#!/bin/bash -l
#SBATCH --job-name=Micro

#SBATCH --mem-per-cpu=1G   
#SBATCH --time=0-01:00:00   
#SBATCH --ntasks=10
#SBATCH --nodes=1

#SBATCH --mail-user=<INSERT YOUR EMAIL ADDRESS>
#SBATCH --mail-type=ALL   


wget \
  -O "sample-metadata.tsv" \
  "https://data.qiime2.org/2021.2/tutorials/moving-pictures/sample_metadata.tsv"

mkdir emp-single-end-sequences

wget \
  -O "emp-single-end-sequences/barcodes.fastq.gz" \
  "https://data.qiime2.org/2021.2/tutorials/moving-pictures/emp-single-end-sequences/barcodes.fastq.gz"
  
  wget \
  -O "emp-single-end-sequences/sequences.fastq.gz" \
  "https://data.qiime2.org/2021.2/tutorials/moving-pictures/emp-single-end-sequences/sequences.fastq.gz"