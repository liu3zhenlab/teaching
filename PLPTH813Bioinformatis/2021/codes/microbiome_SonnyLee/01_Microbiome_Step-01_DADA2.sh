#!/bin/bash -l
#SBATCH --job-name=Micro

#SBATCH --mem-per-cpu=1G   
#SBATCH --time=0-01:00:00   
#SBATCH --ntasks=10
#SBATCH --nodes=1

#SBATCH --mail-user=<INSERT YOUR EMAIL ADDRESS>
#SBATCH --mail-type=ALL   

module load QIIME2/2019.7

#Import data into Qiime 2
qiime tools import \
  --type EMPSingleEndSequences \
  --input-path emp-single-end-sequences \
  --output-path emp-single-end-sequences.qza

#Demultiplexing sequences
qiime demux emp-single \
  --i-seqs emp-single-end-sequences.qza \
  --m-barcodes-file sample-metadata.tsv \
  --m-barcodes-column barcode-sequence \
  --o-per-sample-sequences demux.qza \
  --o-error-correction-details demux-details.qza

qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux.qzv

#DADA2
qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left 0 \
  --p-trunc-len 120 \
  --o-representative-sequences rep-seqs-dada2.qza \
  --o-table table-dada2.qza \
  --o-denoising-stats stats-dada2.qza
