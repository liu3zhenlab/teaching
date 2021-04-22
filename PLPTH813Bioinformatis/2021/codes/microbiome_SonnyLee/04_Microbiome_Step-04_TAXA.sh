#!/bin/bash -l
#SBATCH --job-name=Micro

#SBATCH --mem-per-cpu=1G   
#SBATCH --time=0-01:00:00   
#SBATCH --ntasks=10
#SBATCH --nodes=1

#SBATCH --mail-user=<INSERT YOUR EMAIL ADDRESS>
#SBATCH --mail-type=ALL   

module load QIIME2/2019.7

qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization taxa-bar-plots.qzv
