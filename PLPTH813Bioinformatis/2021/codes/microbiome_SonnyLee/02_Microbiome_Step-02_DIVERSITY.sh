#!/bin/bash -l
#SBATCH --job-name=Micro

#SBATCH --mem-per-cpu=1G   
#SBATCH --time=0-01:00:00   
#SBATCH --ntasks=10
#SBATCH --nodes=1

#SBATCH --mail-user=<INSERT YOUR EMAIL ADDRESS>
#SBATCH --mail-type=ALL   

module load QIIME2/2019.7

qiime metadata tabulate \
  --m-input-file stats-dada2.qza \
  --o-visualization stats-dada2.qzv
  
mv rep-seqs-dada2.qza rep-seqs.qza
mv table-dada2.qza table.qza

#FeatureTable and FeatureData summaries
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file sample-metadata.tsv
  
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv
  
#Generate a tree for phylogenetic diversity analyses
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

#Alpha and beta diversity analysis
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 1103 \
  --m-metadata-file sample-metadata.tsv \
  --output-dir core-metrics-results
