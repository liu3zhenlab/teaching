#!/bin/bash -l

module load Java/1.8.0_192

ref=1_db/Reference_CN_Wuhan_Jan012020.fasta
bamlist=2o_bamlist
ls -1 *bam | sed 's/^/-I /g' > $bamlist
vcfout=3o_gatk.vcf
txtout=3o_gatk.txt

# gatk
/homes/liu3zhen/software/GATK/gatk4/gatk-4.1.9.0/gatk HaplotypeCaller \
	-R $ref \
	--arguments_file $bamlist \
	--sample-ploidy 1 \
	--min-pruning 1 \
	-O $vcfout

# clean_up
rm $bamlist

