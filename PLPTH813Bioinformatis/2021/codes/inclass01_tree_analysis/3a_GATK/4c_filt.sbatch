#!/bin/bash -l

module load Java/1.8.0_192

ref=1_db/Reference_CN_Wuhan_Jan012020.fasta
invcf=3o_gatk.vcf
outvcf=4o_gatk.filt.vcf
table=4o_gatk.filt.txt

# gatk
/homes/liu3zhen/software/GATK/gatk4/gatk-4.1.9.0/gatk SelectVariants \
	-V $invcf \
	--restrict-alleles-to BIALLELIC \
	--select-type SNP \
	-O $outvcf

# convert to table format
perl ../scripts/vcfbox.pl recode $outvcf > $table

#/homes/liu3zhen/software/GATK/gatk4/gatk-4.1.9.0/gatk VariantsToTable \
#	-R $ref \
#	-V $outvcf \
#	-F CHROM -F POS -F REF -F ALT \
#	--genotype-fields "GT" \
#	-O $table

