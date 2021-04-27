#!/bin/bash
#SBATCH --mem-per-cpu=1G
#SBATCH --time=10-00:00:00
#SBATCH --cpus-per-task=16

module load BWA
module load SAMtools
ncpu=$SLURM_CPUS_PER_TASK
ref=1_db/Reference_CN_Wuhan_Jan012020.fasta
indir=../1_selectStrains/

for qry in $indir/[MR]*fasta; do
	echo $qry
	out=`echo $qry | sed 's/.*\///'g | sed 's/.fasta//g'`
	group="@RG\tID:${out}\tSM:${out}"
	echo $group
	perl ../scripts/fa2fq.pl --fa ${qry} --qual 40 > ${out}.fastq
	
	bwa mem -R $group $ref ${out}.fastq > ${out}.sam
	# bam and sort
	samtools view -b -@ ${ncpu} ${out}.sam | samtools sort -@ ${ncpu} -o ${out}.bam
	samtools index ${out}.bam
	rm ${out}.fastq
	rm ${out}.sam
done

