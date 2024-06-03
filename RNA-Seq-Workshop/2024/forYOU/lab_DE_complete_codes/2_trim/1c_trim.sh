#!/bin/bash
perl ../utils/scripts/trimmomatic.sbatch.pl \
--mem 4G \
--time 1-00:00:00 \
--trim_shell "../utils/scripts/trimmomatic.pe.sh" \
--trimmomatic "../utils/Trimmomatic-0.38/trimmomatic-0.38.jar" \
--adaptor_file "../utils/Trimmomatic-0.38/adapters/TruSeq3-PE.fa" \
--indir "../1_raw" \
--outdir "." \
--fq1feature "_1.fastq" \
--fq2feature "_2.fastq" \
--threads 4 \
--min_len 40

