#!/bin/bash

module load BWA

### specify input files
ref=../db/bwa/MG1655.fasta
pe1=/homes/liu3zhen/teaching/datasets/MG1655_illumina/MG1655.pair1.fq
pe2=/homes/liu3zhen/teaching/datasets/MG1655_illumina/MG1655.pair2.fq

### alignment
bwa mem -T 30 $ref $pe1 $pe2 1>aln.sam 2>aln.log

