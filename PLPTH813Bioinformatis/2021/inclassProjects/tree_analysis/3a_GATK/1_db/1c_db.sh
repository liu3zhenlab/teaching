#!/bin/bash

# bwa index
bwa index Reference_CN_Wuhan_Jan012020.fasta 

# GATK index
module load SAMtools 
samtools faidx Reference_CN_Wuhan_Jan012020.fasta
samtools dict Reference_CN_Wuhan_Jan012020.fasta > Reference_CN_Wuhan_Jan012020.dict

