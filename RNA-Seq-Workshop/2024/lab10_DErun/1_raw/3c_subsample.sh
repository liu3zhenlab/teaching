#!/bin/bash
head -n 10000 norm1_1.fastq > norm1sub_1.fastq
head -n 10000 norm1_2.fastq > norm1sub_2.fastq
rm norm1_1.fastq
rm norm1_2.fastq

