#!/bin/bash
# cd to db directory

mkdir bwa
cd bwa
ln -s ../MG1655.fasta .
bwa index MG1655.fasta

