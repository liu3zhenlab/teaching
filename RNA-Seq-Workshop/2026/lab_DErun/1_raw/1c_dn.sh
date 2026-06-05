#!/bin/bash
# must check
meta_file=dataset.sample.txt
srr_col=2
rename_col=4

# running
perl ../utils/scripts/fasterq_dump.sbatch.pl \
  --fdpath ../utils/sratoolkit/bin/fasterq-dump \
  --in $meta_file \
  --srrcol $srr_col

# create a script to rename downloaded files
cut $meta_file -f $srr_col,$rename_col | \
  grep "^[EDS]RR" | \
  sed 's/^/rename /g' | sed 's/\t/ /g' | \
  sed 's/$/ *fastq/g' > 2c_rename.sh

