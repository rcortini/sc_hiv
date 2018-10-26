#!/bin/bash

fastq_dir="$HOME/work/CRG/projects/sc_hiv/data/fastq"
out_dir=$fastq_dir/stats
mkdir -p $out_dir
log=$out_dir/read_numbers.dat
rm -rf $log

for fastq in $(find $fastq_dir -name "*1.fq.gz"); do
  sample_name=$(echo ${fastq%%_1.fq.gz} | sed -e s,$fastq_dir/,,)
  nlines=$(zcat $fastq | wc -l)
  nreads=$(echo "$nlines/4" | bc)
  echo $sample_name $nreads >> $log
done
