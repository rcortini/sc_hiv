#!/bin/bash

sc_hiv_rootdir=$HOME/work/CRG/projects/sc_hiv
fastq_dir=$sc_hiv_rootdir/data/fastq
out_dir=$fastq_dir/quality_control

# clear all contents of output directory
rm -rf $out_dir

# make output directory
mkdir -p $out_dir

# start with getting all the fastq names
fastq_list=$(find $fastq_dir -name "*.fq.gz")

for fastq in $fastq_list; do
  echo "Processing $fastq"
  fastqc --extract -q -o $out_dir $fastq
done
