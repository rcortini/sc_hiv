#!/bin/bash

# check for proper invocation
if [ $# -ne 3 ]; then
  echo >&2 "Usage: fastq_process <fastq_dir> <salmon_index> <gtf>"
  exit 1
fi

# get parameters from command line
fastq_dir=$1
salmon_index=$2
gtf=$3

# iterate over all the fastq files
for fq_1 in $(find $fastq_dir -maxdepth 1 -name "*_1.fastq.gz"); do
  name=${fq_1%%_1.fastq.gz}
  name=${name##$fastq_dir/}
  fq_2="$fastq_dir/$name"_2.fastq.gz
  out_dir=$fastq_dir/postprocess/$name
  salmon quant -i $salmon_index -l A -1 $fq_1 -2 $fq_2 -p 8 -g $gtf --seqBias \
  --gcBias --posBias --output $out_dir --validateMappings
done
