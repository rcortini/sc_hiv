#!/bin/bash

sc_hiv_home="/users/gfilion/rcortini/work/CRG/projects/sc_hiv"

for sample in $(find ../data/fastqs -name "*_1.fastq.gz"); do
  sample_id=${sample%%_1.fastq.gz}
  sample_id=${sample_id##../data/fastqs/}
  # echo $sample_id

  pbs_out="../data/fastqs/pbs_scripts/$sample_id.pbs"
  # create the script
  cat salmon_quant.pbs.in |\
    sed -e s,@SC_HIV_HOME@,$sc_hiv_home,g |\
    sed -e s,@SAMPLE_ID@,$sample_id,g |\
  tee > $pbs_out
done
