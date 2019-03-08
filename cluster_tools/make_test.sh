#!/bin/bash

sc_hiv_home="/users/gfilion/rcortini/work/CRG/projects/sc_hiv"
sample_id="CD2ALANXX_6_N714-S511-NX-xt"

cat salmon_quant.pbs.in |\
  sed -e s,@SC_HIV_HOME@,$sc_hiv_home,g |\
  sed -e s,@SAMPLE_ID@,$sample_id,g |\
tee > salmon_quant.pbs
