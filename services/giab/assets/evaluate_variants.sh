#!/bin/bash

# Basically from https://github.com/khayer/aligner_benchmark/blob/master/templates/star.sh

happy_directory="/opt/hap.py/bin"

path_prefix="$1"  # path needs to be absolute
vcf_file_name="$2"
reference_id="$3"
cd $happy_directory

./hap.py /giab/$reference_id/confidence_calls.vcf \
  $path_prefix/$vcf_file_name \
  -f /giab/$reference_id/confidence_calls.bed \
  -o $path_prefix/Evaluation \
  -r /data/references/$reference_id.fa
