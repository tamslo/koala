#!/bin/bash

# Basically from https://github.com/khayer/aligner_benchmark/blob/master/templates/star.sh

happy_directory="/opt/hap.py/bin"

path_prefix="$1"  # path needs to be absolute
vcf_file_name="$2"
reference_id="$3"
filter_postfix="$4"
cd $happy_directory

HG19=/data/references/$reference_id.fa \
  ./hap.py /giab/$reference_id/confidence_calls"$filter_postfix".vcf \
  $path_prefix/$vcf_file_name \
  -f /giab/$reference_id/confidence_calls"$filter_postfix".bed \
  -o $path_prefix/Evaluation
