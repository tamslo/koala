#!/bin/bash

happy_directory="/opt/hap.py/bin"

path_prefix="$1"  # path needs to be absolute
vcf_file_name="$2"
reference_id="$3"
filter_postfix="$4"
cd $happy_directory

confidence_prefix="confidence_calls$filter_postfix"

HG19=/data/references/$reference_id.fa \
  ./hap.py /giab/$reference_id/$confidence_prefix.vcf \
  $path_prefix/$vcf_file_name \
  -f /giab/$reference_id/$confidence_prefix.bed \
  -o $path_prefix/Evaluation
