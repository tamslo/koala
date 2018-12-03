#!/bin/bash

happy_directory="/opt/hap.py/bin"

path_prefix="$1"  # path needs to be absolute
vcf_file_name="$2"
reference_id="$3"
additional_commands="$4"
additional_commands="${additional_commands//%%/ }"
cd $happy_directory

HGREF=/data/references/$reference_id.fa \
  ./hap.py /giab/$reference_id/confidence_calls.vcf \
  $path_prefix/$vcf_file_name \
  -f /giab/$reference_id/confidence_calls.bed \
  -o $path_prefix/Evaluation \
  $additional_commands
