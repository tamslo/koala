#!/bin/bash

# Basically from https://github.com/khayer/aligner_benchmark/blob/master/templates/star.sh

read_length="$1"
alignment_directory="$2"
truth_file="$3"
cd $alignment_directory

grep -v "^@" ./*Aligned.out.sam | sort -t'.' -k 2n > output.sam
ruby fix_sam.rb -r $read_length output.sam > fixed.sam
ruby compare2truth_multi_mappers.rb -r $read_length $truth_file fixed.sam > comp_res_multi_mappers.txt
ruby compare2truth.rb -r $read_length $truth_file fixed.sam > comp_res.txt
