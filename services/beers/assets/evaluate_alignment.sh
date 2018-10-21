#!/bin/bash

# Basically from https://github.com/khayer/aligner_benchmark/blob/master/templates/star.sh

alignment_tool_directory="/opt/beers_evaluation"

read_length="$1"
alignment_directory="$2"
truth_file="$3" # path needs to be absolute
cd $alignment_directory

grep -v "^@" ./Out.sam | sort -t'.' -k 2n > Sorted.sam
ruby $alignment_tool_directory/fix_sam.rb -r $read_length Sorted.sam > Fixed.sam
ruby $alignment_tool_directory/compare2truth_multi_mappers.rb -r $read_length $truth_file Fixed.sam > Evaluation.multi.txt
ruby $alignment_tool_directory/compare2truth.rb -r $read_length $truth_file Fixed.sam > Evaluation.txt
