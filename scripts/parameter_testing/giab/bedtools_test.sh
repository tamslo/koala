#!/bin/bash

# Script to be run in GIAB container to create different BED intersections
# and run evaluations on them
# To run in background and write logfile, run ./bedtools_test.sh > logfile 2>&1 &
# Required tree:
# giab-playground/
# ├── annotations
# │   ├── hg38_coding_exons.bed
# │   ├── hg38_editing_sites.bed
# │   ├── hg38_exons.bed
# │   └── hg38_exons_overhang_10.bed
# ├── base_coverage_to_bed.py
# ├── bedtools_test.sh
# ├── giab
# │   ├── confidence_calls.bed
# │   └── confidence_calls.vcf
# ├── hg38.fa
# ├── hg38.fa.fai
# ├── Star.bam (alignment the BED files are based on)
# ├── Star.bam.bed
# ├── Star.bam.bed.merged
# ├── Star.bam.coverage
# ├── Star.bam.coverage.10-filtered.bed
# ├── Star.bam.coverage.2-filtered.bed
# ├── StarOpossumChr3.vcf (output of pipeline)
# ├── StarOpossumChr4.vcf (output of pipeline)
# └── StarOpossumChr5.vcf (output of pipeline)

giab_confidence_regions_path=giab/confidence_calls.bed
editing_sites_path=annotations/hg38_editing_sites.bed
declare -a chromosomes=("3" "4" "5")
tests_directory=tests
mkdir -p $tests_directory

function evaluate() {
  local confidence_calls_path=giab/confidence_calls.vcf
  local confidence_regions_path=$1
  local test_directory=$2
  for chromosome in "${chromosomes[@]}"
  do
    local query_calls_path=StarOpossumChr$chromosome.vcf
    /opt/hap.py/bin/hap.py $confidence_calls_path $query_calls_path \
        -f $confidence_regions_path -o $test_directory/Evaluation_chr$chromosome \
        -r hg38.fa --location chr$chromosome
  done
}

function intersect() {
  local intersection_regions_path=$1
  local confidence_regions_path=$2
  local sorted=$3
  bedtools intersect -a $giab_confidence_regions_path \
    -b $intersection_regions_path $sorted > $confidence_regions_path
}

function subtract() {
  local subtraction_regions_path=$1
  local subtracted_regions_path=$2
  bedtools subtract -a $subtraction_regions_path \
    -b $editing_sites_path > $subtracted_regions_path
}

function run_no_editing_sites_task() {
  echo -e "\n$1 without RNA editing sites\n"
  no_editing_sites_directory="$2"_no_editing_sites
  mkdir $no_editing_sites_directory
  cleaned_confidence_regions_path=$no_editing_sites_directory/confidence_regions.bed
  subtract $3 $cleaned_confidence_regions_path
  evaluate $cleaned_confidence_regions_path $no_editing_sites_directory
}

function run_intersection_task() {
  echo -e "\n$1\n"
  test_directory=$tests_directory/$2
  if [ ! -d $test_directory ]; then
    intersection_regions_path=$3
      if [ -f $intersection_regions_path ]; then
        mkdir $test_directory
        confidence_regions_path=$test_directory/confidence_regions.bed
        intersect $intersection_regions_path $confidence_regions_path $4
        evaluate $confidence_regions_path $test_directory
        run_no_editing_sites_task "$1" $test_directory $confidence_regions_path
      else
        echo -e "$intersection_regions_path is not present"
      fi
  else
    echo -e "$test_directory already present"
  fi
}

task_name="No intersection"
echo -e "\n$task_name\n"
test_directory=$tests_directory/no_intersection
if [ ! -d $test_directory ]; then
  mkdir $test_directory
  evaluate $giab_confidence_regions_path $test_directory
  run_no_editing_sites_task "$task_name" $test_directory $giab_confidence_regions_path
else
  echo -e "$test_directory already present"
fi

run_intersection_task "Intersection with coding exons" \
  coding_exons annotations/hg38_coding_exons.bed

run_intersection_task "Intersection with exons" \
  exons annotations/hg38_exons.bed

run_intersection_task "Intersection with exons overhang" \
  exons_overhang annotations/hg38_exons_overhang_10.bed

run_intersection_task "Intersection with BAM regions (from BED)" \
  bam_bed Star.bam.bed "-sorted"

run_intersection_task "Intersection with BAM regions (from merged BED)" \
  bam_bed_merged Star.bam.bed.merged

run_intersection_task "Intersection with BAM regions with coverage (2) (from merged BED)" \
  bam_bed_merged_coverage_2 Star.bam.coverage.2-filtered.bed

run_intersection_task "Intersection with BAM regions with coverage (10) (from merged BED)" \
  bam_bed_merged_coverage_10 Star.bam.coverage.10-filtered.bed


echo -e "\nConcatenating evaluation summaries\n"
concatenated_file="Evaluation.summary.all.csv"
first_file=true
for evaluation_file in $tests_directory/*/Evaluation_*.summary.csv; do
    if [ $first_file = true ]; then
      head -1 $evaluation_file | sed -e 's/^/File,/' > $concatenated_file
      first_file=false
    fi

    tail --lines=+2 $evaluation_file | awk -v file="$evaluation_file", '{print file $0}' >> $concatenated_file
done
