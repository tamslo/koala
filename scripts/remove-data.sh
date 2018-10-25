#!/bin/bash

rm -rf data/experiments
rm -rf data/errored

cd data/datasets
for dataset_json in *.json; do
  dataset_directory=${dataset_json%".json"}
  for cache_directory in $dataset_directory/*; do
    directory_name=${cache_directory#"$dataset_directory/"}
    [ "$directory_name" == "pfal" ] && is_pfal=true || is_pfal=false
    [ "$directory_name" == "hg19" ] && is_hg19=true || is_hg19=false
    [ "$directory_name" == "hg38" ] && is_hg38=true || is_hg38=false
    if $is_pfal || $is_hg19 || $is_hg38; then
      rm -rf $cache_directory
    fi
  done
done

cd ../references
for file in ./*; do
  if [[ ($file == *_*) || ($file == *.fai) || ($file == .dict)  ]]
  then
    rm -r $file
  fi
done
