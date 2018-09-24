#!/bin/bash

rm -r data/experiments

cd data/datasets
for dataset_json in *.json; do
  dataset_directory=${dataset_json%".json"}
  if [[ $dataset_directory =~ ^prepared.* ]]; then
    for cache_directory in $dataset_directory/*; do
      directory_name=${cache_directory#"$dataset_directory/"}
      [ "$directory_name" == "pfal" ] && is_pfal=true || is_pfal=false
      [ "$directory_name" == "hg19" ] && is_hg19=true || is_hg19=false
      [ "$directory_name" == "hg38" ] && is_hg38=true || is_hg38=false
      if $is_pfal || $is_hg19 || $is_hg38; then
        rm -rf $cache_directory
      fi
    done
  else
    rm $dataset_json
    rm -r $dataset_directory
  fi
done
