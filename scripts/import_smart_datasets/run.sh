#!/bin/bash

intermediate_directory=data/intermediate

# Get files from Venerea (better do it manually to better select what should be there)
# mkdir -p $intermediate_directory
# rsync -r -P Tamara.Slosarek@192.168.31.56:/opt/hana/pool/user_files/89/RNA_afterTrimming/* $intermediate_directory

concat_file_path() {
  local directory=$1
  local dataset_id=$2
  local direction_indicator=$3
  path=${directory}/${dataset_id}_${direction_indicator}.fq
  echo $path
}

maybe_exctract() {
  local target_path=$(concat_file_path $1 $2 $3)
  if [ ! -f $target_path ]; then
    compressed_path=$target_path.gz
    echo "Unzipping $compressed_path"
    gunzip -c $compressed_path > $target_path
  fi
}

clean_single_files() {
  local directory=$1
  local dataset_id=$2
  local forward_file=$(concat_file_path $1 $2 "1")
  local reverse_file=$(concat_file_path $1 $2 "2")
  for file in $directory/*; do
    if [ $file != $forward_file ] && [ $file != $reverse_file ]; then
      echo "Removing $file"
      rm -f $file
    fi
  done
}

# Concat files that are not concatenated and remove
echo -e "\nPreparing data set import\n"
for directory in $intermediate_directory/*; do
  dataset_id=$(basename $directory)
  maybe_exctract $directory $dataset_id "1"
  maybe_exctract $directory $dataset_id "2"
  clean_single_files $directory $dataset_id
done

# Run python script
echo -e "\nStarting data set import\n"
python3 scripts/import_smart_datasets/import.py
