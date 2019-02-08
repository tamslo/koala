#!/bin/bash

intermediate_directory=data/intermediate/

# Get files from Venerea
if [ ! -d $intermediate_directory ]; then
  mkdir $intermediate_directory
  rsync -r -P Tamara.Slosarek@192.168.31.56:/opt/hana/pool/user_files/89/RNA_afterTrimming/*
fi

# Run python script
python scripts/import_datasets/import.py
