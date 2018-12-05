# Extract annotations (only for standard chromosomes)
if [ ! -f data/annotations/hg19.txt ]; then
  echo "Parsing annotations"
  gunzip -c data/annotations/refFlat.hg19.txt.gz | grep -v 'chr.*_' > data/annotations/hg19.txt
  gunzip -c data/annotations/refFlat.hg38.txt.gz | grep -v 'chr.*_' > data/annotations/hg38.txt
  echo "Done"
else
  echo "Annotations already parsed"
fi

# Download references for single chromosomes (and delete non-standard)
if [ ! -d data/references/hg19 ]; then
  echo "Downloading reference FASTAs per chromosome"
  current_dir=$pwd
  wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/*.fa.gz -P data/references/hg19
  wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/*.fa.gz -P data/references/hg38

  cd ~/code/data/references/hg19
  rm -f *_*.fa.gz

  cd ~/code/data/references/hg38
  rm -f *_*.fa.gz
  cd $current_dir
  echo "Done"
else
  echo "FASTAs per chromosome already present"
fi
