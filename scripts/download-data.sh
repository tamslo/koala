# Extract annotations
if [ ! -f data/annotations/hg19.txt ]; then
  gunzip -c data/annotations/refFlat.hg19.txt.gz > data/annotations/hg19.txt
  gunzip -c data/annotations/refFlat.hg38.txt.gz > data/annotations/hg38.txt
fi

# Download references for single chromosomes
if [ ! -d data/references/hg19 ]; then
  wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/* -P data/references/hg19
  wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/* -P data/references/hg38
fi
