mkdir -p data/annotations
current_dir=$pwd
cd data/annotations

function minimize_bed() {
  echo "Getting first three columns"
  local in_file=$1
  local out_file=$2
  awk -v OFS='\t' '{print $1,$2,$3}' $in_file > $out_file
}

function fix_hg19_chromosomes() {
  echo "Fixing chromosome names"
  local in_file=$1
  local out_file=$2
  cat $in_file | sed 's/^chr//' > $out_file
}

# Get annotations
if [ ! -f hg19.gtf ]; then
  echo "Getting annotations for hg19..."
  wget -q ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
  gunzip -c gencode.v19.annotation.gtf.gz > hg19.gtf
  rm gencode.v19.annotation.gtf.gz
  echo "Done."
else
  echo "Annotations for hg19 already present"
fi
if [ ! -f hg38.gtf ]; then
  echo "Getting annotations for hg38..."
  wget -q ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz
  gunzip -c gencode.v29.annotation.gtf.gz > hg38.gtf
  rm gencode.v29.annotation.gtf.gz
  echo "Done."
else
  echo "Annotations for hg38 already present"
fi

# Get RNA editing sites for hg19 (none available so far for hg38)
if [ ! -f hg19_editing_sites.bed ]; then
  echo "Getting RNA editing sites for hg19..."
  wget -q http://lilab.stanford.edu/GokulR/database/Human_AG_all_hg19.bed
  fix_hg19_chromosomes Human_AG_all_hg19.bed hg19_editing_sites.bed
  rm Human_AG_all_hg19.bed
  echo "Done."
else
  echo "RNA editing sites for hg19 already present"
fi

# Extract present BED files from UCSC Table Browser
for zipped_file in *.gz; do
  unzipped_file=${zipped_file:0:${#zipped_file}-3}
  if [ ! -f $unzipped_file ]; then
    echo "Extracting $unzipped_file"
    gunzip -c $zipped_file > $unzipped_file
    if [[ $unzipped_file == *.bed ]]; then
      tmp_file="bed.tmp"
      mv $unzipped_file $tmp_file
      minimize_bed $tmp_file $unzipped_file
      rm $tmp_file
      if [[ $unzipped_file == hg19_* ]]; then
        tmp_file="hg19_coding_exons.tmp"
        mv $unzipped_file $tmp_file
        fix_hg19_chromosomes $tmp_file $unzipped_file
        rm $tmp_file
      fi
    fi
    echo "Done."
  else
    echo "Already extracted $unzipped_file"
  fi
done

cd $pwd
