mkdir -p data/annotations
current_dir=$pwd
cd data/annotations

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

# Extract present BED files from UCSC Table Browser
for zipped_file in *.gz; do
  unzipped_file=${zipped_file:0:${#zipped_file}-3}
  if [ ! -f $unzipped_file ]; then
    echo "Extracting $unzipped_file"
    gunzip -c $zipped_file > $unzipped_file
    if [[ $unzipped_file == *_coding_exons.bed ]]; then
      echo "Getting first three columns"
      tmp_file="bed.tmp"
      mv $unzipped_file $tmp_file
      awk -v OFS='\t' '{print $1,$2,$3}' $tmp_file > $unzipped_file
      rm $tmp_file
    fi
    if [[ $unzipped_file == "hg19_coding_exons.bed" ]]; then
      echo "Fixing chromosome names"
      tmp_file="hg19_coding_exons.tmp"
      mv $unzipped_file $tmp_file
      cat $tmp_file | sed 's/^chr//' > $unzipped_file
      rm $tmp_file
    fi
    echo "Done."
  else
    echo "Already extracted $unzipped_file"
  fi
done

cd $pwd
