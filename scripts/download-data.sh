mkdir -p data/annotations

# Extract annotations
if [ ! -f data/annotations/hg19.txt ]; then
  echo "Extracting annotations"
  gunzip -c data/annotations/hg19.gtf.gz > data/annotations/hg19.gtf
  gunzip -c data/annotations/hg38.gtf.gz > data/annotations/hg38.gtf
  echo "Done"
else
  echo "Annotations already extracted"
fi
