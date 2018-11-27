reference_base_path="$1"
index_path="$2"
reference_id="$3"
read_length="$4"
reference_path="$5"

if [ -f /annotations/$reference_id.txt ]; then
  # Build masked genome
  masked_genome_path=/$reference_base_path/"$reference_id"_masked
  if [ ! -f $masked_genome_path ]; then
    java -jar /opt/useq/Apps/MaskExonsInFastaFiles -f /references/$reference_id -u /annotations/$reference_id.txt -s $masked_genome_path
  fi

  # Build annotation files
  rad=$((read_length - 4))
  num=60000
  min=10 # default
  prefix="$reference_id"Rad"$rad"Num$((num / 1000))kMin"$min"
  known_junctions_file="$prefix"Splices.fasta.gz
  known_junctions_path=/$reference_base_path/$known_junctions_file
  theoretical_junctions_file="$prefix"Transcripts.fasta.gz
  theoretical_junctions_path=/$reference_base_path/$theoretical_junctions_file

  if [ ! -f $known_junctions_path ]; then
    # -r sequence lenght radius, set to read length - 4bp
    # -n reduces the maximal number of splices per gene to be faster
    # -m sets maximal minutes to process each gene's splices before interrupting
    # -s skips subsequent occurrences of splices with the same coordinates, memory intensive
    java -jar /opt/useq/Apps/MakeTranscriptome -f /references/$reference_id -u /annotations/$reference_id.txt -r $rad -n $num -m $min -s
    mv /annotations/$known_junctions_file $known_junctions_path
    mv /annotations/$theoretical_junctions_file $theoretical_junctions_path
  fi

  # Build index
  novoindex -n $reference_id /$index_path $masked_genome_path $known_junctions_path $theoretical_junctions_path
else
  novoindex -n $reference_id /$index_path /$reference_path
fi
