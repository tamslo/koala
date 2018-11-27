reference_base_path="$1"
index_path="$2" # needs to be absolute
reference_id="$3"
read_length="$4"

# Build masked genome
masked_genome_path=/$reference_base_path"$reference_id"_masked
if [ ! -f $masked_genome_path ]; then
  java -jar /opt/useq/Apps/MaskExonsInFastaFiles -f /references/$reference_id -u /annotations/$reference_id.txt -s $masked_genome_path
fi
# TODO how are files called?
known_junctions_file=""
theoretical_junctions_file=""
# if [ ! -f $known_junctions_file ]; then
  # -n reduces the maximal number of splices per gene to be faster
  java -jar /opt/useq/Apps/MakeTranscriptome -f /references/$reference_id -u /annotations/$reference_id.txt -r $((read_length - 4)) -n 60000 -s
# fi

# novoindex -n $reference_id $index_path $masked_genome_path $known_junctions_file $theoretical_junctions_file
