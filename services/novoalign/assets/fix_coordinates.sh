sam_file_path="$1"

# -a sets the maximum alignment score
# -n sets the number of locations to which each read may align
# -u saves unmapped reads
java -jar /opt/useq/Apps/SamTranscriptomeParser -f  $sam_file_path -a 50000 -n 100 -u
