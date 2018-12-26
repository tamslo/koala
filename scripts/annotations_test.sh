# Script to create different annotations with USeq and run novoindex on them

root_directory=$pwd

tests_directory=tests
mkdir -p $tests_directory

function build_annotations() {
  local useq_folder=$1
  local genome_path=$2
  local annotation_path=$3
  local task_directory=$4
  local useq_num=$5
  local transcripts_file=$6
  local splices_file=$7
  local deduplicate=$8
  local radius="96" # read length is 100, minus 4
  java -Xmx48G -jar /opt/$useq_folder/Apps/MakeTranscriptome -f $genome_path \
    -u $annotation_path -n $useq_num -r $radius $deduplicate
  gunzip $splices_file.gz
  gunzip $transcripts_file.gz
  mv $splices_file $task_directory
  mv $transcripts_file $task_directory
}

function mask_genome() {
  local task_directory=$1
  local genome_path=$2
  local annotation_path=$3
  local masked_genome_path=$4
  java -jar /opt/useq/Apps/MaskExonsInFastaFiles -f $genome_path \
    -u $annotation_path -s $masked_genome_path
  cat $masked_genome_path/chr*.fasta > "$masked_genome_path".fa
}

function run_novoindex() {
  local task_directory=$1
  local splices_path=$2
  local transcripts_path=$3
  local masked_genome_path=$2
  novoindex $task_directory/index.nix $splices_path $transcripts_path $masked_genome_path
}

function run_task() {
  local task_name=$1
  local useq_folder=$2
  local genome_path=$3
  local annotation_path=$4
  local task_directory=$5
  local useq_num=$6
  local transcripts_file=$7
  local splices_file=$8
  local deduplicate=$9
  local masked_genome_path=$10

  echo -e "\n$task_name"
  build_annotations $useq_folder $genome_path $annotation_path $task_directory \
    $useq_num $transcripts_file $splices_file $deduplicate

  splices_path=$task_directory/$splices_file
  transcripts_path=$task_directory/$transcripts_file
  run_novoindex $task_directory $splices_path $transcripts_path $masked_genome_path
}

function run_tasks() {
  task_directory=$tests_directory/$2
  genome_path=$3
  annotation_path=$4
  useq_num=$5
  masked_genome_path=$task_directory/masked_$genome_path
  if [ ! -d $task_directory ]; then
    echo -e "Running tasks with $1"
    mkdir $task_directory

    echo -e "\nMask genome"
    mask_genome $task_directory $genome_path $annotation_path $masked_genome_path
    masked_genome_path="$masked_genome_path".fa

    annotation_prefix=${annotation_path::-4}
    useq_annotation_prefix="$annotation_prefix"Rad96Num"$((useq_num/1000))"kMin10
    splices_file="$useq_annotation_prefix"Splices.fasta
    transcripts_file="$useq_annotation_prefix"Transcripts.fasta

    run_task "Current verison, no -s" "useq" $genome_path $annotation_path \
      $task_directory $useq_num $transcripts_file $splices_file $masked_genome_path

    subtask_directory="$task_directory"_with_s
    mkdir $subtask_directory
    run_task "Current verison, with -s" "useq" $genome_path $annotation_path \
      $subtask_directory $useq_num $transcripts_file $splices_file "-s" $masked_genome_path

    subtask_directory="$task_directory"_baruzzo_useq_version
    mkdir $subtask_directory
    run_task "Baruzzo verison, no -s" "useq_baruzzo" $genome_path $annotation_path \
      $subtask_directory $useq_num $transcripts_file $splices_file $masked_genome_path


    subtask_directory="$task_directory"_baruzzo_useq_version_with_s
    mkdir $subtask_directory
    run_task "Baruzzo verison, with -s" "useq_baruzzo" $genome_path $annotation_path \
      $subtask_directory $useq_num $transcripts_file $splices_file "-s" $masked_genome_path


    subtask_directory="$task_directory"_old_useq_version
    mkdir $subtask_directory
    run_task "Old verison, no -s" "useq_old" $genome_path $annotation_path \
      $subtask_directory $useq_num $transcripts_file $splices_file $masked_genome_path


    subtask_directory="$task_directory"_old_useq_version_with_s
    mkdir $subtask_directory
    run_task "Old verison, with -s" "useq_old" $genome_path $annotation_path \
      $subtask_directory $useq_num $transcripts_file $splices_file "-s" $masked_genome_path

  else
    echo -e "$1 already present"
  fi
}

echo ""
run_tasks "Baruzzo annotations for hg19 genome" baruzzo_hg19 hg19 baruzzo.hg19.txt 100000
echo ""
run_tasks "UCSC annotations for hg19 genome" ucsc_hg19 hg19 refFlat.hg19.txt 60000
echo ""
