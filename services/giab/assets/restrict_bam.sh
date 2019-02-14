dataset_folder=/data

minimum_coverage=$1
folder_prefix=${2:-""}

# Constants

out_bam_name="Out.bam"
full_bam_name="Out.full.bam"

function get_python_command() {
  # From https://stackoverflow.com/questions/6141581/
  # detect-python-version-in-shell-script
  local python_version=$(python -c 'import sys; print(sys.version_info[:][0])')
  if [ $python_version -eq "3" ]; then
    echo python
  else
    echo python3
  fi
}

python=$(get_python_command)

# Helper functions

# Check if Out.full.bam exists (this script already was executed), otherwise
# set original Out.bam as original file
function get_original_bam_path() {
  out_bam_path=$1
  full_bam_path=$2
  if [ -f $full_bam_path ]; then
    echo $full_bam_path
  else
    echo $out_bam_path
  fi
}

# Create coverage BED file (if not already present)
function create_coverage_bed() {
  alignment_path=$1
  bam_path=$2
  bed_path=$3

  if [ ! -f $bed_path ]; then
    coverage_path=${alignment_path}/Out.base_coverage

    if [ ! -f $coverage_path ]; then
      echo "Creating base coverage file"
      bedtools genomecov -d -ibam $bam_path > $coverage_path
    else
      echo "Base coverage file already present"
    fi

    echo "Creating BED file from base coverage"
    $python base_coverage_to_bed.py $coverage_path $minimum_coverage $bed_path
    rm ${coverage_path}.${minimum_coverage}-filtered
  else
    echo "Coverage BED file already present"
  fi
}

# Intersect full BAM with coverage BED
function create_coverage_bam() {
  original_bam_path=$1
  coverage_bed_path=$2
  coverage_bam_path=$3

  bedtools intersect -a $original_bam_path -b $coverage_bed_path > $coverage_bam_path
}

function process_alignment() {
  alignment_path=$1
  coverage_bam_path=${alignment_path}/Out.coverage_${minimum_coverage}.bam

  out_bam_path=$alignment_path/$out_bam_name
  full_bam_path=$alignment_path/$full_bam_name

  if [ ! -f $coverage_bam_path ]; then
    original_bam_path=$(get_original_bam_path $out_bam_path $full_bam_path)
    coverage_bed_path=${alignment_path}/Out.coverage_${minimum_coverage}.bed
    create_coverage_bed $alignment_path $original_bam_path $coverage_bed_path
    create_coverage_bam $original_bam_path $coverage_bed_path $coverage_bam_path
  else
    echo "Coverage BAM file already present"
  fi

  # Rename original BAM if necessary
  if [ ! -f $full_bam_path ]; then
    echo "Archiving original BAM file"
    mv $out_bam_path $full_bam_path
  fi

  # Create link to coverage BAM named Out.bam
  ln -s -f $coverage_bam_path $out_bam_path
}

# Script execution

# Test if minimum_coverage is a number (from https://stackoverflow.com/
# questions/59838/check-if-a-directory-exists-in-a-shell-script)
re='^[0-9]+$'
if ! [[ $minimum_coverage =~ $re ]] ; then
   echo "[ERROR] Please provide a number as input for minimum coverage" >&2
   exit 1
fi

dataset_counter=0
directory=${dataset_folder}/${folder_prefix}*
echo -e "\n[INFO] Starting BAM restriction with minimum_coverage $minimum_coverage for $directory"

# Traverse data sets
for dataset in $directory; do
  if [ -d $dataset ]; then
    ((dataset_counter++))
    echo -e "\n[PROGRESS] Processing data set $(basename $dataset)"
    for reference_genome in ${dataset}/*; do
      for alignment in ${reference_genome}/*; do
        if [ -d $alignment ]; then
          info_text="Processing $(basename $alignment) alignment "
          info_text+="($(basename $reference_genome))"
          echo -e "\n$info_text"
          process_alignment $alignment
        fi
      done
    done
  fi
done

echo -e "\n[INFO] Done. Processed $dataset_counter data set(s)\n"
