#=== Merging the technical replicates ===

echo "********** combine the replicates samples  ************"

for dir in $(find $Samples -mindepth 1 -maxdepth 1 -type d ); do
    echo " the base sample name is $dir"

    a=${dir}/$(basename $dir)-a_*_R1*.fastq.gz
    b=${dir}/$(basename $dir)-b_*_R1*.fastq.gz
    c=${dir}/$(basename $dir)-c_*_R1*.fastq.gz

    echo "merge the forward reads the technical replicates"
    cat $a $b $c > ${dir}/$(basename $dir)_R1.fastq.gz

    aa=${dir}/$(basename $dir)-a_*_R2*.fastq.gz
    bb=${dir}/$(basename $dir)-b_*_R2*.fastq.gz
    cc=${dir}/$(basename $dir)-c_*_R2*.fastq.gz

    echo "merge the reversed reads the technical replicates"
    cat $aa $bb $cc > ${dir}/$(basename $dir)_R2.fastq.gz
done
