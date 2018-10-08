# See https://github.com/khayer/aligner_benchmark/blob/master/templates/star.sh
# See https://github.com/khayer/aligner_benchmark/blob/master/templates/novoalign.sh

# Convert BAM to sorted SAM if needed
# samtools view merged.bam | sort -t'.' -k 2n > output.sam

# Fix SAM
# ruby aligner_benchmark/fix_sam.rb output.sam > fixed.sam

# Compare to truth
# perl aligner_benchmark/perl_scripts/compare2truth.pl CIGAR_PATH fixed.sam -noHtag > comp_res.txt

# Infer junctions
# perl aligner_benchmark/perl_scripts/sam2junctions.pl fixed.sam > inferred_junctions.txt

# Junction stats
# perl /Users/hayer/github/aligner_benchmark/perl_scripts/compare_junctions_in_simulated_to_INFERRED_junctions.pl PATH_TO_TRANSCRIPTS JUNCTIONS_CROSSED_PATH inferred_junctions.txt junctions > junctions_stats.txt
