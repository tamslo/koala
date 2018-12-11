# Basically from http://www.novocraft.com/documentation/novoalign-2/novoalign-user-guide/rnaseq-analysis-mrna-and-the-spliceosome/
#!/usr/bin/perl

use strict;
use Bio::SeqIO;
my %unique;
my $file = $ARGV[0];
my $seqio = Bio::SeqIO->new(-file => $file, -format => "fasta");
my $outseq = Bio::SeqIO->new(-file => ">$file.uniq", -format => "fasta");

while(my $seqs = $seqio->next_seq) {
    my $id = $seqs->display_id;
    my $seq = $seqs->seq;
    unless(exists($unique{$seq})) {
        $outseq->write_seq($seqs);
        $unique{$seq} +=1;
    }
}
