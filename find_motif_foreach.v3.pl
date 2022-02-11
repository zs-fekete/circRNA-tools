#!/usr/bin/perl -w
use strict;

my $file=shift or die "Missing input!\n\nUsage:\nfind_motif_foreach.pl circRNA-star.output.txt";

my $gtf = "OryCun2.104.exon-only.gtf";
my $fasta = "/molbio/flatfiles/genomes/oc/Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa";
my $treshold = 3;
my $kmer = 7;

print "ID\tk-mer\tk-mer_count\tk-mer-exon_length_ratio\texon_length\tpotential_length\n";

open (MAIN, $file);
while(<MAIN>){
    chomp;
# handle header
    my @main_cell = split (/\t/,$_);
    my $max_length;
# converting into .bed format
    my ($chr, $tmp) = split (/:/, $main_cell[0]);
    my ($start, $end) = split (/-/, $tmp);
    open(TMP1, '>', "tmp1.bed") or die $!;
    if ($start > $end){
	print TMP1 "$chr\t$end\t$start";
	$max_length = $start-$end;
    }
    else {
	print TMP1 "$chr\t$start\t$end";
	$max_length = $end-$start;
    }
    close TMP1;
	
    system "bedtools intersect -b tmp1.bed -a $gtf | cut -f 1,4,5 | bedtools sort -i - | bedtools merge -i - > tmp2.bed";
    system "bedtools getfasta -fi $fasta -bed tmp2.bed >tmp3.fa";
    system "jellyfish count -m $kmer -t 4 tmp3.fa -s 1M";
    my @array =  `jellyfish dump -c mer_counts.jf`;
    my $hwm = 0;
    my $seq;
    foreach my $element (@array){ 
	chomp;
	my @cell = split (" ", $element);
	if ($cell[1] > $hwm){
	    $seq = $cell[0];
	    $hwm = $cell[1];
	}
    }
# calculating length
    my $exon_length =0;
    open (LEN, "tmp3.fa");
    while (<LEN>){
	next if /^>/;
	chomp;
	$exon_length += length($_);
    }
    close LEN;
    my $ratio;
# remove low count
    if ($hwm < $treshold){
	$seq = "NA";
	$hwm = "NA";
	$ratio = "NA";
    }
    else {
	$ratio = $hwm * $kmer / $exon_length;
    }

# print output
    foreach my $i (@main_cell){
	print "$i\t";
    }
    print "$seq\t$hwm\t$ratio\t$exon_length\t$max_length\n";
}
