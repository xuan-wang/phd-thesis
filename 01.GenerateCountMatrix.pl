#! /usr/bin/env perl
use warnings;
use strict;
use autodie;
## Author: Wang Xuan
## Date: 2018-04-03 Tue 16:17:23 CST

use Bio::Seq;
use Bio::SeqIO;
my @samples;
open my $sample,"bcftools query -l $ARGV[1]|";
while (<$sample>) {
    chomp;
    push @samples,$_;
}
close $sample;
my %base=("A"=>"T","G"=>"C","T"=>"A","C"=>"G");
my $fasta_file="$ARGV[0]";
my $fasta_in=Bio::SeqIO->new(-format=>'Fasta',-file=>"$fasta_file");
my %table;
my $chr=$ARGV[2];
while (my $seqobj=$fasta_in->next_seq()) {
    my $id=$seqobj->id;
    last unless $id=~/\d/;
    next if $id <$chr;
    last if $id >$chr;
    print "$id\n";
    open my $vcf,"bcftools view -H -r $id $ARGV[1]|";
    while (<$vcf>) {
        my @tmp=split /\t/;
        my $base3m5;
        if ($tmp[3] eq "A" or $tmp[3] eq "G") {
            $tmp[3]=$base{$tmp[3]};
            $tmp[4]=$base{$tmp[4]};
            $base3m5=$seqobj->trunc($tmp[1]-1,$tmp[1]+1)->revcom->seq();
        }else {
            $base3m5=$seqobj->subseq($tmp[1]-1,$tmp[1]+1);
        }
        for (9..$#tmp){
            $table{$samples[$_-9]}{"$tmp[3]>$tmp[4]:$base3m5"}++ if $tmp[$_]=~/[0|1]\/1/;
        }
    }
    close $vcf;
}
open my $out,">$0.$ARGV[1].$ARGV[2]";
for (("C","T")) {
    for my $alt (sort keys %base) {
        next if $alt eq $_;
        for my $base3 (sort keys %base) {
            for my $base5 (sort keys %base) {
	print $out "\t$_>$alt:${base3}$_${base5}";
            }
        }
    }
}
print $out "\n";
for my $sample (sort keys %table) {
    print $out $sample,"\t";
    for (("C","T")) {
        for my $alt (sort keys %base) {
            next if $alt eq $_;
            for my $base3 (sort keys %base) {
	for my $base5 (sort keys %base) {
	    if (exists $table{$sample}{"$_>$alt:${base3}$_${base5}"}) {
	        print $out "\t$table{$sample}{\"$_>$alt:${base3}$_${base5}\"}";
	    }else {
	        print $out "\t0";
	    }
	}
            }
        }
    }
    print $out "\n";
}
close $out;
