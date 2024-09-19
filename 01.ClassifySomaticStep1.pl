#! /usr/bin/env perl
use warnings;
use strict;
use autodie;
## Author: Wang Xuan
## Date: 2018-01-23 Tue 21:40:32 CST

# usage: $0 *baseCN chr
my @baseCN=@ARGV;
my $chr=pop @baseCN;
print "$chr\n";
my %base;
my %samples;
my %ref;
my %code=(3=>"A",4=>"C",5=>"G",6=>"T");
for my $baseCN (@baseCN) {
    my $sample=`basename $baseCN .$chr.baseCN.gz`;
    chomp $sample;
    $samples{$sample}++;
    open my $file,"pigz -dc $baseCN|";
    while (<$file>) {
        my @tmp=split /\t/;
        my $gt;
        chomp $tmp[7];
        if ($tmp[7]>1) {
            for (3..6) {
	$gt.=$code{$_} if $tmp[$_]>0;
            }
            $base{$tmp[1]}{$sample}=$gt if ($gt);
            $ref{$tmp[1]}=$tmp[2];
        }
    }
    close $file;
}
open my $s1,"|bgzip >$0.$chr.s1.gz";
open my $mono,"|bgzip >$0.$chr.mono.gz";
open my $miss2,"|bgzip >$0.$chr.miss2plus.gz";
my $head="chr\tpos\tref";
for my $sample (sort keys %samples) {
    $head.="\t$sample";
}
print $s1 $head,"\n";
print $mono $head,"\n";
print $miss2 $head,"\n";
for my $pos (sort{$a<=>$b} keys %base) {
    my %count;
    my $line;
    my $miss=0;
    for my $sample (sort keys %samples) {
        if (exists $base{$pos}{$sample}) {
#            print $pos,$sample,"\n" unless($base{$pos}{$sample});
            $line.="\t$base{$pos}{$sample}";
            $count{$base{$pos}{$sample}}++;
        }else {
            $line.="\t.";
            $miss++;
        }
    }
    if ($miss<2) {
        if (keys %count == 1) {
            print $mono "chr$chr","\t",$pos,"\t",$ref{$pos},$line,"\n"; # monomorphism genotype
        }elsif (keys %count >1) {
            print $s1 "chr$chr","\t",$pos,"\t",$ref{$pos},$line,"\n"; # somatic step1
        }
    }else {
        print $miss2 "chr$chr","\t",$pos,"\t",$ref{$pos},$line,"\n"; # at least 2 miss
    }
}
close $s1;
close $mono;
close $miss2;
