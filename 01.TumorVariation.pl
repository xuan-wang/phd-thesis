#! /usr/bin/env perl
use warnings;
use strict;
use autodie;
## Author: Wang Xuan
## Date: 2023-09-09 Sat 10:43:15 CST

open my $file,"pigz -dc $ARGV[0]|";
my $chr=$ARGV[1];
my $sample=`basename $ARGV[0] .$chr.baseCN.gz`;
chomp $sample;
my %code=(3=>"A",4=>"C",5=>"G",6=>"T");
open my $out,"|bgzip -@ 3 > $0.$sample.$chr.vcf.gz";
print $out "##fileformat=VCFv4.2
\##FILTER=<ID=PASS,Description=\"All filters passed\">
\##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
\#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sample\n";
open my $multi,">$0.$sample.$chr.multi";
while (<$file>) {
     my @tmp=split /\t/;
     my %tumorallele;
     chomp $tmp[7];
     if ($tmp[7]>1) {
         for (3..6) {
             $tumorallele{$code{$_}}++ if $tmp[$_]>0;
         }
         my @tumorallele=keys %tumorallele;
         if (exists $tumorallele{$tmp[2]} and scalar(@tumorallele)==1) {
             #next;
             print $out "chr$chr\t$tmp[1]\t.\t$tmp[2]\t.\t.\tPASS\t.\tGT\t0/0\n";
         }elsif (exists $tumorallele{$tmp[2]} and scalar(@tumorallele)==2) {
             for my $allele (@tumorallele) {
	 next if $allele eq $tmp[2];
	 print $out "chr$chr\t$tmp[1]\t.\t$tmp[2]\t$allele\t.\tPASS\t.\tGT\t0/1\n"
             }
         }elsif (scalar(@tumorallele)==1 and $tumorallele[0] ne $tmp[2]) {
             print $out "chr$chr\t$tmp[1]\t.\t$tmp[2]\t$tumorallele[0]\t.\tPASS\t.\tGT\t1/1\n";
         }else {
             my $gt=join("",@tumorallele);
             print $multi "chr$chr\t$tmp[1]\t.\t$tmp[2]\t$gt\n";
         }
     }
 }
close $file;
close $out;
close $multi;
