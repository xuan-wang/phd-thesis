#! /usr/bin/env perl
use warnings;
use strict;
use autodie;
## Author: Wang Xuan
## Date: 2018-01-29 Mon 16:41:50 CST

open my $vcf,"bcftools view -r $ARGV[2] -HG $ARGV[0]|cut -f 1-5|";
open my $germline,">$0.$ARGV[2].GL.vcf";
open my $s2,">$0.$ARGV[2].s2.vcf";
open my $s2multi,">$0.$ARGV[2].s2.multi";
while (<$vcf>) {
    chomp;
    s/chr//;
    my @tmp=split /\t/;
    my $alt=$tmp[4];
    my $mono=`tabix $ARGV[1] $ARGV[2]:$tmp[1]-$tmp[1]`;
    next if $mono eq '';
    chomp $mono;
    @tmp=split /\t/,$mono;
    my %monogeno;
    $monogeno{$_}++ for @tmp[3..$#tmp];
    for my $monogeno (keys %monogeno) {
        next if $monogeno eq ".";
        my @monoallele=split //,$monogeno;
        if (@monoallele == 1) {
            if ($monoallele[0] eq $tmp[2]) {
	print $germline $tmp[0],"\t",$tmp[1],"\t.\t",$tmp[2],"\t.\t.\tPASS\t.\tGT\t0/0\n";
            }elsif ($monoallele[0] eq $alt) {
	print $germline $tmp[0],"\t",$tmp[1],"\t.\t",$tmp[2],"\t",$alt,"\t.\tPASS\t.\tGT\t1/1\n";
            }else {
	print $s2 $tmp[0],"\t",$tmp[1],"\t.\t",$tmp[2],"\t",$monoallele[0],"\t.\tPASS\t.\tGT\t1/1\n";	# candidate somatic 2
            }
        }elsif (@monoallele == 2) {
            my %monoallele;
            $monoallele{$_}++ for @monoallele;
            if (exists $monoallele{$tmp[2]}) {
	for my $allele (@monoallele) {
	    next if $allele eq $tmp[2];
	    if ($allele eq $alt) {
	        print $germline $tmp[0],"\t",$tmp[1],"\t.\t",$tmp[2],"\t",$alt,"\t.\tPASS\t.\tGT\t0/1\n";
	    }else {
	        print $s2 $tmp[0],"\t",$tmp[1],"\t.\t",$tmp[2],"\t",$allele,"\t.\tPASS\t.\tGT\t0/1\n"; # candidate somatic 2
	    }
	}
            }else {
	my $twos2alt=join(",",@monoallele);
	print $s2 $tmp[0],"\t",$tmp[1],"\t.\t",$tmp[2],"\t",$twos2alt,"\t.\tPASS\t.\tGT\t1/2\n"; # candidate somatic 2
            }
        }else {
            print $s2multi $tmp[0],"\t",$tmp[1],"\t",$tmp[2],"\t",$monogeno,"\n"; # candidate multiallelic somatic 2
        }
    }
}
close $vcf;
close $s2;
close $s2multi;
close $germline;
