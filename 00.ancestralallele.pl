#! /usr/bin/env perl
use warnings;
use strict;
use autodie;
## Author: Wang Xuan
## Date: 2022-08-13 Sat 15:59:09 CST

open my $in,"<$ARGV[0]";
while (<$in>) {
    if (/#/){
        print;
        next;
    }
    chomp;
    my $line=$_;
    my @tmp=split /\t/,$line;
    my $out=`bcftools view -H -r $tmp[0]:$tmp[1] ../../AndeanFox.vcf.gz`;
    chomp $out;
    my @tmp_out=split /\t/,$out;
    if ($tmp_out[9]=~/1\/1/) {
        $line=~s/0\/0/A\/A/g if($line=~/0\/0/);
        $line=~s/1\/1/0\/0/g if($line=~/1\/1/);
        $line=~s/A\/A/1\/1/g if($line=~/A\/A/);
    }
    print $line,"\n";
}
close $in;
