#! /usr/bin/env perl
use warnings;
use strict;
use autodie;
## Author: Wang Xuan
## Date: 2018-08-14 Tue 22:01:45 CST

# my @results=<*.Q.withid.structure>;
open my $out,">$0.out";
for my $k (2..10) {
    open my $in,"<GDJCYTWLFDOGaDNACTVT_F152inds.SnpGap3.snp.filter.PASS.biallelic.miss10.unphase.concCTVT_F.vcf.gz.ld02.$k.Q.withid.structure";
    while (<$in>) {
        chomp;
        my @tmp=split /\s+/;
        for my $v (1..$k) {
            my $n=$tmp[$v]*100000;
            for my $m (1..$n) {
	print $out "$tmp[0]\tA$v\tK=$k\n";
            }
        }
    }
    close $in;
}
close $out;
