#! /usr/bin/env perl
use warnings;
use strict;
use autodie;
## Author: Wang Xuan
## Date: 2022-08-02 Tue 14:28:21 CST

die("usage: $0 pop\n") unless $ARGV[0];
my $pop=$ARGV[0];
open my $out,">$0.$pop.out";
print $out "chr\tpos\tLikelihood\talpha\tStartPos\tEndPos\n";
for my $i (1..38) {
    open my $in,"<SweeD_Report.$pop.chr$i.10k";
    <$in>;<$in>;
    while (<$in>) {
        my @tmp=split /\t/;
        next if $tmp[0]<1000000;
        print $out "$i\t$_";
    }
    close $in;
}
close $out;
