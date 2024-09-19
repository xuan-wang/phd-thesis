#! /usr/bin/env perl
use warnings;
use strict;
use autodie;
## Author: Wang Xuan
## Date: 2018-04-03 Tue 22:37:01 CST

my @result_chr=@ARGV;
my %table;
my $head;
for my $chr (@result_chr) {
    open my $in,"<$chr";
    $head=<$in>;
    while (<$in>) {
        chomp;
        my @tmp=split /\t/;
        for (2..$#tmp) {
            $table{$tmp[0]}{$_}+=$tmp[$_];
        }
    }
    close $in;
}
print $head;
for my $sample (sort keys %table) {
    print $sample;
    for my $type (sort{$a<=>$b} keys $table{$sample}) {
        print "\t$table{$sample}{$type}";
    }
    print "\n";
}
