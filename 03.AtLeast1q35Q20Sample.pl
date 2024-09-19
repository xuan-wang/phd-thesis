#! /usr/bin/env perl
use warnings;
use strict;
use autodie;
## Author: Wang Xuan
## Date: 2023-12-27 Wed 20:45:24 CST

my $chr=$ARGV[0];
open my $samplein,"bcftools query -l ../02.TumorVariation.pl.5CTVT.v2-v5.noHostPanel.vcf.gz|";
my @samples;
while (<$samplein>) {
    chomp;
    push (@samples,$_);
}
close $samplein;
open my $vcf,"bcftools view -r $chr -H ../02.TumorVariation.pl.5CTVT.v2-v5.noHostPanel.vcf.gz|";
my $head=`bcftools view -h ../02.TumorVariation.pl.5CTVT.v2-v5.noHostPanel.vcf.gz`;
open my $out,"|bgzip > $0.$chr.vcf.gz";
print $out $head;
my %bases;
$bases{A}=4;
$bases{C}=5;
$bases{G}=6;
$bases{T}=7;
while (my $line=<$vcf>) {
    chomp $line;
    my @tmp=split /\t/,$line;
    $tmp[0]=~/chr(.*)/;
    my $nochr=$1;
    for (my $i=9;$i<=$#tmp;$i++) {
        my $acgt=`tabix /data/wanggd_group/wangxuan/wangxuan/projects/ctvt/14.DNAandRNA2CanFam/07.Sequenza.AutosomeXY/04.makeacgt/$samples[$i-9].q35Q20.acgt.gz $nochr:$tmp[1]-$tmp[1]`;
        next unless $acgt;
        my @tmpacgt=split /\t/,$acgt;
        if ($tmpacgt[$bases{$tmp[4]}]>=2) {
            print $out $line,"\n";
            last;
        }
    }
}
close $vcf;
close $out;
