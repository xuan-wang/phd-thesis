#! /usr/bin/env perl
use warnings;
use strict;
use autodie;
## Author: Wang Xuan
## Date: 2023-12-28 Thu 15:37:42 CST

my $chr=$ARGV[0];
open my $samplein,"bcftools query -l ../04.StrandBias/04.StrandBias.pl.vcf.gz|";
my @samples;
while (<$samplein>) {
    chomp;
    push (@samples,$_);
}
close $samplein;
open my $vcf,"bcftools view -r $chr -H ../04.StrandBias/04.StrandBias.pl.vcf.gz|";
open my $out,"|bgzip > $0.$chr.vcf.gz";
my $head=`bcftools view -h ../02.TumorVariation.pl.5CTVT.v2-v5.noHostPanel.vcf.gz`;
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
    my $existhostless20=0;
    for (my $i=9;$i<=$#tmp;$i++) {
        next unless ($tmp[$i] eq "1/1" or $tmp[$i] eq "0/1");
        $samples[$i-9]=~/T\.(.*)/;
        my $id=$1;
        my $acgt=`tabix /data/wanggd_group/wangxuan/wangxuan/projects/ctvt/14.DNAandRNA2CanFam/07.Sequenza.AutosomeXY/04.makeacgt/CFAM.H$id.acgt.gz $nochr:$tmp[1]-$tmp[1]`;
        next unless $acgt;
        my @tmpacgt=split /\t/,$acgt;
        if ($id eq "609" and $tmpacgt[3]<20) {
            $existhostless20=1;
            last;
        }elsif ($id ne "609" and $tmpacgt[3]<10) {
            $existhostless20=1;
            last;
        }
    }
    print $out $line,"\n" if $existhostless20==0;
}
close $vcf;
close $out;
