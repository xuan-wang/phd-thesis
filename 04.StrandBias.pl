#! /usr/bin/env perl
use warnings;
use strict;
use autodie;
## Author: Wang Xuan
## Date: 2023-12-28 Thu 15:37:42 CST

my $chr=$ARGV[0];
open my $samplein,"bcftools query -l ../03.AtLeast1q35Q20Sample/03.AtLeast1q35Q20Sample.pl.vcf.gz|";
my @samples;
while (<$samplein>) {
    chomp;
    push (@samples,$_);
}
close $samplein;
open my $vcf,"bcftools view -r $chr -H ../03.AtLeast1q35Q20Sample/03.AtLeast1q35Q20Sample.pl.vcf.gz|";
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
    my $totalreadsnumberalltumors;
    my $reversereadsnumberalltumors;
    for (my $i=9;$i<=$#tmp;$i++) {
        my $acgt=`tabix /data/wanggd_group/wangxuan/wangxuan/projects/ctvt/14.DNAandRNA2CanFam/07.Sequenza.AutosomeXY/04.makeacgt/$samples[$i-9].acgt.gz $nochr:$tmp[1]-$tmp[1]`;
        next unless $acgt;
        my @tmpacgt=split /\t/,$acgt;
        $totalreadsnumberalltumors+=$tmpacgt[$bases{$tmp[4]}];
        my @rv=split /:/,$tmpacgt[8];
        my $rvindex=$bases{$tmp[4]}-4;
        $reversereadsnumberalltumors+=$rv[$rvindex];
    }
    if ($totalreadsnumberalltumors <= 15) {
        if ($reversereadsnumberalltumors >=3 and ($totalreadsnumberalltumors-$reversereadsnumberalltumors)>=3) {
            print $out $line,"\n";
        }
    }else {
        if (($reversereadsnumberalltumors/$totalreadsnumberalltumors)>=0.2 and (($totalreadsnumberalltumors-$reversereadsnumberalltumors)/$totalreadsnumberalltumors)>=0.2) {
            print $out $line,"\n";
        }
    }
}
close $vcf;
close $out;
