#! /usr/bin/env perl
use warnings;
use strict;
use autodie;
## Author: Wang Xuan
## Date: 2017-12-27 Wed 21:21:14 CST
use Math::Round;

my $chrnum=$ARGV[-1];
print $chrnum,"\n";

open my $hostdel,"bcftools view -H -v indels -r $chrnum $ARGV[0]|vcf2bed -n -d|cut -d \";\" -f 1|";
print "Host del\n";
my %hostdel;
while (<$hostdel>) {
    chomp;
    my @tmp=split /\t/;
    if ($tmp[-1] eq "AC=2"){
        for (my $i=($tmp[1]+1+length($tmp[6]));$i<($tmp[1]+1+length($tmp[5]));$i++) {
            $hostdel{$i}{c}+=2;
            $hostdel{$i}{g}=1;
        }
    }else {
        for (my $i=($tmp[1]+1+length($tmp[6]));$i<($tmp[1]+1+length($tmp[5]));$i++) {
            $hostdel{$i}{c}++;
            $hostdel{$i}{g}=1;
        }
    }
}
close $hostdel;
open my $tdel,"bcftools view -H -r $chrnum -v indels $ARGV[1]|vcf2bed -n -d|" or die ("Can't open $ARGV[1]\n");
my %tumordel;
print "Tumor m2.del\n";
while (<$tdel>) {
    my @tmp=split /\t/;
    my @format=split /:/,$tmp[10];
    $format[1]=~/,(\d+)/;
    my $tumor_num=$1;
    @format=split /:/,$tmp[11];
    $format[1]=~/(\d+),(\d+)/;
    my $host_c=0;
    $host_c=round($2*2/($1+$2)) unless ($1+$2)==0;
    for (my $i=($tmp[1]+1+length($tmp[6]));$i<($tmp[1]+1+length($tmp[5]));$i++) {
        $tumordel{$i}{c}+=$tumor_num;
        $tumordel{$i}{g}=1;
        if ($host_c != 0) {
            $hostdel{$tmp[1]}{c}+=$host_c unless exists $hostdel{$tmp[1]}{g};
        }
    }
}
close $tdel;
open my $tbdel2,"bcftools view -H -M 2 -r $chrnum $ARGV[2]|vcf2bed -n -d|" or die ("Can't open $ARGV[2]\n");
print "Tumor bcftools.del2\n";
while (<$tbdel2>) {
    my @tmp=split /\t/;
    my @format=split /:/,$tmp[10];
    $format[2]=~/,(\d+)/;
    my $tumor_num=$1;
    for (my $i=($tmp[1]+1+length($tmp[6]));$i<($tmp[1]+1+length($tmp[5]));$i++) {
        unless (exists $tumordel{$i}) {
            $tumordel{$i}{c}+=$tumor_num;
            $tumordel{$i}{g}=1;
        }
    }
}
close $tbdel2;
open my $tbdelm,"bcftools view -H -m 3 -r $chrnum $ARGV[2]|" or die ("Can't open $ARGV[2]\n");
print "Tumor bcftools.delm\n";
while (<$tbdelm>) {
    my @tmp=split /\t/;
    my $reflength=length($tmp[3]);
    my @alt=split /,/,$tmp[4];
    my %index;
    for (my $i=0;$i<$#alt;$i++) {
        my $altlength=length($alt[$i]);
        $index{$i}=$altlength if ($altlength<$reflength);
    }
    my @format=split /:/,$tmp[9];
    my @AD=split /,/,$format[2];
    for my $index (sort keys %index) {
        for (my $i=($tmp[1]+$index{$index});$i<($tmp[1]+$reflength);$i++) {
            unless (exists $tumordel{$i}{g}) {
	$tumordel{$i}{c}+=$AD[$index+1];
            }
        }
    }
}
close $tbdelm;
open my $hcnv,"tabix $ARGV[3] $chrnum|" or die("Can't open $ARGV[3]\n");
print "Host cnv\n";
my %hostcnv;
while (<$hcnv>) {
    my @tmp=split /\t/;
    chomp $tmp[3];
    for (my $i=($tmp[1]+1);$i<=$tmp[2];$i++) {
        $hostcnv{$i}=$tmp[3];
    }
}
close $hcnv;
open my $cnv,"tabix $ARGV[4] $chrnum|" or die ("Can't open $ARGV[4]\n");
print "Tumor cnv\n";
my %cnv;
while (<$cnv>) {
    my @tmp=split /\t/;
    for (my $i=$tmp[1];$i<$tmp[2];$i++) {
        $cnv{$i}=$tmp[9];
    }
}
close $cnv;
open my $hatcg,"tabix $ARGV[5] $chrnum|" or die ("Can't open $ARGV[5]\n");
print "Host atcg\n";
my %amb;
my %hostgeno;
while (<$hatcg>) {
    my @tmp=split /\t/;
    my $sum=$tmp[4]+$tmp[5]+$tmp[6]+$tmp[7];
    if ($sum>1) {
        $hostgeno{$tmp[1]}{sum}=$sum;
    }else {
        $amb{$tmp[1]}="host total raw reads num less two";
        next;
    }
    $hostgeno{$tmp[1]}{base}{A}=$tmp[4] if $tmp[4]!=0;
    $hostgeno{$tmp[1]}{base}{C}=$tmp[5] if $tmp[5]!=0;
    $hostgeno{$tmp[1]}{base}{G}=$tmp[6] if $tmp[6]!=0;
    $hostgeno{$tmp[1]}{base}{T}=$tmp[7] if $tmp[7]!=0;
}
close $hatcg;
open my $atcg,"tabix $ARGV[6] $chrnum|" or die ("Can't open $ARGV[6]\n");
my $sample=`basename $ARGV[6] .acgt.gz|cut -d "." -f 2`;
chomp $sample;
open my $out,"|bgzip >T.$sample.$chrnum.baseCN.gz";
print "Tumor atcg\n";
my $cellularity=$ARGV[7];
while (<$atcg>) {
    my @tmp=split /\t/;
    if (exists $cnv{$tmp[1]}) {
        my $total=$tmp[4]+$tmp[5]+$tmp[6]+$tmp[7];
        if ($total==0){
            print $out $chrnum,"\t",$tmp[1],"\t",$tmp[2],"\t0\t0\t0\t0\t0\n";
            next;
        }elsif ($total > 0 and $total < 2) {
            $amb{$tmp[1]}="tumor total raw reads num less two";
            next;
        }
        next if (exists $amb{$tmp[1]});
        unless (exists $hostgeno{$tmp[1]}) {
            unless (exists $hostdel{$tmp[1]} and $hostdel{$tmp[1]}{c}==2) {
	$amb{$tmp[1]}="host didn't cover";
	next;
            }
        }
        if (exists $tumordel{$tmp[1]}) {
            $total+=$tumordel{$tmp[1]}{c};
        }
        my $hostcnv;
        if (exists $hostcnv{$tmp[1]}) {
            $hostcnv=$hostcnv{$tmp[1]};
        }else {
            $hostcnv=2;
        }
        my $hostRN=$total*$hostcnv*(1-$cellularity)/($hostcnv*(1-$cellularity)+$cnv{$tmp[1]}*$cellularity);
        my %tumor_rawRN;
        $tumor_rawRN{A}=$tmp[4];
        $tumor_rawRN{C}=$tmp[5];
        $tumor_rawRN{G}=$tmp[6];
        $tumor_rawRN{T}=$tmp[7];
        my $typecount=0;
        my $ExistAltAlleleInHost=0;
        if (exists $hostdel{$tmp[1]} and $hostdel{$tmp[1]}{c}==2) {
            $tumordel{$tmp[1]}{c}=$tumordel{$tmp[1]}{c}-$hostRN if (exists $tumordel{$tmp[1]});
        }else {
            for my $allele (keys $hostgeno{$tmp[1]}{base}) {
	if ($tumor_rawRN{$allele}==0) {
	    if ($hostgeno{$tmp[1]}{base}{$allele}/$hostgeno{$tmp[1]}{sum} < 0.1) {
	        delete $hostgeno{$tmp[1]}{base}{$allele};
	    }else {
	        $typecount++;
	    }
	}else {
	    $typecount++;
	}
            }
            my @hostallels=keys $hostgeno{$tmp[1]}{base};
            $ExistAltAlleleInHost=1 unless (scalar(@hostallels)==1 and $hostallels[0] eq $tmp[2]);
            if (exists $hostdel{$tmp[1]} and $hostdel{$tmp[1]}{c}==1){
	$typecount*=2;
	$tumordel{$tmp[1]}{c}=$tumordel{$tmp[1]}{c}-($hostRN/2) if (exists $tumordel{$tmp[1]});
            }
        }
        my $sum;
        my $ExistAlleleNotInHost=0;
        for my $allele (sort keys %tumor_rawRN) {
            $tumor_rawRN{$allele}-=($hostRN/$typecount) if ($typecount!=0 and exists $hostgeno{$tmp[1]}{base}{$allele});
            $ExistAlleleNotInHost=1 if ($tumor_rawRN{$allele}>0 and not exists $hostgeno{$tmp[1]}{base}{$allele});
            $tumor_rawRN{$allele}=0 if $tumor_rawRN{$allele}<0;
            $sum+=$tumor_rawRN{$allele};
        }
        my $realcnv=$cnv{$tmp[1]};
        if ($sum==0) {
            $realcnv=0;
        }
        if (exists $tumordel{$tmp[1]}) {
            $realcnv=round($realcnv*$sum/($sum+$tumordel{$tmp[1]}{c}));
        }
        next if ($ExistAltAlleleInHost==1 and $ExistAlleleNotInHost==0);
        print $out $chrnum,"\t",$tmp[1],"\t",$tmp[2],"\t";
        if ($sum != 0) {
            for my $allele (sort keys %tumor_rawRN) {
	if ($tumor_rawRN{$allele}/$sum*$realcnv<0.5 and $tumor_rawRN{$allele}/$sum*$realcnv >= 0.3 and $tumor_rawRN{$allele} > 2) {
	    $tumor_rawRN{$allele}=1;
	}elsif ($tumor_rawRN{$allele}/$sum*$realcnv>=0.5 and $tumor_rawRN{$allele} > 2) {
	    $tumor_rawRN{$allele}=round($tumor_rawRN{$allele}/$sum*$realcnv);
	}else {
	    $tumor_rawRN{$allele}=0;
	}
	print $out $tumor_rawRN{$allele},"\t";
            }
        }else {
            print $out "0\t0\t0\t0\t";
        }
        print $out $realcnv,"\n";
    }else {
        print $out $chrnum,"\t",$tmp[1],"\t",$tmp[2],"\t0\t0\t0\t0\t0\n";
    }
}
close $atcg;
close $out;
open my $amb,"|bgzip >T.$sample.$chrnum.amb.gz";
print "amb\n";
for my $pos (sort{$a<=>$b} keys %amb) {
    print $amb $chrnum,"\t",$pos,"\t",$amb{$pos},"\n";
}
close $amb;
