#!/usr/bin/env perl
use warnings;
use strict;
my ($reads,$allel,$sam,@tmp,%atgc,@fas,@z,@y,@x);
%atgc=("AA"=>"A","CC"=>"C","GG"=>"G","TT"=>"T","AG" => "R","GA" => "R","CT" => "Y","TC" => "Y","AC" => "M","CA" => "M","GT" => "K","TG" => "K","GC" => "S","CG" => "S","AT" => "W","TA" => "W");
while($reads=<>)
  {
      chomp $reads;
      last if($reads=~/^#CHROM\s+POS\s+ID\s+REF/);
  }
#$reads=~s/000//g;$reads=~s/FAM//g;$reads=~s/LUP//g;
@tmp=split/\s+/,$reads;
for(my $i=9;$i<=$#tmp;$i++)
  {   
      push @fas,">$tmp[$i]\n";
  }
while($reads=<>)
  {
      chomp $reads;
      @tmp=split/\s+/,$reads;
      @z=split/\,/,$tmp[4];
      for(my $i=9;$i<=$#tmp;$i++)
        {
            $sam=$i-9;$allel=0;
            if($tmp[4]=~/,/){
	if($tmp[$i]=~/\.\/\./){
	    $fas[$sam].="N";
	}else{
	    @y=split/:/,$tmp[$i];
	    @x=split/\//,$y[0];
	    if($x[0]==0){
	        $allel=$tmp[3];
	    }else{$allel=$z[$x[0]-1];}
	    if($x[1]==0){
	        $allel.=$tmp[3];
	    }else{$allel.=$z[$x[1]-1];}
	    $fas[$sam].=$atgc{$allel};
	}
            }else{
	if($tmp[$i]=~/\.\/\./){
	    $fas[$sam].="N";
	}elsif($tmp[$i]=~/1\/1/){
	    $fas[$sam].=$tmp[4];
	}elsif($tmp[$i]=~/0\/0/){
	    $fas[$sam].=$tmp[3];
	}else{
	    $allel=$tmp[3].$tmp[4];
	    $fas[$sam].=$atgc{$allel};
	}
            }
        }
  }
print "$_\n" for @fas;

