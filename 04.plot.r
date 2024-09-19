pdf("04.plot.pdf",width=18,height=6)
par(mar=c(2,4,0.2,0.2))
#par(mfrow=c(3,1))
chr1size = 122678785
chr2size = 85426708
chr3size = 91889043
chr4size = 88276631
chr5size = 88915250
chr6size = 77573801
chr7size = 80974532
chr8size = 74330416
chr9size = 61074082
chr10size =69331447
chr11size =74389097
chr12size =72498081
chr13size =63241923
chr14size =60966679
chr15size =64190966
chr16size =59632846
chr17size =64289059
chr18size =55844845
chr19size =53741614
chr20size =58134056
chr21size =50858623
chr22size =61439934
chr23size =52294480
chr24size =47698779
chr25size =51628933
chr26size =38964690
chr27size =45876710
chr28size =41182112
chr29size =41845238
chr30size =40214260
chr31size =39895921
chr32size =38810281
chr33size =31377067
chr34size =42124431
chr35size =26524999
chr36size =30810995
chr37size =30902991
chr38size =23914537
chrsize =c(chr1size,chr2size,chr3size,chr4size,chr5size,chr6size,chr7size,chr8size,chr9size,chr10size,chr11size,chr12size,chr13size,chr14size,chr15size,chr16size,chr17size,chr18size,chr19size,chr20size,chr21size,chr22size,chr23size,chr24size,chr25size,chr26size,chr27size,chr28size,chr29size,chr30size,chr31size,chr32size,chr33size,chr34size,chr35size,chr36size,chr37size,chr38size)
totalsize=2203764842

cols = c("forestgreen","darkblue","brown2","black","cornflowerblue",
         "forestgreen","darkblue","brown2","black","cornflowerblue",
         "forestgreen","darkblue","brown2","black","cornflowerblue",
         "forestgreen","darkblue","brown2","black","cornflowerblue",
         "forestgreen","darkblue","brown2","black","cornflowerblue",
         "forestgreen","darkblue","brown2","black","cornflowerblue",
         "forestgreen","darkblue","brown2","black","cornflowerblue",
         "forestgreen","darkblue","brown2","black","cornflowerblue")
fst<- read.table("03.catresult.pl.GreenlandDogs.out",header=T)

for(i in 1:38)
{
	chrdata = subset(fst,chr==i)
	if(i==1)
	{
		relative_pos=0
	}
	else
	{
		relative_pos = sum(chrsize[1:(i-1)])
	}
	if(i<38)
	{
		plot((as.numeric(chrdata[[2]])+relative_pos+5000),chrdata[[3]],pch=16,xlim=c(0,totalsize),ylim=c(0,2000),col=cols[i],cex=0.5,cex.axis=1,xaxt="n",yaxt="n",xlab="",ylab="")
		axis(1,at=relative_pos+(chrsize[i]/2),labels=i,line=-0.5,tick=F)
		#segments(1,quantile(fst[[3]], 0.999),2203764842, quantile(fst[[3]], 0.999),lty=3)
	}
	else
	{
		plot((as.numeric(chrdata[[2]])+relative_pos+5000),chrdata[[3]],pch=16,xlim=c(0,totalsize),ylim=c(0,2000),col=cols[i],cex=0.5,cex.axis=1,xaxt="n",xlab="",ylab="Likelihood")
		axis(1,at=relative_pos+(chrsize[i]/2),labels=i,line=-0.5,tick=F)
		segments(1,quantile(fst[[3]], 0.99),2203764842, quantile(fst[[3]], 0.99),lty=3)
	}
	par(new=T)
}
par(new=F)