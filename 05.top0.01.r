a<-read.table("03.catresult.pl.GreenlandDogs.out",header=T)
a[a$Likelihood>quantile(a$Likelihood,0.99),]->b
round(b$pos-5000,0)->b$start
round(b$pos+5000,0)->b$end
paste("chr",b$chr,sep="")->b$chr
write.table(b[,c(1,7,8)],"05.top0.01.r.10k.region.bed",sep="\t",row.names=F,col.names=F,quote=F)
