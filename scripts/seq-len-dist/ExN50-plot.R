#!/usr/bin/RScript


# source:
#https://github.com/yuragal/fradians-rnaseq/blob/e0dc282dbf16956468ee5376c519ed199e729ce2/2.7.annotation.R

#############################################
#Figure 2b
#############################################
#ExN50 plot
#See https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Contig-Nx-and-ExN50-stats#contig-ex90n50-statistic-and-ex90-transcript-count

dir="/Volumes/Mac_HD2/proj_dir/docs/"
exn50=read.table(paste0(dir,'ExN50.stats'),sep='\t',header=TRUE)
gl=seq(0,43,5)
pdf('ExN50.pdf',width=6,height=5)
par(mar=c(5,4,4,5)+.1)
plot(exn50$E,exn50$ExN50,type='l',col='red',xlab="Expression percentile",ylab="N50, bp")
par(new=TRUE)
plot(exn50$E,exn50$num_transcripts,type='l',col='blue',xaxt="n",yaxt="n",xlab="",ylab="")
#axis(4,labels=paste0(gl,'K'),at=gl*1000)
axis(4, xaxp = c(2, 9, 7))


mtext("Number of transcripts",side=4,line=2)
legend("topleft",col=c("red","blue"),lty=1,legend=c("ExN50","Transcripts"))
abline(v=97, lty=3)
dev.off()
#axis(1,labels=c(96),at=c(96))


