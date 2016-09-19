####Copy number variation using cn.mops
rm(list=ls()) 

setwd("C:/Users/hji3/Desktop/Thesis Work/Panama Isolate Genome Project")


library(cn.mops)
source("https://bioconductor.org/biocLite.R")
biocLite("DNAcopy")
library(DNAcopy)

#########
#cn.mops
#for merged bams
setwd("C:/Users/hji3/Desktop/Thesis Work/Panama Isolate Genome Project/mergedBAMS")

BAMFiles1 <- list.files(pattern=".bam$")
#for only PE bams
setwd("C:/Users/hji3/Desktop/Thesis Work/Panama Isolate Genome Project/PE_BAMS")


BAMFiles2 <- list.files(pattern=".bam$")


#need refSeqName
refseqnames <- c("Supercontig_1.1","Supercontig_1.2","Supercontig_1.3","Supercontig_1.4","Supercontig_1.5","Supercontig_1.6","Supercontig_1.7","Supercontig_1.8","Supercontig_1.9","Supercontig_1.10","Supercontig_1.11","Supercontig_1.12","Supercontig_1.13","Supercontig_1.14","Supercontig_1.15","Supercontig_1.16","Supercontig_1.17","Supercontig_1.18","Supercontig_1.19","Supercontig_1.20","Supercontig_1.21","Supercontig_1.22","Supercontig_1.23","Supercontig_1.24","Supercontig_1.25","Supercontig_1.26","Supercontig_1.27","Supercontig_1.28","Supercontig_1.29","Supercontig_1.30","Supercontig_1.31","Supercontig_1.32","Supercontig_1.33","Supercontig_1.34","Supercontig_1.35","Supercontig_1.36","Supercontig_1.37","Supercontig_1.38","Supercontig_1.39","Supercontig_1.40","Supercontig_1.41","Supercontig_1.42","Supercontig_1.43","Supercontig_1.44","Supercontig_1.45","Supercontig_1.46","Supercontig_1.47","Supercontig_1.48","Supercontig_1.49","Supercontig_1.50","Supercontig_1.51","Supercontig_1.52","Supercontig_1.53","Supercontig_1.54","Supercontig_1.55","Supercontig_1.56","Supercontig_1.57","Supercontig_1.58","Supercontig_1.59","Supercontig_1.60","Supercontig_1.61","Supercontig_1.62","Supercontig_1.63","Supercontig_1.64","Supercontig_1.65","Supercontig_1.66","Supercontig_1.67","Supercontig_1.68","Supercontig_1.69")
## run to get Granges object
##mode should be paired
setwd("C:/Users/hji3/Desktop/Thesis Work/Panama Isolate Genome Project/mergedBAMS")
bamDataRanges1 <- getReadCountsFromBAM(BAMFiles1, sampleNames=c('Campana','Sora','JEL410','RioMaria','JEL412','JEL413'), refSeqName=refseqnames, mode="paired")

res1 <- cn.mops(bamDataRanges1)
rescount1 <- calcIntegerCopyNumbers(res1)

##for only PE Bams
##
setwd("C:/Users/hji3/Desktop/Thesis Work/Panama Isolate Genome Project/PE_BAMS")
bamDataRanges2 <- getReadCountsFromBAM(BAMFiles2, sampleNames=c('Campana','Sora','JEL410','RioMaria','JEL412','JEL413'), refSeqName=refseqnames, mode="paired")

res2 <- cn.mops(bamDataRanges2)
rescount2 <- calcIntegerCopyNumbers(res2)

##CNV per individual are stored in cnvs


################
##JEL410


#get index for JEL410
JEL410idx <- which(cnvs(rescount1)$sampleName=="JEL410")
#pull out only JEL410
JEL410CN <- (cnvs(rescount1)[JEL410idx])
##make a matrix with supercontig number and CN count
JEL410M <- matrix(nrow=length(JEL410idx), ncol=2)
#puts in the contig name
for(i in 1:length(JEL410idx)){
  JEL410M[i,1] <- as.character(JEL410CN@seqnames[i]@values[1])}
#puts in the CN number
JEL410M[,2] <- JEL410CN$CN

################
##JEL412

#get index for JEL412
JEL412idx <- which(cnvs(rescount1)$sampleName=="JEL412")
#pull out only JEL412
JEL412CN <- (cnvs(rescount1)[JEL412idx])
##make a matrix with supercontig number and CN count
JEL412M <- matrix(nrow=length(JEL412idx), ncol=2)
#puts in the contig name
for(i in 1:length(JEL412idx)){
  JEL412M[i,1] <- as.character(JEL412CN@seqnames[i]@values[1])}
#puts in the CN number
JEL412M[,2] <- JEL412CN$CN

################
##JEL413

#get index for JEL413
JEL413idx <- which(cnvs(rescount1)$sampleName=="JEL413")
#pull out only JEL413
JEL413CN <- (cnvs(rescount1)[JEL413idx])
##make a matrix with supercontig number and CN count
JEL413M <- matrix(nrow=length(JEL413idx), ncol=2)
#puts in the contig name
for(i in 1:length(JEL413idx)){
  JEL413M[i,1] <- as.character(JEL413CN@seqnames[i]@values[1])}
#puts in the CN number
JEL413M[,2] <- JEL413CN$CN

################
##Campana

#get index for Campana
Campanaidx <- which(cnvs(rescount1)$sampleName=="Campana")
#pull out only Campana
CampanaCN <- (cnvs(rescount1)[Campanaidx])
##make a matrix with supercontig number and CN count
CampanaM <- matrix(nrow=length(Campanaidx), ncol=2)
#puts in the contig name
for(i in 1:length(Campanaidx)){
  CampanaM[i,1] <- as.character(CampanaCN@seqnames[i]@values[1])}
#puts in the CN number
CampanaM[,2] <- CampanaCN$CN

################
##RioMaria

#get index for RioMaria
RioMariaidx <- which(cnvs(rescount1)$sampleName=="RioMaria")
#pull out only RioMaria
RioMariaCN <- (cnvs(rescount1)[RioMariaidx])
##make a matrix with supercontig number and CN count
RioMariaM <- matrix(nrow=length(RioMariaidx), ncol=2)
#puts in the contig name
for(i in 1:length(RioMariaidx)){
  RioMariaM[i,1] <- as.character(RioMariaCN@seqnames[i]@values[1])}
#puts in the CN number
RioMariaM[,2] <- RioMariaCN$CN

################
##Sora

#get index for Sora
Soraidx <- which(cnvs(rescount1)$sampleName=="Sora")
#pull out only Sora
SoraCN <- (cnvs(rescount1)[Soraidx])
##make a matrix with supercontig number and CN count
SoraM <- matrix(nrow=length(Soraidx), ncol=2)
#puts in the contig name
for(i in 1:length(Soraidx)){
  SoraM[i,1] <- as.character(SoraCN@seqnames[i]@values[1])}
#puts in the CN number
SoraM[,2] <- SoraCN$CN

setwd("C:/Users/hji3/Desktop/Thesis Work/Panama Isolate Genome Project/CN count from cn.mops")

dir.create("CN count from cn.mops")
write.table(CampanaM,file="CampanaRC.tab",sep="\t",quote=F,row.names=F)
write.table(SoraM,file="SoraRC.tab",sep="\t",quote=F,row.names=F)
write.table(RioMariaM,file="RioMariaRC.tab",sep="\t",quote=F,row.names=F)
write.table(JEL410M,file="JEL410RC.tab",sep="\t",quote=F,row.names=F)
write.table(JEL412M,file="JEL412RC.tab",sep="\t",quote=F,row.names=F)
write.table(JEL413M,file="JEL413RC.tab",sep="\t",quote=F,row.names=F)


##CNV regions are stored in cnvr
regions <- cnvr(rescount1)
regions1 <- as.data.frame(regions)
write.table(regions1,file="cnvr_regions.tab",sep="\t",quote=F,row.names=F)


plot(res1,which=1)
plot(res2,which=1)

##plots chromosome plot for all samples for only chrom1

segplot(rescount1,seqnames='Supercontig_1.1', sampleIdx=5, ylim=c(-3,3))
segplot(rescount1,seqnames='Supercontig_1.2', sampleIdx=1)
segplot(rescount1,seqnames='Supercontig_1.3', sampleIdx=1)
segplot(rescount1,seqnames='Supercontig_1.4', sampleIdx=1)
segplot(rescount1,seqnames='Supercontig_1.5', sampleIdx=1)
segplot(rescount1,seqnames='Supercontig_1.6', sampleIdx=1)
segplot(rescount1,seqnames='Supercontig_1.7', sampleIdx=1)
segplot(rescount1,seqnames='Supercontig_1.8', sampleIdx=1)
segplot(rescount1,seqnames='Supercontig_1.9', sampleIdx=1)
segplot(rescount1,seqnames='Supercontig_1.10', sampleIdx=1)
segplot(rescount1,seqnames='Supercontig_1.11', sampleIdx=1)
segplot(rescount1,seqnames='Supercontig_1.12', sampleIdx=1)
segplot(rescount1,seqnames='Supercontig_1.13', sampleIdx=1)











segplot(rescount1,seqnames='Supercontig_1.2', sampleIdx=2)
segplot(rescount2,seqnames='Supercontig_1.2', sampleIdx=2)

segplot(rescount,seqnames=refseqnames[1:10], sampleIdx=1)
segplot(rescount,cbys.layout=layout)

