setwd("~/Google Drive/Experiments/yose_popgen/fluidigm_data/")
##############
##### SPLIT FASTQ FILE INTO SEPARATE SAMPLE FASTQ FILES
# read in input files      
suppressPackageStartupMessages(library("ShortRead"))
fq <- readFastq((("mergeOnly/yoseMin25/flash2.extendedFrags.fastq.gz")))
#fq <- readFastq(file.path(output,paste0(flash_prefix,".extendedFrags.fastq.gz")))
nms <- as.character(id(fq))
id <- sapply(strsplit(sapply(strsplit(nms,split=" "),"[[",2L),split=":"),"[[",4L)
split_tt <- split(fq,id)
dir.create("split_samples")
procs = 2
mclapply(names(split_tt), function(x){
  writeFastq(split_tt[[x]],file.path(paste("split_samples",sep="."),paste("Sample",x,"fastq.gz",sep=".")))
}, mc.cores = procs)
#######
# parse reference file - fix names
getwd() 
ref = readDNAStringSet("testing/consensus.Amplicons.refListEdit.fa")
nm = names(ref)
nm1 = sapply(nm, function(x) unlist(strsplit(x, split="\\|"))[1])
nm2 = sapply(nm1, function(x) unlist(strsplit(x, split=":"))[1])
nm1 = sapply(nm1, function(x) unlist(strsplit(x, split=":"))[2])
nm3 = paste(nm1,"REF",nm2,sep = "|")
names(ref) = nm3
writeXStringSet(ref, "testing/consensus.Amplicons.refListEdit2.fa")
##############################################################
##############################################################
##############################################################
##############################################################
# DO ALIGNMENTS
# system("cd /clusterfs/vector/scratch/tpoorten/yose/align/denovoTargetCapturePopGen/")
dir.create("align")
setwd("~/Google Drive/Experiments/yose_popgen/fluidigm_data/align")
dir.create("bams")
fqFiles = list.files("../split_samples/")
threads = 2
i=1
if(!("refList.fa.bwt" %in% list.files())){
  print("Indexing Reference for bwa")
  system("bwa index refList.fa")
}
if(!("refList.fa.fai" %in% list.files())){
  print("Indexing Reference for samtools")
  system("samtools faidx refList.fa")
}
for(i in 1:length(fqFiles)){
  # for(i in 5:10){
  print(i)
  fq = fqFiles[i]
  fqBase = sub(".fastq.gz","",fq)
  print(c(fqBase))
  print("Aligning reads")
  system(paste0("bwa mem -t ", threads," refList.fa ../split_samples/",fq," > bams/",fqBase,".sam"))
  print("Convert to Bam")
  system(paste0("samtools view -bt refList.fa.fai -o bams/", fqBase,".bam bams/",fqBase,".sam"))
  print("Sorting reads")
  system(paste0("samtools sort bams/",fqBase,".bam bams/",fqBase,".sort"))
  print("Clean up")
  system(paste0("rm bams/*.sam")); 
  #   system(paste0("rm bams/*.bam")) this will delete the final bam file too ..? ##*/
}
### clean up anything in bams folder that doesn't end with sort.bam - BE CAREFUL !!!!!!
cd bams
ls | grep -v sort.bam
ls | grep -v sort.bam | xargs rm
#######################################################
# Evaluate bam files - Coverage
module load bedtools
mkdir aln_eval
#
setwd("/clusterfs/vector/scratch/tpoorten/yose/align/denovoTargetCapturePopGen/")
rm(list=ls())
bamFiles = list.files("bams/")
bamFiles = bamFiles[grep("sort.bam",bamFiles)]
i=1
for(i in 1:length(bamFiles)){
  print(i)
  bm = bamFiles[i]
  bmBase = sub(".sort.bam","",bm)
  system(paste0("bedtools genomecov -ibam bams/",bm," -d -split > aln_eval/",bmBase,".cov.txt"))
  #   system(paste0("bedtools genomecov -ibam bams/",bm," -bg > aln_eval/",bmBase,".cov.txt"))
  #   system(paste0("bedtools merge -i aln_eval/",bmBase,".cov.txt -d 10 > aln_eval/",bmBase,".covMerge.txt "))  
}
####
# make Table of results
setwd("/clusterfs/vector/scratch/tpoorten/yose/align/denovoTargetCapturePopGen/")
rm(list=ls())
covFiles = list.files("aln_eval/")
i=1
ref = system("grep '>' refList.fa ", intern=T)
ref = sub(">","",ref)
tmp = rep(NA,length(ref))
resMed = data.frame(empty = tmp, stringsAsFactors = F, row.names=ref)
resLen = data.frame(empty = tmp, stringsAsFactors = F, row.names=ref)
for(i in 1:length(covFiles)){
  print(i)
  cover = read.table(paste0("aln_eval/",covFiles[i]),sep="\t", stringsAsFactors=F)
  covName = sub(".cov.txt","",covFiles[i])
  coverLen = tapply(cover[,1],cover[,1],length)
  coverMedian = tapply(cover[,3],cover[,1],median)  
  if(i==1){
    #     resLen = coverLen[match(ref,names(coverLen))]
    #     resLen
    #     resMed = coverMedian[match(ref,names(coverMedian))]
    resMed[,i] = coverMedian[match(ref,names(coverMedian))]
    colnames(resMed)[i] = covName
    resLen[,i] = coverLen[match(ref,names(coverLen))]
    colnames(resLen)[i] = covName
  } else{
    resMed[,i] = coverMedian[match(ref,names(coverMedian))] 
    colnames(resMed)[i] = covName
    resLen[,i] = coverLen[match(ref,names(coverLen))]
    colnames(resLen)[i] = covName
  }
}
write.table(resMed, "bwamem_medianCoverage.txt", sep="\t", quote=F)
write.table(resLen, "bwamem_lengthCoverage.txt", sep="\t", quote=F)
#
xx = read.table("bwamem_medianCoverage.txt",sep="\t", header=T, stringsAsFactors=F)
max(xx, na.rm = T)
# 2072
hist(as.matrix(xx))
#


########################################################################
# ADD READ GROUPS IF NOT DONE ALREADY - using bamaddrg program
setwd("/clusterfs/vector/scratch/tpoorten/yose/align/")
rm(list=ls())
bamFiles = list.files("bams/")
bamFiles = bamFiles[grep("sort.bam$",bamFiles,perl=T)]
i=1
dir.create("bamsRG")
for(i in 1:length(bamFiles)){
  print(i)
  bam = bamFiles[i]
  bamB = sub("Sample.","",bamFiles[i])
  bamB = sub(".sort.bam","",bamB)
#   system(paste0("$BAMADDRG -b bams/",bam[i]," -s ",bamB," > bamsRG/Sample.",bamB,".sortRG.bam"))
  system(paste0("/global/home/users/tpoorten/bin/bamaddrg/bamaddrg -b bams/",bam," -s ",bamB," > bamsRG/Sample.",bamB,".sortRG.bam"))
}
#
#
#/global/home/users/tpoorten/bin/bamaddrg/bamaddrg -b bams/ -s  > bamsRG/Sample.
#
########################################################################
#  Do the indel realignments
#### do this after adding read group tags, then link bam files to gatk/bams/, and indexing the bam files
module load java
R
setwd("/clusterfs/vector/scratch/tpoorten/yose/align/gatk")
rm(list=ls())
GATK = "/global/home/users/tpoorten/bin/GenomeAnalysisTK.jar"
#
bamFiles = list.files("bams/")
bamFiles = bamFiles[grep("sortRG.bam$",bamFiles,perl=T)]
i=1
#
dir.create("bamsReAlign")
dir.create("intervals")

#for(i in 2:3){
for(i in 1:length(bamFiles)){
  print(i);
  bam = bamFiles[i]
  bamB = sub("Sample.","",bamFiles[i])
  bamB = sub(".sortRG.bam","",bamB)
  print(bamB)
  system(paste0("java -Xmx1g -jar ",GATK," -T RealignerTargetCreator -R refTargets.fa -o intervals/",bamB,".intervals -I bams/",bam))
  system(paste0("java -Xmx1g -jar ",GATK," -T IndelRealigner -R refTargets.fa -I bams/",bam," -targetIntervals intervals/",bamB,".intervals -o bamsReAlign/",bam))
}
#######################
# get active region list
library(Biostrings)
xx = readDNAStringSet("refTargets.fa")
write.table(names(xx),"activeRegions.txt",quote=F,row.names=F,col.names=F)
write.table(paste(names(xx),paste(1,width(xx),sep="-"),sep=":"),"activeRegions2.txt",quote=F,row.names=F,col.names=F)
write.table(paste(names(xx),paste(1,width(xx),sep=" "),sep=" "),"activeRegions.bed",quote=F,row.names=F,col.names=F)
###########################################################################
###############

#################################################################################
#################################################################################
#################################################################################
#  Run HAPLOTYPE CALLER v2
# add: define active regions; set min read length
cd /clusterfs/vector/scratch/tpoorten/yose/align/gatk
module load java
R
#
setwd("/clusterfs/vector/scratch/tpoorten/yose/align/gatk")
rm(list=ls())
GATK = "/global/home/users/tpoorten/bin/GenomeAnalysisTK.jar"
bamFiles = list.files("bamsReAlign/")
bamFiles = bamFiles[grep("sortRG.bam$",bamFiles,perl=T)]
i=2
dir.create("logs2")
dir.create("vcfs2")
#for(i in 1:10){
#for(i in 11:30){
#for(i in 31:50){
#for(i in 51:100){
#for(i in 101:150){
#for(i in 151:189){
#for(i in 1:80){
#for(i in 81:150){
#for(i in 61:70){
#for(i in 79:80){
  print(i);
  bam = bamFiles[i]
  bamB = sub("Sample.","",bamFiles[i])
  bamB = sub(".sortRG.bam","",bamB)
  print(bamB)
#  system(paste0("java -Xmx10g -jar ",GATK," -T HaplotypeCaller -R refTargets.fa -I bamsReAlign/",bam," --emitRefConfidence GVCF -o vcfs2/",bamB,".g.vcf -AR activeRegions.list --pcr_indel_model CONSERVATIVE -dfrac 1 -mmq 50 -rf NotPrimaryAlignment -rf BadCigar -rf ReadLength -minRead 250 -maxRead 600 -log logs2/",bamB,".log"))
#  system(paste0("java -Xmx25g -jar ",GATK," -T HaplotypeCaller -nct 8 -R refTargets.fa -I bamsReAlign/",bam," --emitRefConfidence GVCF -o vcfs2/",bamB,".g.vcf -AR activeRegions.list --pcr_indel_model CONSERVATIVE -dfrac 1 -mmq 50 -rf NotPrimaryAlignment -rf BadCigar -rf ReadLength -minRead 250 -maxRead 600 -log logs2/",bamB,".log"))

#	system(paste0("java -Xmx25g -jar ",GATK," -T HaplotypeCaller -R refTargets.fa -I bamsReAlign/",bam," --emitRefConfidence GVCF -o vcfs2/",bamB,".g.vcf -AR activeRegions.bed --pcr_indel_model CONSERVATIVE -dfrac 1 -mmq 50 -rf NotPrimaryAlignment -rf BadCigar -rf ReadLength -minRead 250 -maxRead 600 -L RS1_REF_RKS7981:1-299 --activeRegionMaxSize 600 --activeRegionOut AR.out --dontTrimActiveRegions -maxNumHaplotypesInPopulation 2 --bamOutput haplotypes.bam --bamWriterType CALLED_HAPLOTYPES -log logs2/",bamB,".log"))
	system(paste0("java -Xmx25g -jar ",GATK," -T HaplotypeCaller -R refTargets.fa -I bamsReAlign/",bam," --emitRefConfidence GVCF -o vcfs2/",bamB,".g.vcf -AR activeRegions.bed --pcr_indel_model CONSERVATIVE -dfrac 1 -mmq 50 -rf NotPrimaryAlignment -rf BadCigar -rf ReadLength -minRead 250 -maxRead 600 -L RS1_REF_RKS7981:1-299 --activeRegionMaxSize 600 --activeRegionOut AR.out --dontTrimActiveRegions -maxNumHaplotypesInPopulation 2 --bamOutput haplotypes.bam --bamWriterType CALLED_HAPLOTYPES -log logs2/",bamB,".log"))
}
-AR activeRegions.txt
-minRead 250
--activeRegionOut
-maxNumHaplotypesinPopulation 2 --bamOutput haplotypes.bam --bamWriterType CALLED_HAPLOTYPES


#################################################################################
#################################################################################
#  Read backed phasing
cd /clusterfs/vector/scratch/tpoorten/yose/align/gatk
module load java
R
#
setwd("/clusterfs/vector/scratch/tpoorten/yose/align/gatk")
rm(list=ls())
GATK = "/global/home/users/tpoorten/bin/GenomeAnalysisTK.jar"
samples = readLines("filteredGenoSNP_v2_HQsamples.txt")
i=1
dir.create("phased")
for(i in 5:length(samples)){
	system(paste0("java -jar ",GATK," -T ReadBackedPhasing -R refTargets.fa --variant filteredGenoSNP_v2_HQsamples.vcf -o phased/phased3_SNPs_",samples[i],".vcf --phaseQualityThresh 20.0 -dfrac 1 -mmq 50 -rf NotPrimaryAlignment -rf BadCigar --sampleToPhase ",samples[i]," -I bamsReAlign/Sample.",samples[i],".sortRG.bam"))
}
#########################################################################
###########################################################################
