setwd("~/Google Drive/Experiments/bd_penny/reCall0815/")
setwd("~/Google Drive/Experiments/bd_penny/reCall0815_newRef/")
library("VariantAnnotation")
rm(list=ls())
# vcf = readVcf("rawGeno.vcf","Bdpenny")
# vcf = readVcf("filter.selectGenoMIX.vcf","Bdpenny")
# vcf = readVcf("rawGenoSNP.UCL.vcf","Bdpenny")
# vcf = readVcf("rawGenoSNP.Cal.vcf","Bdpenny")
# #
# vcf = readVcf("rawGenoINDEL.UCL.vcf","Bdpenny")
# vcf = readVcf("rawGenoINDEL.Cal.vcf","Bdpenny")
#
vcf = readVcf("filter.selectGenoSNP.Cal.vcf","Bdpenny")
vcf = readVcf("filter.selectGenoINDEL.Cal.vcf","Bdpenny")

vcf = readVcf("filter.selectGenoSNP.Calmmq60.vcf","Bdpenny")
vcf = readVcf("filter.selectGenoINDEL.Calmmq60.vcf","Bdpenny")

vcf = readVcf("filter.selectGenoSNP.UCL.vcf","Bdpenny")
vcf = readVcf("filter.selectGenoINDEL.UCL.vcf","Bdpenny")

#writeVcf(vcf,"output/rawMultiAllel_subset2.vcf")
####
# mat <- genotypeToSnpMatrix(vcf)
# mat = mat$genotypes
# mat.summ = col.summary(mat)
# hist(mat.summ$P.AB, ylim=c(0,100))
# # cut out contigs that are likely paralogs
contigs = as.character(seqnames(vcf))
# head(contigs)
contigs1 = sapply(contigs, function(x) unlist(strsplit(x, split=":"))[1])
head(contigs1)
####
# allele count - get rid of snps with only one allele
# dim(vcf)
# maxAllelecount = lapply(info(vcf)$AC, max)
# vcf = vcf[which(maxAllelecount == 1),]
# dim(vcf)
####
vcf.geno = geno(vcf)$GT
head(vcf.geno); dim(vcf.geno)
###########
# turn multi-allelic calls into "." missing data <- THIS IS NEW 7/19/15
# apply(vcf.geno, 2, table)
# vcf.geno[which(vcf.geno != "0/0" & vcf.geno != "0/1" & vcf.geno != "1/1" & vcf.geno != ".")] = "."

# turn low GQ into "." missing data <- THIS IS NEW 7/19/15
vcf.gq = geno(vcf)$GQ
hist(vcf.gq, xlim=c(0,100), breaks=20)
apply(vcf.gq,2,function(x) table(x<99))
# table(vcf.gq)
table(vcf.geno)
vcf.geno[which(vcf.gq < 20)] = "."
vcf.geno[which(vcf.gq < 40)] = "."

geno(vcf)$GT = vcf.geno
##
vcf.dp = geno(vcf)$DP
apply(vcf.dp, 2, function(x) mean(x, na.rm = T))
apply(vcf.dp, 2, function(x) median(x, na.rm = T))
apply(vcf.dp, 2, function(x) sd(x, na.rm = T))*5
# cutoff of max DP
apply(vcf.dp, 2, function(x) median(x, na.rm = T)) + apply(vcf.dp, 2, function(x) sd(x, na.rm = T))*6
# apply(vcf.dp, 2, table)
hist(vcf.dp,xlim=c(0,30),breaks=200)
summary(vcf.dp)
vcf.geno[which(vcf.dp < 4)] = "."
geno(vcf)$GT = vcf.geno
table(geno(vcf)$GT)
############
# count polymorphic sites between samples (ignore ref)
head(apply(vcf.geno,1,function(x) x[1] != x[2]))
table(apply(vcf.geno,1,function(x) x[1] != x[2]))
table(apply(vcf.geno,1,function(x) x[1] != x[2] & x[1] != "." & x[2] != "."))
############
#
vcf.rowData = rowData(vcf)
vcf.info = info(vcf)
vcf.infoHeader = info(header(vcf))
vcf.genoHeader = geno(header(vcf))
#
as.matrix(vcf.infoHeader)
head(vcf.info)
contigs = rownames(vcf.info)
head(contigs)
contigs1 = sapply(contigs, function(x) unlist(strsplit(x, split=":"))[1])
table(contigs1)
vcf.geno[541,]
vcf.info[541,]
#
as.matrix(vcf.infoHeader)
hist(vcf.info$DP, breaks =20)
median(vcf.info$DP) + (sd(vcf.info$DP)*6)
vcf.rowData[which(vcf.info$DP > 125000),]
hist(vcf.info$DP / vcf.info$AN)
xx=vcf.rowData[which(vcf.info$DP / vcf.info$AN > 40000),]
as.character(seqnames(xx))

hist(vcf.info$BaseQRankSum, breaks =20) # < -5
hist(vcf.info$ClippingRankSum, breaks =20)
hist(vcf.info$FS, breaks =20) # nothing (maxlenRef: > 5)
hist(vcf.info$FS, breaks =200,xlim=c(0,100)) # nothing (maxlenRef: > 5)
hist(vcf.info$MQ, breaks=20) # nothing (maxlenRef: < 50)
table(vcf.info$MQ<60)
# hist(vcf.info$GQ, breaks=20) #
#
as.matrix(vcf.infoHeader)
hist(vcf.info$InbreedingCoeff, breaks=20) # 
hist(vcf.info$MQRankSum, breaks=200) # < -3 or -2 (maxlenRef: < -5)
hist(vcf.info$MQRankSum, breaks=200,xlim=c(-5,5)) # < -3 or -2 (maxlenRef: < -5)
hist(vcf.info$QD, breaks=40) # MAYBE TRY A RANGE OF VALUES HERE - FROM 2 TO 10
summary(vcf.info$QD)
hist(vcf.info$DP, breaks=400) # 
hist(vcf.info$DP, breaks=400,xlim=c(0,100)) # 
hist(vcf.info$ReadPosRankSum, breaks=20) # < -8
hist(vcf.info$ReadPosRankSum, breaks=200,xlim=c(-10,10)) # < -8
hist(vcf.info$AN, breaks =50) # < 100 or 200
hist(vcf.info$SOR, breaks =50) # 
#
hist(vcf.rowData$QUAL, breaks=1000, xlim=c(0,400))
summary(vcf.rowData$QUAL)

########
filter1 = which(vcf.info$DP > 80)
filter1 = which(vcf.info$DP > 150)
filter2 = which(vcf.info$QD < 8)
filter3 = which(vcf.info$ReadPosRankSum < -8)
filters = unique(c(filter1, filter2, filter3))#, filter3, filter4, filter5, filter6, filter7, filter8, filter9))
length(filters)
# new object of filtered sites
dim(vcf.geno)
vcf = vcf[-filters,]
vcf.geno = geno(vcf)$GT
head(vcf.geno); dim(vcf.geno)
##################################################
# get median coverage per chromosome
head(vcf.dp)
#
chrom.cov = apply(vcf.dp, 2, function(x) tapply(x, contigs1, function(x) median(x, na.rm = T)))
chrom.covMean = apply(chrom.cov,2,function(x) median(x, na.rm=T))
chrom.covNorm = t(apply(chrom.cov, 1, function(x) x / chrom.covMean))
#
head(contigs1)
length(contigs1)
dim(chrom.cov)
#chrom.covMeanChrom1 = apply(chrom.cov[which(contigs1 == "bden_JEL423_supercont1.1"),],2,function(x) mean(x, na.rm=T))
chrom.covNormByChrom1 = t(apply(chrom.cov, 1, function(x) x / chrom.cov[2,]))
# RESULTS - THIS GENERALLY CONFIRMS CCNV IN FARRER ET AL
#   ACON is putatively triploid across the largest six supercontigs, 
#   whereas CON2A has lost a copy of supercontig IV and gained a copy of supercontigs V. 
#   APEP has gained a copy of supercontigs V.
#
chrom.covNorm = as.data.frame(chrom.covNorm)
#CON2A
cbind(rownames(chrom.covNorm), fold= round((chrom.covNorm$CON2A*3/2) - (chrom.covNorm$ACON*3/2),digits = 1))
#APEP
cbind(rownames(chrom.covNorm), fold= round((chrom.covNorm$APEP*3/2) - (chrom.covNorm$ACON*3/2),digits = 1))
############
acon.af = geno(vcf)
acon.af = acon.af$AD
# choose strain
j = 3
acon.af = acon.af[which(vcf.geno[,j] == "0/1" | vcf.geno[,j] == "0/2" | vcf.geno[,j] == "1/2"),]
acon.af = acon.af[,j]
head(acon.af.contigs)
acon.af = lapply(acon.af, function(x) max(x / sum(x)))
acon.af.contigs = sapply(names(acon.af), function(x) unlist(strsplit(x, split=":"))[1])
hist(unlist(acon.af[which(acon.af.contigs == "bden_JEL423_supercont1.1")]), breaks=100)
for(i in 1:20){
  hist(unlist(acon.af[which(acon.af.contigs == paste0("bden_JEL423_supercont1.",i))]), breaks=40, xlab=i, main=i)
  cat ("Press [enter] to continue")
  line <- readline()
}
1/5
###################################################
# subset polymorphic sites between samples (ignore ref)
head(apply(vcf.geno,1,function(x) x[1] != x[2]))
table(apply(vcf.geno,1,function(x) x[1] != x[2]))
table(apply(vcf.geno,1,function(x) x[1] != x[2] & x[1] != "." & x[2] != "."))
#
dim(vcf)
vcf = vcf[which(apply(vcf.geno,1,function(x) x[1] != x[2] & x[1] != "." & x[2] != ".")),]
dim(vcf)
vcf = vcf[,c(2,1)]
head(geno(vcf)$GT)
#
writeVcf(vcf,"filter.selectGenoINDEL.diffs.vcf")
writeVcf(vcf,"filter.selectGenoSNP.Cal.diffs.vcf")
writeVcf(vcf,"filter2.selectGenoSNP.Cal.diffs.vcf")
writeVcf(vcf,"filter3.selectGenoSNP.Calmmq60.diffs.vcf")

# extra filter1: filter out sample-wise: GQ < 20; DP < 4; dataset: DP < 150; 
writeVcf(vcf,"filter1.selectGenoSNP.Calmmq60.diffs.vcf")
writeVcf(vcf,"filter2.selectGenoINDEL.Calmmq60.diffs.vcf")
###########
# subset for UCL
#sites with same genotype in all
invariant1 = which(apply(vcf.geno,1,function(x) x[1] == x[2] & x[2] == x[3]))
length(invariant1)
#sites with 2 missing data
invariant2 = which(apply(vcf.geno,1,function(x) length(which(x=="."))) == 2)
invar = unique(c(invariant1, invariant2))
vcf = vcf[-invar,]
vcf.geno=geno(vcf)$GT
dim(vcf)
# cut sites with 1 missing data AND invariant genotype calls
#missing1 = which(apply(vcf.geno,1,function(x) length(which(x=="."))) == 1)
missing1 = which(apply(vcf.geno,1,function(x) length(which(x==".")) == 1 & length(unique(x)) < 3 ))
vcf = vcf[-missing1,]
dim(vcf)
##
#writeVcf(vcf,"filter.selectGenoINDEL.diffs.vcf")
writeVcf(vcf,"filter.selectGenoSNP.UCL.diffs.vcf")
#############
# subset for UCL - separately for each isolate
#sites with same genotype in all
colnames(vcf)
apep = vcf[,c(1,2)]
apep.geno = geno(apep)$GT
dim(apep)
apep = apep[which(apply(apep.geno,1,function(x) x[1] != x[2] & x[1] != "." & x[2] != ".")),]
dim(apep)
head(geno(apep)$GT)
# writeVcf(apep,"filter.selectGenoINDEL.diffs.vcf")
writeVcf(apep,"filter.selectGenoSNP.UCL_APEP.diffs.vcf")
#####
con2a = vcf[,c(1,3)]
con2a.geno = geno(con2a)$GT
dim(con2a)
con2a = con2a[which(apply(con2a.geno,1,function(x) x[1] != x[2] & x[1] != "." & x[2] != ".")),]
dim(con2a)
head(geno(con2a)$GT)
# writeVcf(con2a,"filter.selectGenoINDEL.diffs.vcf")
writeVcf(con2a,"filter.selectGenoSNP.UCL_CON2A.diffs.vcf")
#
#
#
###############################################################################
###############################################################################
#  run bedtools annotate
#bedtools annotate -counts -files filter.selectGenoINDEL.diffs.vcf -i batrachochytrium_dendrobatidis_1_transcripts.CDS.gtf > indel.CDS.countsGTF.txt

bedtools annotate -counts -files filter1.selectGenoINDEL.Calmmq60.diffs.vcf -i batrachochytrium_dendrobatidis_1_transcripts_CDS.gff3 > indel.Calmmq60.CDS.countsGFF3.txt
#### indels, newRef
# filter 1: DP > 150; filter0: nothing; filter 2: DP > 80
bedtools annotate -counts -files filter1.selectGenoINDEL.Calmmq60.diffs.vcf -i batrachochytrium_dendrobatidis_1_transcripts_CDS2.gtf > indel.Calmmq60.CDS.countsGTF.txt
bedtools annotate -counts -files filter0.selectGenoINDEL.Calmmq60.diffs.vcf -i batrachochytrium_dendrobatidis_1_transcripts_CDS2.gtf > filter0.indel.Calmmq60.CDS.countsGTF.txt
bedtools annotate -counts -files filter2.selectGenoINDEL.Calmmq60.diffs.vcf -i batrachochytrium_dendrobatidis_1_transcripts_CDS2.gtf > filter2.indel.Calmmq60.CDS.countsGTF.txt
###############################################################################
###############################################################################
setwd("~/Google Drive/Experiments/bd_penny/reanalyzeHapCall_Aug15/")
#library("VariantAnnotation")
rm(list=ls())
#
gtf = read.table("batrachochytrium_dendrobatidis_1_transcripts.CDS.gtf", sep="\t")
head(gtf)
gtf$V1 = gsub("Supercontig_1","bden_JEL423_supercont1",gtf$V1)
write.table(gtf,"batrachochytrium_dendrobatidis_1_transcripts.CDS.gtf", sep="\t",quote=F,col.names=F,row.names=F)
#
gtf = gtf[,c(1,4,5)]
########################
#
setwd("~/Google Drive/Experiments/bd_penny/reCall0815_newRef/")
rm(list=ls())
###############
# g=read.table(file="~/Experiments/bd_reseq/batrachochytrium_dendrobatidis_1_genome_summary_per_gene.txt",sep="\t",header=T,quote="",row.names=NULL, stringsAsFactors = F)
# cnames = colnames(g); cnames = cnames[-1] ; cnames = c(cnames, "null"); colnames(g) = cnames;
# # head(g); 
# dim(g)
# g.names = g[,1]
# g=g[,c(5,6,9)]
# g = as.matrix(g)
anno = read.table("~/Google Drive/Experiments/bd_carter/Bd_annotation_files/batrachochytrium_dendrobatidis_1_genome_summary_per_gene_annotated.txt", sep="\t", header=T, stringsAsFactors = F,quote="")
anno2 = read.table("~/Google Drive/Experiments/bd_carter/Bd_annotation_files/Bd_GeneNames_annotation.txt", sep="\t", header=T, stringsAsFactors = F,quote="")
anno3 = read.table("~/Google Drive/Experiments/bd_penny/reanalyzeHapCall_Aug15/JEL423.interpro.tab", sep="\t", header=F, stringsAsFactors = F,fill=T,quote="")
transcriptGene = read.csv("~/Google Drive/Experiments/bd_penny/reanalyzeHapCall_Aug15/BD_trans_to_gene.tab",
                          header=F,stringsAsFactors=F,sep=" ",quote="")
anno3a = cbind(GeneID =transcriptGene$V2[match(anno3$V1,transcriptGene$V1)],anno3)
#  

###############
res = read.table("indel.Calmmq60.CDS.countsGFF3.txt",sep="\t", stringsAsFactors = F)
res = read.table("indel.Calmmq60.CDS.countsGTF.txt",sep="\t", stringsAsFactors = F)
res = read.table("filter2.indel.Calmmq60.CDS.countsGTF.txt",sep="\t", stringsAsFactors = F)
head(res)
# order is old, young
table(res$V10)
#
res2 = res[which(res$V10 > 0),]
res.genes = res$V9[which(res$V10 > 0)]
res.genes = sapply(res.genes, function(x) unlist(strsplit(x, split=";"))[1])
res.genes = sapply(res.genes, function(x) unlist(strsplit(x, split=" "))[2])
length(res.genes)
head(res.genes)
#
anno.sub = anno[which(anno$name %in% res.genes),]
anno2.sub = anno2[which(anno2$Broad_Gene_ID %in% res.genes),] 
#
anno.sub$name.1 = gsub("predicted protein","",anno.sub$name.1)
anno.sub$name.1 = gsub("conserved hypothetical protein","",anno.sub$name.1)
anno.sub$name.1 = gsub("hypothetical protein","",anno.sub$name.1)
#
write.table(anno.sub,"indel.Calmmq60.CDS.anno.txt", sep="\t",row.names=F)
write.table(anno2.sub,"indel.Calmmq60.CDS.anno2.txt", sep="\t",row.names=F)
#
############################
############################
############################
# read in the table
# nonSynSNPtable <-read.table("49strains_Nonsynonymous_UM142_anno_top10quan.txt",header=T,sep="\t",stringsAsFactors=F, quote="")
# # get ride of genes which empy lines
# genes <- subset(nonSynSNPtable$Broad_Gene_ID,nonSynSNPtable$Broad_Gene_ID != 'NA')

#get the list of all genes
# allgenes <- read.csv("BD_trans_to_gene.tab",
#                      header=F,stringsAsFactors=F,sep=" ",quote="")
# # just want the names of the universe of possible genes
# universe <- allgenes$V2


#universe
#genes
genes = unique(res.genes)
length(unique(genes))
length((genes))
# problem matching mode of this before
mode(universe)
mode(genes)

library("AnnotationDbi")
# define the table of GO terms for genes from this tab delimited file
godat <- read.table("~/Google Drive/Experiments/bd_penny/reanalyzeHapCall_Aug15/JEL423.IPR.GO",header=F, stringsAsFactors = F);
goframeData <- data.frame(godat$V1, godat$V2, godat$V3)
goFrame <- GOFrame(goframeData,organism="Batrachochytrium dendrobatidis")
goAllFrame=GOAllFrame(goFrame)

universe = unique(unique(godat$V3))

library("GSEABase")
library("GOstats")

gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())

params <- GSEAGOHyperGParams(name="My Custom GSEA based annot Params",
                             geneSetCollection=gsc,
                             geneIds = genes,
                             universeGeneIds = universe,
                             ontology = "MF",
                             pvalueCutoff = 0.05,
                             conditional = T,
                             testDirection = "over")


Over <- hyperGTest(params)
summary(Over)
Over.summ = summary(Over)
byGO = unlist(lapply(geneIdsByCategory(Over), function(x) paste(x, collapse = ";")) )
OverMF.summ = cbind(Over.summ, GeneIds = byGO[match(Over.summ$GOMFID,names(byGO))])

paramsCC <- GSEAGOHyperGParams(name="My Custom GSEA based annot Params",
                               geneSetCollection=gsc,
                               geneIds = genes,
                               universeGeneIds = universe,
                               ontology = "CC",
                               pvalueCutoff = 0.05,
                               conditional = T,
                               testDirection = "over")

OverCC <- hyperGTest(paramsCC)
summary(OverCC)
OverCC.summ = summary(OverCC)

paramsBP <- GSEAGOHyperGParams(name="My Custom GSEA based annot Params",
                               geneSetCollection=gsc,
                               geneIds = genes,
                               universeGeneIds = universe,
                               ontology = "BP",
                               pvalueCutoff = 0.05,
                               conditional = T,
                               testDirection = "over")

OverBP <- hyperGTest(paramsBP)
summary(OverBP)
OverBP
OverBP.summ = summary(OverBP)
byGO = unlist(lapply(geneIdsByCategory(OverBP), function(x) paste(x, collapse = ";")) )
OverBP.summ = cbind(OverBP.summ, GeneIds = byGO[match(OverBP.summ$GOBPID,names(byGO))])
#
write.table(OverBP.summ,"indel.Calmmq60.CDS.countsGTF_GOenrich.txt",sep="\t", row.names=F)
write.table(OverMF.summ,"indel.Calmmq60.CDS.countsGTF_GOMFenrich.txt",sep="\t", row.names=F)
#
write.table(OverBP.summ,"filter0.indel.Calmmq60.CDS.countsGTF_GOenrich.txt",sep="\t", row.names=F)
write.table(OverMF.summ,"filter0.indel.Calmmq60.CDS.countsGTF_GOMFenrich.txt",sep="\t", row.names=F)
####
# check out genes in interesting GO
# row 12, 8, 13
i=13
OverBP.summ[i,]
goi = OverBP.summ$GeneIds[i]
goi = unlist(strsplit(as.character(goi),split=";"))
#
annoSub = anno[match(goi,anno$name),]
anno2Sub = anno2[match(goi,anno2$Broad_Gene_ID),]
anno3aSub = anno3a[match(goi,anno3a$GeneID),]
# you can look for underrepresented too with the 'testDirection="under"'

##########


