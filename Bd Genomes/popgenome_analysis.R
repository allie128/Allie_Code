##Pipeline for checking and filtering raw SNPs and indels from GATK

source("http://bioconductor.org/biocLite.R")
biocLite()

setwd("C:/Users/hji3/Desktop/Thesis Work/Panama Isolate Genome Project")
library(VariantAnnotation)

#try out for analyzing filtered vcf
require(PopGenome)

#import the filtered vcf
vcfsnp <- readVcf(file.choose(), "BdSNP")

vcfsnp.class.slide <- sliding.window.transform(vcfsnp,1000,1000,type=2)






#pull out the genotype information
vcfsnp.geno <- geno(vcfsnp)$GT
#check the number of snps before any filtering
dim(vcfsnp.geno)
#63237
#check distribution of genotypes
apply(vcfsnp.geno, 2, table)

### turn low GQ into "." missing data <- THIS IS NEW 7/19/15
vcfsnp.gq = geno(vcfsnp)$GQ
hist(vcfsnp.gq, xlim=c(0,100), breaks=20)
apply(vcfsnp.gq,2,function(x) table(x<99))
apply(vcfsnp.gq,2,function(x) table(x<20))
apply(vcfsnp.gq,2,function(x) table(x<40))
table(vcfsnp.gq)
table(vcf.geno)
#filters out genotype quality less than 20
vcfsnp.geno[which(vcfsnp.gq < 20)] = "."
#changes the original vcf genotype
geno(vcfsnp)$GT = vcfsnp.geno

### filter by read depth
vcfsnp.dp = geno(vcfsnp)$DP
apply(vcfsnp.dp, 2, function(x) mean(x, na.rm = T))
apply(vcfsnp.dp, 2, function(x) median(x, na.rm = T))
apply(vcfsnp.dp, 2, function(x) sd(x, na.rm = T))*5

# cutoff of max DP
apply(vcfsnp.dp, 2, function(x) median(x, na.rm = T)) + apply(vcfsnp.dp, 2, function(x) sd(x, na.rm = T))*6
cutoffs <- apply(vcfsnp.dp, 2, function(x) median(x, na.rm = T)) + apply(vcfsnp.dp, 2, function(x) sd(x, na.rm = T))*6
apply(vcfsnp.dp, 2, table)
hist(vcfsnp.dp,xlim=c(0,30),breaks=200)
summary(vcfsnp.dp)

#filters out based on max DP calculated earlier as (median+ SD)*6
vcfsnp.geno[which(vcfsnp.dp[,1] > cutoffs[1])] = "."
vcfsnp.geno[which(vcfsnp.dp[,2] > cutoffs[2])] = "."
vcfsnp.geno[which(vcfsnp.dp[,3] > cutoffs[3])] = "."
vcfsnp.geno[which(vcfsnp.dp[,4] > cutoffs[4])] = "."
vcfsnp.geno[which(vcfsnp.dp[,5] > cutoffs[5])] = "."
vcfsnp.geno[which(vcfsnp.dp[,6] > cutoffs[6])] = "."


#filters out genotypes with a read depth less than 4
vcfsnp.geno[which(vcfsnp.dp < 4)] = "."

#changes original vcf genotype
table(geno(vcfsnp)$GT) #check that numbers change
geno(vcfsnp)$GT = vcfsnp.geno
table(geno(vcfsnp)$GT) #check that numbers change


## explore the distribution of variables for filter creation
vcf.info <- info(vcfsnp)
hist(vcf.info$InbreedingCoeff, breaks=20) # 
hist(vcf.info$MQRankSum, breaks=200) # < -3 or -2 (maxlenRef: < -5)
hist(vcf.info$MQRankSum, breaks=200,xlim=c(-5,5)) # < -3 or -2 (maxlenRef: < -5)
hist(vcf.info$QD, breaks=100) # MAYBE TRY A RANGE OF VALUES HERE - FROM 2 TO 10
summary(vcf.info$QD)
hist(vcf.info$DP, breaks=400) # 
hist(vcf.info$DP, breaks=400,xlim=c(0,100)) # 
hist(vcf.info$ReadPosRankSum, breaks=20) # < -8
hist(vcf.info$ReadPosRankSum, breaks=200,xlim=c(-10,10)) # < -8
hist(vcf.info$MQ, breaks =50) # < 100 or 200
hist(vcf.info$SOR, breaks =50) # 
#
hist(vcf.rowData$QUAL, breaks=1000, xlim=c(0,400))
summary(vcf.rowData$QUAL)

########
filter1 = which(vcf.info$FS > 60)
filter2 = which(vcf.info$QD < 2)
filter3 = which(vcf.info$ReadPosRankSum < -8)
filter4 = which(vcf.info$MQ < 50)  #already filtered
filter5 = which(vcf.info$SOR > 4)
filter6 = which(vcf.info$MQRankSum < -3)
filters = unique(c(filter1, filter2, filter3, filter3, filter4, filter5, filter6))#, filter7, filter8, filter9))
length(filters)
# new object of filtered sites
dim(vcfsnp.geno)
vcfsnp = vcfsnp[-filters,]
vcfsnp.geno = geno(vcfsnp)$GT
head(vcfsnp.geno); dim(vcfsnp.geno)

#write out filtered snps
writeVcf(vcfsnp, "snp_filtered.vcf")


###subetting snps for pairwise comparisons
test <- subset(vcfsnp.geno,  vcfsnp.geno[,1] != ".")
test2 <- subset(test, test[,2] != ".")
Campana_410 <- test2[,1:2]
###
test3 <- subset(test, test[,3] != ".")
Campana_412 <- test3[,c(1,3)]
###
test4 <- subset(test, test[,4] != ".")
Campana_413 <- test4[,c(1,4)]
###
test5 <- subset(test, test[,5] != ".")
Campana_RMaria <- test5[,c(1,5)]
###
test6 <- subset(test, test[,6] != ".")
Campana_Sora <- test6[,c(1,6)]
#############
test <- subset(vcfsnp.geno,  vcfsnp.geno[,2] != ".")
test2 <- subset(test, test[,3] != ".")
J410_J412 <- test2[,c(2,3)]
###
test3 <- subset(test, test[,4] != ".")
J410_J413 <- test3[,c(2,4)]
###
test4 <- subset(test, test[,5] != ".")
J410_RMaria <- test4[,c(2,5)]
###
test5 <- subset(test, test[,6] != ".")
J410_Sora <- test5[,c(2,6)]
###########
test <- subset(vcfsnp.geno,  vcfsnp.geno[,3] != ".")
test2 <- subset(test, test[,4] != ".")
J412_J413 <- test2[,c(3,4)]
###
test3 <- subset(test, test[,5] != ".")
J412_RMaria <- test3[,c(3,5)]
###
test4 <- subset(test, test[,6] != ".")
J412_Sora <- test4[,c(3,6)]
##########
test <- subset(vcfsnp.geno,  vcfsnp.geno[,4] != ".")
test2 <- subset(test, test[,5] != ".")
J413_RMaria <- test2[,c(4,5)]
###
test3 <- subset(test, test[,6] != ".")
J413_Sora <- test3[,c(4,6)]
##########
test <- subset(vcfsnp.geno,  vcfsnp.geno[,5] != ".")
test2 <- subset(test, test[,6] != ".")
RMaria_Sora <- test2[,c(5,6)]


##% segregating
comparison1 <- which(J412_J413[,1] != J412_J413[,2])
comparison2 <- which(J412_RMaria[,1] != J412_RMaria[,2])
comparison3 <- which(J412_Sora[,1] != J412_Sora[,2])
comparison4 <- which(J413_RMaria[,1] != J413_RMaria[,2])
comparison5 <- which(J413_Sora[,1] != J413_Sora[,2])
comparison6 <- which(RMaria_Sora[,1] != RMaria_Sora[,2])


#########
#cn.mops
biocLite("cn.mops")
library(cn.mops)
setwd("C:/Users/hji3/Desktop/Thesis Work/Panama Isolate Genome Project/mergedBAMS")
BAMFiles <- list.files(pattern=".bam$")
bamDataRanges <- getReadCountsFromBAM(BAMFiles, mode="paired")
res <- cn.mops(bamDataRanges)

plot(res,which=1)


############
# count polymorphic sites between samples (ignore ref)
head(apply(vcfsnp.geno,1,function(x) x[1] != x[2]))

#try to pull out indicies of SNPS with all of the same GT
uninterestingsnps <- which(apply(vcfsnp.geno, 1, function(x) x[1] == x[2] & x[2] == x[3] & x[3]== x[4] & x[4] == x[5] & x[5] == x[6] ))


table(apply(vcfsnp.geno,1,function(x) x[1] != x[2]))

table(apply(vcfsnp.geno,1,function(x) x[1] != x[2] & x[1] != "." & x[2] != "."))
############

############
#
vcf.rowData = rowData(vcfsnp)
vcf.info = info(vcfsnp)
vcf.infoHeader = info(header(vcfsnp))
vcf.genoHeader = geno(header(vcfsnp))
#
as.matrix(vcf.infoHeader)
head(vcf.info)
contigs = rownames(vcf.info)
head(contigs)
contigs1 = sapply(contigs, function(x) unlist(strsplit(x, split=":"))[1])
#how many snps are on each contig?
table(contigs1)

# get median coverage per chromosome
head(vcfsnp.dp)
#
chrom.cov = apply(vcfsnp.dp, 2, function(x) tapply(x, contigs1, function(x) median(x, na.rm = T)))
chrom.covMean = apply(chrom.cov,2,function(x) median(x, na.rm=T))
chrom.covNorm = t(apply(chrom.cov, 1, function(x) x / chrom.covMean))
#
head(contigs1)
length(contigs1)
dim(chrom.cov)
#chrom.covMeanChrom1 = apply(chrom.cov[which(contigs1 == "bden_JEL423_supercont1.1"),],2,function(x) mean(x, na.rm=T))
chrom.covNormByChrom1 = t(apply(chrom.cov, 1, function(x) x / chrom.cov[2,]))




#####################
#SNPEFF

# add the old snp data set
#      MAYBE NOT GOOD WAY
old=read.table(file.choose(),header=T)
head(old)
dim(old)
old.sites = sapply(1:nrow(old), function(x) paste(c(old$CHROM[x],old$POS[x]), collapse="-") )
head(old.sites)
tab.sites = sapply(1:nrow(tab), function(x) paste(c(tab$CHROM[x],tab$POS[x]), collapse="-") )
head(tab.sites)
###
# find overlapping SNP sites
tab2 = cbind(tab,old[match(tab.sites, old.sites),])
head(tab2)
tail(tab2)
tab2 = tab2[,-c(9,10)]
table(is.na(tab2$CJB4))
write.table(tab2,file="snps_6strains.select.2.GT.oldData.tab",sep="\t",quote=F,row.names=F)




