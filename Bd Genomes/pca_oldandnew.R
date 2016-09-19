## pca for old and new snps together


setwd("C:/Users/hji3/Desktop/Thesis Work/Panama Isolate Genome Project")
library(VariantAnnotation)

#import the raw vcfs from this experiment

vcfsnp <- readVcf(file.choose(), "BdSNP")
#parse the data to match that of the old dataset
vcfsnp.geno <- geno(vcfsnp)$GT


##reformatting my data to match the old data
test2 <- sapply(1:nrow(vcfsnp.geno), function(x) strsplit(rownames(vcfsnp.geno)[x], split=":"))
test3 <- sapply(1:length(test2), function(x) strsplit(test2[[x]][2], split="_"))
test4 <- sapply(1:length(test2), function(x) sub('Supercontig_', 'bden_JEL423_supercont', test2[[x]][1]))
output <- matrix(unlist(test3), ncol = 2, byrow = TRUE)
newdata <- cbind(test4,output[,1],vcfsnp.geno[,1],vcfsnp.geno[,2],vcfsnp.geno[,3],vcfsnp.geno[,4],vcfsnp.geno[,5],vcfsnp.geno[,6])
#add column names
colnames(newdata) <- c("CHROM","POS","Campana","JEL410","JEL412","JEL413","Rio_Maria","Sora")
#rename genotypes to integers to match old data

newdata1 <- t(sapply(1:nrow(newdata), function(x) sub('0/0', '0', newdata[x,])))
newdata1 <- t(sapply(1:nrow(newdata1), function(x) sub('0/1', '1', newdata1[x,])))
newdata1 <- t(sapply(1:nrow(newdata1), function(x) sub('1/1', '2', newdata1[x,])))
#gets rid of the other alt genotypes
newdata1 <- t(sapply(1:nrow(newdata1), function(x) sub('0/2', '.', newdata1[x,])))
newdata1 <- t(sapply(1:nrow(newdata1), function(x) sub('1/2', '.', newdata1[x,])))
newdata1 <- t(sapply(1:nrow(newdata1), function(x) sub('2/2', '.', newdata1[x,])))
#turns all missing data into NA
newdata1[newdata1=="."] <- NA

newdatadf <- data.frame(newdata1)

##write out formatted new snps file
write.table(newdata1,file="snps_6strains_oldformat_altrem.tab",sep="\t",quote=F,row.names=F)

newdata1 <- read.table(file.choose())

#check the genotypes
apply(newdata1[,3:7], 2, function(x) table(x))
apply(old[,3:31], 2, function(x) table(x))


write.table(tab2,file="testjoin1.tab",sep="\t",quote=F,row.names=F)

# add the old snp data set

old=read.table(file="Bd_29.selectedSNPs.5.GT.4X.tab",sep="\t",header=T)
head(old)
dim(old)
#old.sites = t(sapply(1:nrow(old), function(x) paste(c(old$CHROM[x],old$POS[x]) )))
oldsites <- paste(old$CHROM, old$POS)
newsites <- paste(newdatadf$CHROM, newdatadf$POS)

newdatadf <- data.frame(newdata1)

###
# find overlapping SNP sites

matches1 <- which(newsites %in% oldsites)
matches2 <- which(oldsites %in% newsites)
##join the overlapping sites
tab4 = cbind(newdata1[matches1,], old[matches2,])
#gets rid of the extra CHROM and POS columns
tab4 = tab4[,-c(9,10)]

#write it out
write.table(tab4,file="snps_oldandnewdata.tab",sep="\t",quote=F,row.names=F)



###############################################
###############################################
# PCA
rm(list=ls()) 
install.packages('gplots')
library(gplots)
setwd("C:/Users/hji3/Desktop/Thesis Work/Panama Isolate Genome Project")
snp=read.table(file="snps_oldandnewdata.tab",sep="\t",header=T, stringsAsFactors=F)
tail(snp)

## take out divergent strains 
colnames(snp)
#take out LFT001_10, UM142, CLFT024.02
snp = snp[,-c(c(31,37,15))]

snp2 = newdatadf
snp.omit = na.omit(snp2)
head(snp.omit); dim(snp.omit)
as.data.frame(apply(snp.omit[,-c(1,2)], 2, function(x) table((x))) )

y.snp=t(snp.omit[,3:ncol(snp.omit)])
colnames(y.snp)=apply(snp.omit,1,function(x) paste(x[1],x[2]))


y.snp.scale=scale(y.snp)

sc <- attr(y.snp.scale, "scaled:scale")
#remove the invariant loci to each other - these may be present if some strains are removed (eg JEL423, JAM81)
length(which(sc==0)) # num of invariant loci
#195
y.snp.scale=y.snp.scale[,which(sc!=0)]

y.snp.pc=prcomp(y.snp.scale)





# pdf(file="pca_snps_6strains.select.2.GT.oldData.pdf",height=7,width=7,onefile=T)

par(mfrow=c(1,1))

plot(predict(y.snp.pc)[,1],predict(y.snp.pc)[,2],pch=" ",xlab="pc1",ylab="pc2",main="PCA for Historic and Contemporary Isolates",sub=paste(ncol(y.snp.scale),"loci"))
#plot(predict(y.snp.pc)[,1],predict(y.snp.pc)[,2],pch=" ",xlab="pc1",ylab="pc2",xlim=c(15,30), ylim=c(-20,20), main="PCA",sub=paste(ncol(y.snp.scale),"loci"))
text(predict(y.snp.pc)[,1],predict(y.snp.pc)[,2],labels=rownames(predict(y.snp.pc)),cex=.7, font=2)
#text(predict(y.snp.pc)[,1],predict(y.snp.pc)[,2],labels=rownames(predict(y.snp.pc)),cex=.3)
#dev.off()

# check the PC scores
y.snp.pc$x
#

plot(y.snp.pc)
textplot(round(predict(y.snp.pc)[sort.list(predict(y.snp.pc)[,1]),1:2],digits=4))



