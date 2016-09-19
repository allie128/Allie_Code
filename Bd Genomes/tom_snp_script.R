############################################
# add the old snp data set
#      MAYBE NOT GOOD WAY
old=read.table(file="/mnt/home/tpoorten/bd_reseq/round49/Bd_29.selectedSNPs.5.GT.4X.tab",sep="\t",header=T)
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

###############################################
###############################################
# PCA
rm(list=ls()) 
library(gplots)
setwd("~/Dropbox/bd_carter/")
snp=read.table(file="snps_6strains.select.2.GT.oldData.tab",sep="\t",header=T, stringsAsFactors=F)
tail(snp)
#remove mitochondria
snp = snp[which(snp$CHROM != "bden_JEL423_MT"),]
head(snp);dim(snp)
## take out divergent strains 
colnames(snp)
snp = snp[,-c(c(31,37,15))]

snp2 = snp
snp.omit = na.omit(snp2)
head(snp.omit); dim(snp.omit)
as.data.frame(apply(snp.omit[,-c(1,2)], 2, function(x) table((x))) )

y.snp=t(snp.omit[,3:ncol(snp.omit)])
colnames(y.snp)=apply(snp.omit,1,function(x) paste(x[1],x[2]))
y.snp.scale=scale(y.snp)
sc <- attr(y.snp.scale, "scaled:scale")
#remove the invariant loci to each other - these may be present if some strains are removed (eg JEL423, JAM81)
length(which(sc==0)) # num of invariant loci
y.snp.scale=y.snp.scale[,which(sc!=0)]

y.snp.pc=prcomp(y.snp.scale)

# pdf(file="pca_snps_6strains.select.2.GT.oldData.pdf",height=7,width=7,onefile=T)

par(mfrow=c(1,1))

plot(predict(y.snp.pc)[,1],predict(y.snp.pc)[,2],pch=" ",xlab="pc1",ylab="pc2",main="PCA",sub=paste(ncol(y.snp.scale),"loci"))
text(predict(y.snp.pc)[,1],predict(y.snp.pc)[,2],labels=rownames(predict(y.snp.pc)),cex=.7)
#dev.off()

# check the loadings
y.snp.pc$x
#
plot(y.snp.pc)
textplot(round(predict(y.snp.pc)[sort.list(predict(y.snp.pc)[,1]),1:2],digits=4))

################################################
########################################################################
## MAKE NJ AND MP TREES
rm(list=ls()) 
setwd("~/Google Drive/Experiments/bd_carter/")
tab=read.table(file="snps_6strains.select.2.GT.oldData.tab",sep="\t",header=T, stringsAsFactors=F)
apply(tab[,3:7], 2, function(x) table(x))

library(phangorn)
save.tab = tab
tab = na.omit(tab)

chrom = sapply(tab$CHROM, function(x) unlist(strsplit(x, split="\\."))[2])
chrom = as.numeric(chrom)

dir.create("Figs")
dir.create("Figs/njContigs")
dir.create("Figs/mpContigs")

# GET NJ TREE, ROOTED WITH UM142
colnames(tab)
snp <- as.matrix(tab[,-c(1:2,3:4,8)])
colnames(snp)
snp.dist = dist(t(snp))
nuc.nj = nj(snp.dist)
nuc.nj.root = root(phy=nuc.nj,outgroup="UM142", resolve.root=F)

write.tree(nuc.nj.root, file="Figs/nj_penny.tree")
pdf(paste("Figs/njContigs/nj_penny","-",nrow(snp),".pdf",sep=""), width=10)
plot.phylo(nuc.nj.root, cex=.8)
dev.off()


#MP - nuclear
snp.pD = as.phyDat(t(snp), type="USER", levels=c(0,1,2))
cost.m = matrix(data=c(0,1,2,1,0,1,2,1,0),nrow=3)
treePR = pratchet(data=snp.pD, start=nuc.nj.root, k=100, cost=cost.m, method="sankoff", trace=F)
treeAC = acctran(tree=treePR, data=snp.pD)
pdf(paste("Figs/mp_penny","-",nrow(snp),".pdf",sep=""), width=10)
plot.phylo(treeAC, cex=.8)
dev.off()
write.tree(treeAC, file="Figs/mp_penny.tree")
### ADD BOOTSTRAP 
boot.trees = bootstrap.phyDat(snp.pD, FUN=function(x)pratchet(data=x, start=nuc.nj.root, k=100, cost=cost.m, method="sankoff", trace=F), bs=100, multicore=TRUE)
pdf(paste("Figs/mp_boot_penny","-",nrow(snp),".pdf",sep=""), width=10)
plotBS(treeAC, boot.trees, type="ph", cex=.8, bs.adj=1)
dev.off()
boot = plotBS(treeAC, boot.trees, type="ph", cex=.8, bs.adj=1)
write.tree(boot, file="Figs/mp_boot_penny.tree")


#

i = 1
for(i in 1:15){
  #NJ
  snp <- as.matrix(tab[which(i == chrom),-c(1:2)])
  nuc.nj = nj(dist(t(snp)))
  nuc.nj.root = root(phy=nuc.nj,outgroup="UM142", resolve.root=F)
  pdf(paste("Figs/njContigs/nj_penny_chrom",i,"-",nrow(snp),".pdf",sep=""), width=10)
  plot.phylo(nuc.nj.root, cex=.8)
  dev.off()
  
  # MP
  snp.pD = as.phyDat(t(snp), type="USER", levels=c(0,1,2))
  cost.m = matrix(data=c(0,1,2,1,0,1,2,1,0),nrow=3)
  treePR = pratchet(data=snp.pD, start=nuc.nj.root, k=100, cost=cost.m, method="sankoff", trace=F)
  treeAC = acctran(tree=treePR, data=snp.pD)
  pdf(paste("Figs/mpContigs/mp_penny_chrom",i,"-",nrow(snp),".pdf",sep=""), width=10)
  plot.phylo(treeAC, cex=.8)
  dev.off()
}

########################
# ADD BOOTSTRAP VALUES

#MP - nuclear
snp.pD = as.phyDat(t(snp), type="USER", levels=c(0,1,2))
cost.m = matrix(data=c(0,1,2,1,0,1,2,1,0),nrow=3)
treePR = pratchet(data=snp.pD, start=nuc.nj.root, k=100, cost=cost.m, method="sankoff", trace=F)
treeAC = acctran(tree=treePR, data=snp.pD)
### ADD BOOTSTRAP 
boot.trees = bootstrap.phyDat(snp.pD, FUN=function(x)pratchet(data=x, start=nuc.nj.root, k=100, cost=cost.m, method="sankoff", trace=F), bs=2)
plotBS(treeAC, boot.trees, type="ph", cex=.6)




########################################################################
########################################################################
## Drop superyoung 4/15/15
rm(list=ls()) 
setwd("~/Google Drive/Experiments/bd_penny//")
library(phangorn)
tree = read.tree("Figs/mp_boot_penny.tree")
plot.phylo(tree)
tree = drop.tip(phy = tree, tip = "JEL427Superyoung")
plot.phylo(tree)
write.tree(tree, file="Figs/mp_boot_penny_2strains.tree")
#
#
#####################################################
# SNPEFF

####
# PREP snpeff database
grep CDS batrachochytrium_dendrobatidis_1_transcripts.gtf > batrachochytrium_dendrobatidis_1_transcripts_CDS.gtf 
grep -v exon batrachochytrium_dendrobatidis_1_transcripts_CDS.gtf > batrachochytrium_dendrobatidis_1_transcripts_CDS2.gtf 
#
###
# downloaded new version of snpEff (8/7/15)
# put gff3 file and genome seq file in the dir: snpeff/data/Bd1.1/ OR USE GTF FILE
# gff was used by Jason for the snpeff reference build
java -jar snpEff.jar build -gff3 -v Bd1.1  
#java -jar snpEff.jar build -gtf22 -v Bd1.1

java -jar snpEff.jar -v -d -c snpEff.config -s snps.youngold.html Bd1.1 ../filter.selectGenoSNP.Cal.diffs.vcf > snpeff.youngold.out.xt
#####################