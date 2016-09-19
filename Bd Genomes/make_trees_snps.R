################################################
########################################################################
## MAKE NJ AND MP TREES
rm(list=ls()) 
setwd("C:/Users/hji3/Desktop/Thesis Work/Panama Isolate Genome Project")

tab=read.table(file="snps_6strains_oldformat_altrem.tab",sep="\t",header=T, stringsAsFactors=F)
taball=read.table(file="snps_oldandnewdata.tab",sep="\t",header=T, stringsAsFactors=F)

apply(tab[,3:7], 2, function(x) table(x))

library(phangorn)
save.tab = tab
save.taball = taball
tab = na.omit(tab)
taball = na.omit(taball)

chrom = sapply(tab$CHROM, function(x) unlist(strsplit(x, split="\\."))[2])
chromall = sapply(taball$CHROM, function(x) unlist(strsplit(x, split="\\."))[2])
chrom = as.numeric(chrom)
chromall = as.numeric(chromall)

dir.create("Figs")
dir.create("Figs/njContigs")
dir.create("Figs/mpContigs")

# GET Neighbor Joining TREE, ROOTED WITH UM142
colnames(tab)
snp <- as.matrix(taball[,-c(1:2)])
colnames(snp)
snp.dist = dist(t(snp))
nuc.nj = nj(snp.dist)
nuc.nj.root = root(phy=nuc.nj,outgroup="UM142", resolve.root=F)

write.tree(nuc.nj.root, file="Figs/nj_1.tree")
pdf(paste("Figs/njContigs/nj_1","-",nrow(snp),".pdf",sep=""), width=10)
plot.phylo(nuc.nj.root, cex=.8)
dev.off()


#MP - nuclear
snp.pD = as.phyDat(t(snp), type="USER", levels=c(0,1,2))
cost.m = matrix(data=c(0,1,2,1,0,1,2,1,0),nrow=3)
treePR = pratchet(data=snp.pD, start=nuc.nj.root, k=100, cost=cost.m, method="sankoff", trace=F)
treeAC = acctran(tree=treePR, data=snp.pD)
pdf(paste("Figs/mp_1","-",nrow(snp),".pdf",sep=""), width=10)
plot.phylo(treeAC, cex=.8)
dev.off()
write.tree(treeAC, file="Figs/mp_1.tree")
### ADD BOOTSTRAP 
boot.trees = bootstrap.phyDat(snp.pD, FUN=function(x)pratchet(data=x, start=nuc.nj.root, k=100, cost=cost.m, method="sankoff", trace=F), bs=100)
pdf(paste("Figs/mp_boot_1","-",nrow(snp),".pdf",sep=""), width=10)
plotBS(treeAC, boot.trees, type="ph", cex=.8, bs.adj=1)
dev.off()
boot = plotBS(treeAC, boot.trees, type="ph", cex=.8, bs.adj=1)
write.tree(boot, file="Figs/mp_boot_1.tree")


#

i = 1
for(i in 1:15){
  #NJ
  snp <- as.matrix(tab[which(i == chrom),-c(1:2)])
  nuc.nj = nj(dist(t(snp)))
  nuc.nj.root = root(phy=nuc.nj,outgroup="UM142", resolve.root=F)
  pdf(paste("Figs/njContigs/nj_chrom",i,"-",nrow(snp),".pdf",sep=""), width=10)
  plot.phylo(nuc.nj.root, cex=.8)
  dev.off()
  
  # MP
  snp.pD = as.phyDat(t(snp), type="USER", levels=c(0,1,2))
  cost.m = matrix(data=c(0,1,2,1,0,1,2,1,0),nrow=3)
  treePR = pratchet(data=snp.pD, start=nuc.nj.root, k=100, cost=cost.m, method="sankoff", trace=F)
  treeAC = acctran(tree=treePR, data=snp.pD)
  pdf(paste("Figs/mpContigs/mp_chrom",i,"-",nrow(snp),".pdf",sep=""), width=10)
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