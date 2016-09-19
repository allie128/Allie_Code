################################################
########################################################################
## MAKE NJ AND MP TREES FOR BD FLUIDIGM PAPER
rm(list=ls()) 
setwd("C:/Users/hji3/Desktop/Thesis Work/Fluidigm")


library(phangorn)
library(seqinr)
library(VariantAnnotation)

################################
###FOR WHOLE GENOME SNP DATA####
################################

#reads in the table with all of the whole genome isolate genotypes by snp
taball=read.table(file="C:/Users/hji3/Desktop/Thesis Work/Panama Isolate Genome Project/snps_oldandnewdata.tab",sep="\t",header=T, stringsAsFactors=F)
#check the data
apply(taball[,3:7], 2, function(x) table(x))

save.taball = taball
taball = na.omit(taball)


chromall = sapply(taball$CHROM, function(x) unlist(strsplit(x, split="\\."))[2])
chromall = as.numeric(chromall)


dir.create("Figs_forBDFl")


#subset to only the 25 shared strains that were included in the Bd fluidigm run and the whole genome sequencing 

taball_25 <- taball[,-c(5,6,8,12,13,14,15,16,22,25)]


# GET Neighbor Joining TREE, ROOTED WITH UM142
colnames(taball_25)
WG_isolate_list <- colnames(taball_25)[3:27]
snp <- as.matrix(taball_25[,-c(1:2)])
colnames(snp)
snp.dist = dist(t(snp))
nuc.nj = nj(snp.dist)
nuc.nj.root = root(phy=nuc.nj,outgroup="UM142", resolve.root=F)
#add bootstrap
snp.pD = as.phyDat(t(snp), type="USER", levels=c(0,1,2))
NJtrees.boot_WG <- bootstrap.phyDat(snp.pD, FUN=function(x)NJ(snp.dist), bs=100)
NJtreeBS_WG <- plotBS(nuc.nj.root,NJtrees.boot_WG, type="phylogram")

write.tree(NJtreeBS_WG, file="Figs_forBDFl/nj_1_boot_WG.tree")
pdf(paste("Figs_forBDFl/njContigs/nj_1_WG","-",nrow(snp),".pdf",sep=""), width=10)



plot.phylo(nuc.nj.root, cex=.8)
dev.off()


#ML for WG SNP data
fit_WG <- pml(nuc.nj.root, snp.pD)
fit_WG <- optim.pml(fit_WG, model="GTR", optInv=TRUE, rearrangements="NNI",control = pml.control(trace = 0))
fit_WG <- optim.pml(fit_WG, model="GTR", optInv=TRUE, optGamma=TRUE,rearrangement = "stochastic", control = pml.control(trace = 0))




bs <- bootstrap.pml(fit_WG, bs=100, optNni=TRUE)

treeBS <- plotBS(fit_WG$tree,bs)
write.tree(treeBS, file="Figs_forBDFl/ML_1_boot_WG.tree")

#MP for WG SNP data
cost.m = matrix(data=c(0,1,2,1,0,1,2,1,0),nrow=3)
treePR = pratchet(data=snp.pD, start=nuc.nj.root, k=10, cost=cost.m, method="sankoff", trace=F)
treeAC = acctran(tree=treePR, data=snp.pD)
pdf(paste("Figs/mp_1","-",nrow(snp),".pdf",sep=""), width=10)
plot.phylo(treePR, cex=.8)
dev.off()
write.tree(treeAC, file="Figs/mp_1_WG.tree")
### ADD BOOTSTRAP 
boot.trees_WG = bootstrap.phyDat(snp.pD, FUN=function(x)pratchet(data=x, start=nuc.nj.root, k=10, cost=cost.m, method="sankoff", trace=F), bs=100)

pdf(paste("Figs/mp_boot_1_WG","-",nrow(snp),".pdf",sep=""), width=10)

MP_WG <- plotBS(treeAC, boot.trees_WG, type="ph", cex=.8, bs.adj=1)

dev.off()

write.tree(MP_WG, file="Figs/mp_boot_1_WG.tree")


#####################
####FLUIDIGM DATA####
#####################

#read in concatenated consensus sequences and convert to phyDat
myseqs <- read.alignment(file.choose(), format = "fasta")
seq.pD <- as.phyDat(myseqs)
#calculate distance matrix for phyDat
mat <- dist.alignment(myseqs, matrix = "identity")

#make NJ tree for Fluidigm data
nuc.nj.mat = nj(mat)
nuc.nj.root_FL = root(phy=nuc.nj.mat,outgroup="UM142", resolve.root=F)
plot.phylo(nuc.nj.root_FL, cex=.8)
#add bootstraps
NJtrees.boot <- bootstrap.phyDat(seq.pD, FUN=function(x)NJ(mat), bs=100)
treeNJ <- plotBS(nuc.nj.root, NJtrees.boot, "phylogram")


#ML for Fluidigm data
fit_FL <- pml(nuc.nj.root_FL, seq.pD)
fit_FL <- optim.pml(fit_FL, rearrangements="NNI")
bs_FL <- bootstrap.pml(fit_FL, bs=100, optNni=TRUE)
treeBS_FL <- plotBS(fit_FL$tree,bs_FL)

write.tree(treeBS, file="Figs_forBDFl/ML_1_boot_FL.tree")

#MP for Fluidigm data
cost.m = matrix(data=c(0,1,2,1,0,1,2,1,0),nrow=3)
treePR = pratchet(data=seq.pD, start=nuc.nj.root_FL, k=100, cost=cost.m, method="sankoff", trace=F)
plot.phylo(treePR, cex=.8)
#add bootstrap
MP.boot.trees = bootstrap.phyDat(seq.pD, FUN=function(x)pratchet(data=x, start=nuc.nj.root_FL, k=10, cost=cost.m, method="sankoff", trace=F), bs=100)

treeMP_FL <- plotBS(treePR, MP.boot.trees, type="ph", cex=.8, bs.adj=1)
write.tree(treeMP_FL, file="Figs_forBDFl/mp_1_boot_FL.tree")

#MP for Fluidigm data


treePR_FL = pratchet(data=seq.pD, start=nuc.nj.root_FL, k=10, cost=cost.m, method="sankoff", trace=F)
treeAC_FL = acctran(tree=treePR_FL, data=seq.pD)

plot.phylo(treeAC_FL, cex=.8)

dev.off()
write.tree(treeAC, file="Figs_forBDFl/mp_1_FL.tree")
### ADD BOOTSTRAP 
boot.trees_MPFL = bootstrap.phyDat(seq.pD, FUN=function(x)pratchet(data=x, start=nuc.nj.root_FL, k=10, method="fitch", trace=F), bs=100)
pdf(paste("Figs_forBDFl/mp_boot_1_FL","-",nrow(snp),".pdf",sep=""), width=10)

FL_MP <- plotBS(treeAC_FL, boot.trees_MPFL, type="ph", cex=.8, bs.adj=1)


write.tree(FL_MP, file="Figs_forBDFl/mp_boot_1_FL.tree")










#########################################
###MAKE COPHYLOGENY USING PHYTOOLS#######
#########################################

install.packages('phytools')
library(phytools)



ML_WG <- as.phylo(ML_boot_WG)
ML_FL <- as.phylo(treeBS_FL)

FL_names <- ML_FL$tip.label[order(ML_FL$tip.label)]
WG_names <- ML_WG$tip.label[order(ML_WG$tip.label)]

assoc <- as.matrix(cbind(WG_names, FL_names))


#ML trees

cophy_ML<-cophylo(ML_WG,ML_FL,assoc=assoc)

plot(cophy_ML, fsize=0.4) 


#make node labels less than 70 blank
ML_WG$node.label[which(ML_WG$node.label < 70)] <- " "
ML_FL$node.label[which(ML_FL$node.label < 70)] <- " "



nodelabels.cophylo(ML_WG$node.label, which="left",adj=c(1,-0.2),frame="none", cex=0.4)
nodelabels.cophylo(ML_FL$node.label, which="right",adj=c(-0.2,-0.7),frame="none", cex=0.4)

#MP trees
plotBS(MP_WG, type="ph")
plotBS(MP_FL, type="ph")

boot.trees_WG


MP_WG <- as.phylo(MP_WG)

MP_FL <- as.phylo(FL_MP)

cophy_MP<-cophylo(MP_WG,MP_FL,assoc=assoc)

MP_WG$node.label[which(MP_WG$node.label < 70)] <- " "
MP_FL$node.label[which(MP_FL$node.label < 70)] <- " "
h1<-max(nodeHeights(MP_WG))
h2<-max(nodeHeights(MP_FL))

plot(cophy_MP, fsize=0.4,scale.bar=round(0.5*c(h1,h2),2)) 

testnex <- as.phylo(file.choose())

nodelabels.cophylo(MP_WG$node.label, which="left",adj=c(1,-0.2),frame="none", cex=0.4)
nodelabels.cophylo(MP_FL$node.label, which="right",adj=c(-0.2,-0.7),frame="none", cex=0.4)
