################################################
########################################################################
## MAKE NJ AND MP TREES
rm(list=ls()) 
setwd("C:/Users/hji3/Desktop/Thesis Work/Fluidigm/Results from run Aug 2015/Bdfluidigm.reduced.tar/Bdfluidigm.reduced/Bdfluidigm.reduced/ambiguities.split_amplicon")

#tab=read.table(file="snps_6strains_oldformat_altrem.tab",sep="\t",header=T, stringsAsFactors=F)
#taball=read.table(file="snps_oldandnewdata.tab",sep="\t",header=T, stringsAsFactors=F)


source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")

library(Biostrings)
library(GenomicAlignments)
library(ShortRead)
library(phangorn)


#read in data for primer 1
#p1data <- read.dna(file.choose(), format = "fasta")
#p4data <- read.dna(file.choose(), format = "fasta")

#labelsp1 <- sapply(names(p1data), function(x) unlist(strsplit(x, split=":"))[1])
#labelsp4 <- sapply(labels(p4data), function(x) unlist(strsplit(x, split=":"))[1])

#matches <- match(labelsp4, labelsp1)
#intersect(labelsp4, labelsp1)
#union(labelsp4, labelsp1)

##write for loop to find max number of samples
setwd("C:/Users/hji3/Desktop/Thesis Work/Fluidigm/Results from run Aug 2015/Bdfluidigm.reduced.tar/Bdfluidigm.reduced/Bdfluidigm.reduced/ambiguities.split_amplicon _edited")
files <- list.files(".")
files <- grep("merged",files,value=TRUE)
##checks that they are all the same length
i=1
for(i in 1:length(files)){
  pdata <- read.dna(file=files[i], format = "fasta")
  #pdata1 <- readDNAStringSet(file=files[i])
  print(files[i])
  print(i)
  print(nrow(pdata))
}

##the longest set(most number of samples)
savepdata1 <- readDNAStringSet(file=files[35])
savelabels <- sapply(names(savepdata1), function(x) unlist(strsplit(x, split=":"))[1])

##
##reads in all the files and concatenates them together
for(i in 1:length(files)){
  pdata <- read.dna(file=files[i], format = "fasta")
  #pdata1 <- readDNAStringSet(file=files[i])
  if(!is.null(nrow(pdata))){
    pdata1 <- readDNAStringSet(file=files[i]) 
    labelsp1 <- sapply(names(pdata1), function(x) unlist(strsplit(x, split=":"))[1])
    if(i!=35){
      newset <- DNAStringSet(rep(paste(rep("-",width(pdata1)[1]),collapse=""),length(savepdata1)))
      newset[match(labelsp1,savelabels),] <- pdata1
      #pdata1 <- pdata1[na.omit(match(savelabels,labelsp1)),]
      savepdata1 <- xscat(savepdata1,newset)
    }
  }
}


names(savepdata1) <- savelabels

allbin <- as.DNAbin(savepdata1)

allbin2 <-as.phyDat(allbin)

alldist <-dist.hamming(allbin2)

testnj <- NJ(alldist)
testup <- upgma(alldist)
plot.phylo(testnj, cex=.8)
plot.phylo(testup, cex=.8)

###subset out bad samples

trimsamples <- savepdata1[-c(48,16,5,47,45,14,26,42),]
trimbin <- as.DNAbin(trimsamples)
trimbin2 <-as.phyDat(trimbin)
trimdist <-dist.hamming(trimbin2)
trimup <- upgma(trimdist)

plot.phylo(trimup, cex=.4)

trimnj <- NJ(trimdist)
plot.phylo(trimnj, cex=.8)

write.dna(trimsamples, file="ambiguities_trimmed.fasta", format="fasta")

write.tree(trimup, file="consensus_trimmed_upgma.tree")







########
paste(rep("N",width(pdata1)[1]),collapse="")



p1data <- read.dna(file="", format = "fasta")

testcat <- xscat(pdata, pdata1)

consensus <- read.dna(file.choose(), format = "fasta")


matp1 <- as.matrix(distp1)

distp1 <-dist.dna(p1data)


distp4 <- dist.dna(p4data)
distp

disttest <-dist.dna(consensus)


testc <- as.dist(c(distp1, distp4))

testnj <- nj(distp1)
testnj4<- nj(distp4)



plot.phylo(testnj, cex=.8)


write.tree(nuc.nj.root, file="Figs/nj_1.tree")

plot.phylo(testnj4, cex=.8)


##how do i make a nj tree using all the primers and samples together?

write.tree(testnj, file="primer1_njtree_test")
pdf(paste("Figs/njContigs/nj_1","-",nrow(snp),".pdf",sep=""), width=10)
plot.phylo(testnj4, cex=.8)
dev.off()


#apply(tab[,3:7], 2, function(x) table(x))

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
