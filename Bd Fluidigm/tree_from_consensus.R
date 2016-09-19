
library(Biostrings)
library(GenomicAlignments)
library(ShortRead)
library(phangorn)

setwd("C:/Users/hji3/Desktop/Thesis Work/Fluidigm")
primerlist=read.table(file="Primer_list_details_half.csv",sep=",",header=T, stringsAsFactors=F)

primerlist = primerlist[-c(3,12,14),]
#rename the Mt primers to chrom number 70 because that is how it is in the ref seq
primerlist[92,4] <- 70
primerlist[93,4] <- 70
rownames(primerlist) <- 1:93

#setwd("C:/Users/hji3/Desktop/Thesis Work/Fluidigm/Results from run Aug 2015/Bdfluidigm.reduced.tar/Bdfluidigm.reduced/Bdfluidigm.reduced/consensus.split_amplicon_edited")
#files <- list.files(".")
#files <- grep("merged",files,value=TRUE)
setwd("C:/Users/hji3/Desktop/Thesis Work/Fluidigm/29 isolate fastas/cns_good")
ref <- readDNAStringSet(file="CJB4trim.fasta")
setwd("C:/Users/hji3/Desktop/Thesis Work/Fluidigm/Results from run Aug 2015/Bdfluidigm.reduced.tar/Bdfluidigm.reduced/Bdfluidigm.reduced/consensus.split_samples")

VACPS279 <- readDNAStringSet(file="Sample.V-ACPS279.fasta")
PH261 <- readDNAStringSet(file="Sample.V-PH261.fasta")
PMB16753 <- readDNAStringSet(file="Sample.V-PMB.16753.fasta")

JEL310 <- readDNAStringSet(file="Sample.Bd-JEL310.fasta")
UM142 <- readDNAStringSet(file="Sample.Bd-UM142.fasta")

Lb-Aber <- readDNAStringSet(file="Sample.Bd-Lb-Aber.fasta")
CJB4 <- readDNAStringSet(file="Sample.Bd-CJB4.fasta")

##great already has 93 primers (the 93 that worked)
names1 <- sapply(names(VACPS279), function(x) unlist(strsplit(x, split="merged"))[1])
names2 <- sapply(names1, function(x) unlist(strsplit(x, split="p"))[2])
names3 <- substr(names2, 1, nchar(names2)-1)
names(VACPS279) <- names3

##PH261
names1 <- sapply(names(PH261), function(x) unlist(strsplit(x, split="merged"))[1])
names2 <- sapply(names1, function(x) unlist(strsplit(x, split="p"))[2])
names3 <- substr(names2, 1, nchar(names2)-1)
names(PH261) <- names3

##PMB16753
names1 <- sapply(names(PMB16753), function(x) unlist(strsplit(x, split="merged"))[1])
names2 <- sapply(names1, function(x) unlist(strsplit(x, split="p"))[2])
names3 <- substr(names2, 1, nchar(names2)-1)
names(PMB16753) <- names3


##JEL310
names1 <- sapply(names(JEL310), function(x) unlist(strsplit(x, split="merged"))[1])
names2 <- sapply(names1, function(x) unlist(strsplit(x, split="p"))[2])
names3 <- substr(names2, 1, nchar(names2)-1)
names(JEL310) <- names3

##UM142
names1 <- sapply(names(UM142), function(x) unlist(strsplit(x, split="merged"))[1])
names2 <- sapply(names1, function(x) unlist(strsplit(x, split="p"))[2])
names3 <- substr(names2, 1, nchar(names2)-1)
names(UM142) <- names3
##CJB4
names1 <- sapply(names(CJB4), function(x) unlist(strsplit(x, split="merged"))[1])
names2 <- sapply(names1, function(x) unlist(strsplit(x, split="p"))[2])
names3 <- substr(names2, 1, nchar(names2)-1)
names(CJB4) <- names3

##remove if primer 21 or 25 are included
which(names(PH261)=="25") ##94
which(names(PH261)=="21")

PH261 <- PH261[-94]

which(names(PMB16753)=="25") ##94
which(names(PMB16753)=="21")

which(names(CJB4)=="25") ##94
which(names(CJB4)=="21")

#UM142 <- UM142[-94]
CJB4 <- CJB4[-c(94,95)]


#renames the names for DNAstringset to match those used for the 29 isolates
names(VACPS279) <- rownames(primerlist[match(names(VACPS279), primerlist$Primer_num),])

#PH261
names(PH261) <- rownames(primerlist[match(names(PH261), primerlist$Primer_num),])

#PMB16753
names(PMB16753) <- rownames(primerlist[match(names(PMB16753), primerlist$Primer_num),])

#UM142
names(UM142) <- rownames(primerlist[match(names(UM142), primerlist$Primer_num),])

#UM142
names(CJB4) <- rownames(primerlist[match(names(CJB4), primerlist$Primer_num),])

##makes a new set with the dimensions that should match the end product 
newset <- ref
##orders them
#newset <- VACPS279[match(names(ref),names(VACPS279))]
#newset <- PH261[match(names(ref),names(PH261))]
newset <- PMB16753[match(names(ref),names(PMB16753))]
newset <- CJB4[match(names(ref),names(CJB4))]

match(names(ref),names(UM142))[51]
match(names(ref),names(CJB4))[51]

#newset <- CJB4[match(names(ref),names(CJB4))[-51]]
##concatenates into a long sequence'
  copy <- newset
  for (x in 1:length(copy)){
    if(x==1){
      savedata <- copy[[x]]
    }
    savedata <- xscat(savedata,copy[[x]])
  }
savedata <- DNAStringSet(savedata[-(1:145)])

#write.dna(savedata, file="VACPS279_cat.fasta",format="fasta")

#write.dna(savedata, file="VPH261_cat.fasta",format="fasta")

write.dna(savedata, file="CJB4_cat.fasta",format="fasta")

width(ref)

width(newset)




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
i=1
for(i in 1:length(files)){
  pdata <- read.dna(file=files[i], format = "fasta")
  #pdata1 <- readDNAStringSet(file=files[i])
  if(!is.null(nrow(pdata))){
    pdata1 <- readDNAStringSet(file=files[i]) 
    labelsp1 <- sapply(names(pdata1), function(x) unlist(strsplit(x, split=":"))[1])
    if(i!=35){
      newset <- DNAStringSet(rep(paste(rep("N",width(pdata1)[1]),collapse=""),length(savepdata1)))
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

write.dna(trimsamples, file="consensus_trimmed.fasta", format="fasta")

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
