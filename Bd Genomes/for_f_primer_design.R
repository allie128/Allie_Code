

setwd("C:/Users/hji3/Desktop/Thesis Work/Panama Isolate Genome Project")

tab=read.table(file="snps_6strains_oldformat_altrem.tab",sep="\t",header=T, stringsAsFactors=F)
taball=read.table(file="snps_oldandnewdata.tab",sep="\t",header=T, stringsAsFactors=F)

apply(tab[,3:7], 2, function(x) table(x))

chrom1idx <- which(tab$CHROM=='bden_JEL423_supercont1.1')
chrom1 <- tab[chrom1idx,]
which(chrom1$Campana!=chrom1$JEL410 & chrom1$JEL410!=chrom1$JEL412 & chrom1$JEL412!=chrom1$JEL413 & chrom1$JEL413!=chrom1$Rio_Maria)
#9576 14444 14445 22016 27181 31936 32334 39555 42993 58264 5887

which(tab$Campana!=tab$JEL410 & tab$JEL410!=tab$JEL412 & tab$JEL412!=tab$JEL413 & tab$JEL413!=tab$Rio_Maria & tab$Rio_Maria!=tab$Sora)


which(tab$Rio_Maria!=tab$JEL410 & tab$JEL410!=tab$JEL412 & tab$JEL412!=tab$JEL413 & tab$Rio_Maria!=tab$JEL412)



chrom2idx <- which(tab$CHROM=='bden_JEL423_supercont1.2')
chrom3idx <- which(tab$CHROM=='bden_JEL423_supercont1.3')
chrom4idx <- which(tab$CHROM=='bden_JEL423_supercont1.4')
chrom5idx <- which(tab$CHROM=='bden_JEL423_supercont1.5')


#find regions with lots of differences between isolates

#remove rows with na values
tabnona <- na.omit(tab)

rowsum <- sum(tabnona[,3:8])

##
testrs <- rowSums(tabnona[,3:8])


threes <- which(testrs==3)
sumthree <- tabnona[threes,]

##JEL423 reference

setwd("C:/Users/hji3/Desktop/Thesis Work/Fluidigm/29 isolate fastas/cns_good")


ref <- readDNAStringSet(file="JEL423_17Jan07.fa") 

ref[[5]][208652:208662]

newprimer1 <- ref[[5]][208452:209052]

newprimer1 <- DNAStringSet(newprimer1)
names(newprimer1) <- "Panama_primer_5"
newprimer2 <- ref[[57]][3093:3583]
newprimer3 <- ref[[39]][663:1361]

newprimer3 <- DNAStringSet(newprimer3)
names(newprimer3) <- "Panama_primer_39"

write.dna(newprimer3, file="panama_primer_39",format="fasta")

