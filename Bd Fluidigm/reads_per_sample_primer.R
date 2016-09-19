###########Bd Fluidigm Data Analysis#####################

setwd("C:/Users/hji3/Desktop/Thesis Work/Fluidigm/Results from run Aug 2015/Bdfluidigm.reduced.tar/Bdfluidigm.reduced/Bdfluidigm.reduced")

source("https://bioconductor.org/biocLite.R")
biocLite("ShortRead")
biocLite("GenomicAlignments")
library(GenomicAlignments)
library(ShortRead)

fqp1 <- readFasta(file.choose())


##To make plot for read count per sample
##read in read count file
count <- read.table(file.choose(),header=T, stringsAsFactors=F, row.names=1)

boxplot(count)

p1 <- count[,1]

readsperprimer <- (sapply(1:ncol(count), function(x) mean(count[,x])))

readspersample <- (sapply(1:nrow(count), function(x) mean(as.numeric(count[x,]))))

plot(readsperprimer, ylab="Number of Amplicons", main="Amplicons per Primer")
plot(readspersample, ylab="Number of Amplicons", main="Amplicons per Sample", col=colors)



samplenames <- read.table(file.choose(), stringsAsFactors=F, sep="\n")
##assumes that the barcodes in the Identified_Barcodes file are in order
rps <- cbind(samplenames, readspersample)
colnames(readspersample) <- samplenames
colors <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6)
barplot(rps, ylab="Number of Amplicons", main="Amplicons per Sample")
#figure out how to tilt the bar titles then

barplot(height=rps[,2], ylab="Number of Amplicons", main="Amplicons per Sample", names.arg=rps[,1], cex.names=.1, col=colors)

write.table(rps,file="reads_per_sample.tab",sep="\t",quote=F,row.names=F)

##compare ze to amplicons

loaddata <- read.table(file.choose(), stringsAsFactors=F, header=T, sep=",")
compare <- data.frame(loaddata[4],loaddata[5])
colors2 <- c(1,1,1,1,1,1,1,1,1,1,1,1,4,4,4,4,5,5,5,5,5,5,2,2,2,2,2,2)
plot(compare, main='Reads per Sample vs ZE', col=colors2, pch=16)
plot(loaddata[4], loaddata[5], ylab="Number of Amplicons", xlab="ZE")

#get rid of outliers
compare2 <- compare[-1,]
colors3 <- c(1,1,1,1,1,1,1,1,1,1,1,4,4,4,4,5,5,5,5,5,5,2,2,2,2,2,2)
plot(compare2, main='Reads per Sample vs ZE', col=colors3, pch=16)
which(compare$ZE > 3*sd(compare$ZE))

##parse by sample type

vance <- compare[1:12,]
plot(vance, main='Samples from Vance', col=1, pch=16)
#remove outlier
plot(vance[-c(1,4),], main='Samples from Vance', col=1, pch=16)
##
reinf <- compare[13:16,]
plot(reinf, main='Rana Reinfection Experiment', col=4, pch=16)
##
ranawild <- compare[17:22,]
plot(ranawild, main='Rana Swabs from Field', col=5, pch=16)
##
ranainf <- compare[23:28,]
plot(ranainf, main='Rana Infection Experiment', col=2, pch=16)





