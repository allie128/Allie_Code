###Pull out target regions for 29 sequenced isolates to compare Bd fluidigm data for Vance

setwd("C:/Users/hji3/Desktop/Thesis Work/Fluidigm")
library(Biostrings)
library(phangorn)

primerlist=read.table(file="Primer_list_details_half.csv",sep=",",header=T, stringsAsFactors=F)

#remove primers 5 (because we didnt use it), 21, and 25 (because they didnt work)
primerlist = primerlist[-c(3,12,14),]
#rename the Mt primers to chrom number 70 because that is how it is in the ref seq
primerlist[92,4] <- 70
primerlist[93,4] <- 70


#read in vcfs 
#library(VariantAnnotation)

#import the raw vcfs
#UM142vcf <- readVcf(file.choose(), "UM142")
#TST75vcf <- readVcf(file.choose(), "TST75")
#allvcf <- readVcf(file.choose(), "all")

#Philipp <- readVcf(file.choose(), "P")
#P.geno <- geno(Philipp)$GT


#all.geno <- geno(allvcf)$GT
#UM142.geno <- geno(UM142vcf)$GT
#TST75.geno <- geno(TST75vcf)$GT



####import the consensus fastas for each strain


##JEL423 reference

setwd("C:/Users/hji3/Desktop/Thesis Work/Fluidigm/29 isolate fastas/cns_good")


ref <- readDNAStringSet(file="JEL423_17Jan07.fa") 

#mt <- readDNAStringSet(file=file.choose()) 

#29 isolates
files <- list.files(".")
files <- grep("cns",files,value=TRUE)



##reads in all the files and concatenates them together
i=92

##takes each isolate fasta and pulls out the target sequences
##should end up with a list of 93 sequences
##takes a bit of time to run
for(i in 1:length(files)){
    isolate <- readDNAStringSet(file=files[i])
    isolatename = unlist(strsplit(files[i], split="_"))[1]
    isolatename = sub("-","_",isolatename)
    if(i==1){
      isolatenamelist <- isolatename
    }
      for(i in 1:nrow(primerlist)){
        chromnum <- as.numeric(primerlist[i,4])
        chrom <- isolate[[chromnum]]
        pref <- chrom[primerlist[i,6]:primerlist[i,7]]
        newset <- DNAStringSet(paste(pref))
          if(i==1){
            savedata=newset
            names(savedata)[i] <- i
          }
          savedata <- append(savedata, newset)
          names(savedata)[i] <- i
          if(i==93){
            savedata <- savedata[-1]
            names(savedata)[93] <- 93
          }}
      assign(paste(isolatename),savedata)
      isolatenamelist <- append(isolatenamelist, isolatename)
}
isolatenamelist <- isolatenamelist[-1]

#pulls out the target sequences from the ref for each primer

#for (i in 1:nrow(primerlist)){
  #if (primerlist[i,4]=="Mt"){
    #chrom <- mt[[1]]
    #pref <- chrom[primerlist[i,6]:primerlist[i,7]]
    #newset <- DNAStringSet(paste(pref))
    #}
      #chromnum <- as.numeric(primerlist[i,4])
      #chrom <- CJB4[[chromnum]]
      #pref <- chrom[primerlist[i,6]:primerlist[i,7]]
      #newset <- DNAStringSet(paste(pref))
    #if(i==1){
      #savedata=newset
    #}
  #savedata <- append(savedata, newset)
}
#gets rid of the extra row created
#savedata = savedata[-1]
##add the mt primers
#chrom <- mt[[1]]
#pref <- chrom[primerlist[92,6]:primerlist[92,7]]
#newset <- DNAStringSet(paste(pref))
#savedata <- append(savedata, newset)
#pref <- chrom[primerlist[93,6]:primerlist[93,7]]
#newset <- DNAStringSet(paste(pref))
#savedata <- append(savedata, newset)
#names(savedata) <- primerlist[,1]



##need to trim the primers off of the sequences
##also writes the data out as isolatetrim.fasta
i=2
for (i in 1:length(isolatenamelist)){
  copy <- get(isolatenamelist[i])
  if(i==1){
    trimisolatenamelist <- paste(isolatenamelist[i],"trim", sep="")
  }
  for (j in 1:nrow(primerlist)){
    nfor <- nchar(primerlist[j,2])
    nrev <- nchar(primerlist[j,3])
    copy[[j]] <- copy[[j]][-(1:nfor)]
    trimend <- width(copy)[j]-(nrev)+1
    copy[[j]] <- copy[[j]][-(trimend:width(copy)[1])]
    names(copy)[j] = j
  }
  assign(paste(isolatenamelist[i],"trim", sep=""),copy)
  trimisolatenamelist <- append(trimisolatenamelist, paste(isolatenamelist[i],"trim", sep=""))
  filename <- paste(isolatenamelist[i],"trim.fasta",sep="")
  write.dna(copy, file=filename,format="fasta")
}

trimisolatenamelist <- trimisolatenamelist[-1]

##need to concatenate all the sequences into one large sequence

##reads in all the files and concatenates them together
##should end up with a sequence length = 13676
i=1
x=1
for (i in 1:length(trimisolatenamelist)){
  copy <- get(trimisolatenamelist[i])
    for (x in 1:length(copy)){
      if(x==1){
        savedata <- copy[[x]]
      }
      savedata <- xscat(savedata,copy[[x]])
    }
  savedata <- DNAStringSet(savedata[-(1:145)])
  names(savedata) <- isolatenamelist[i]
  assign(paste(trimisolatenamelist[i],"_cat", sep=""),savedata)
  filename <- paste(trimisolatenamelist[i],"_cat.fasta",sep="")
  write.dna(savedata, file=filename,format="fasta")
}




