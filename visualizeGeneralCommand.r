
####### in R

### search for relatives with plink --genome
filone<-read.table("plink.genome", header=T)
pdf("distributionPIHAT.pdf")
hist(filone$PI_HAT, breaks=100)
dev.off()

filone[which(filone$PI_HAT>0.24),]->Relat  # set the threshold to spot relatives. Second degree relative correspond to a cousin.

write.table(Relat, "couplesCousins.txt", sep="\t") # open this table to manually screen potential relative pairs

library(ggplot2)

# plot the distribution of pi hat for family groups
# remember when you replaced the .ind file with the names of the populations, before running plink? 
# FID1:	Family ID for first sample should correspond to the population

p <- ggplot(filone, aes(FID1,PI_HAT))
p + geom_boxplot()
dev2bitmap("boxplotDistribuPiHAT.pdf", type="pdfwrite")

## now something simila for the F value, a measure of heterozygozity and consanguineity
het<-read.table("plink.het", header=T)

p <- ggplot(het, aes(FID,F))
p + geom_boxplot()
dev2bitmap("boxplotDistribuF_heterozigosity.pdf", type="pdfwrite")


