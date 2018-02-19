
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

#__________________________________________
### plot Runs of Homozygosity (ROH)

# use an external file with information on population, group, plotting color etc for each individual
roh<-read.table("yourfile_ROH.hom.indiv", header=T, as.is=T)
myinfo<-read.table("infoData.txt",header=T)
roh<-roh[which(roh$IID%in%myinfo$sample_ID),]
rownames(myinfo)<-myinfo$sample_ID
myinfoSort<-myinfo[roh$IID,]

# an external file with characteristic same as above but for each population (is for the legend)
popinfo<-read.table("listpopInfoPlay.txt", header=T)

pdf(file="ROH_plot_numberofROH_length", pointsize=12,width=12)
plot(x=roh[,5]/1000,y=roh[,4],col=as.character(myinfoSort$color),pch=as.numeric(myinfoSort$pch),xlab="total length in ROH (Mb)", ylab="number of ROH")
legend("bottomright", legend=popinfo$population,  pch=as.numeric(popinfo$pch), col=as.character(popinfo$color), text.col="gray", ncol=5,bty="n", cex=0.5)
dev.off()

## another kind of visualization with ggplot
mergina<-cbind(roh,myinfoSort)

library(ggplot2)

pdf(file="ROH_test.pdf", pointsize=12,width=12)

gg <- ggplot(mergina, aes(x=population, y=KB, fill=group) ) + 
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Distribution of ROH lengths in KB per population") +
  scale_fill_manual(labels = c("group 1", "group 2", "group 3"), values = c("brown","orange", "red"))
gg
dev.off()

