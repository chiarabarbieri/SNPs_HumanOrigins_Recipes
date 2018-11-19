### R
#### commands for IBD block analysis
# in this script, there will be several examples of plots aimed to show the IBD sharing patterns between populations.
# information for each population will be used.
# First, create a file that for each ID associates a population, some grouping categories (continent, macro greographic region, language family, etc), and geographic coordinate to plot on a map.


ibd<-read.table("all.refinedIBD.Merged", as.is=T)
colnames(ibd)<-c("firstID","firstHapIndex","secondID","secondHapIndex", "chromosome", "start","end","LOD","length")


infocomplete<-read.csv("infoForPopulation.csv", header=T, as.is=T , comment.char = "", fill=T) #info file each line one population
infoID<-read.csv("InfoForIndividual.csv",header=T, as.is=T , comment.char = "", fill=T).  #info file each line one individual

minimuminfo<-infoID[,c(2,1)]   # select the columns which have the sample ID name and the relative population
colnames(minimuminfo)[1]<-"firstID"
ibdmatch<-merge(ibd,minimuminfo,all.x=TRUE)   # associate the population source for the first sample ID of the couple
colnames(ibdmatch)[10]<-"source1"
colnames(minimuminfo)[1]<-"secondID"
ibdmatch2<-merge(ibdmatch,minimuminfo,all.x=TRUE). # associate the population source for the second sample ID of the couple
colnames(ibdmatch2)[11]<-"source2". 

write.table(ibdmatch2, "refinedIBD_merged_withinfo.txt", row.names = F, quote=F)
ibd<-ibdmatch2

### table sharing per population

poporder<-unique(infoID$PopName)
pops<-table(infoID$PopName)

perpop<-matrix(NA,length(poporder),11)
colnames(perpop)<-c("population","samplesize","numberSharingTot","numbersharingWithin","numberSharingOut","FreqSharingTot","FreqsharingWithin","FreqSharingOut","Mean_lengthsharingWithin","totallenghtsharing","howmanypops")
perpop[,1]<-poporder
perpop[,2]<-pops[poporder]

for (i in 1:nrow(perpop)){
popp<-poporder[i]
within<-which(ibd$source1%in%popp & ibd$source2%in%popp)
tempWithin<-ibd[within,]
tempTOT<-ibd[union(which(ibd$source1%in%popp),which(ibd$source2%in%popp)),]
tempOut<-tempTOT[-which(tempTOT$source1 == tempTOT$source2),]
perpop[i,3]<-nrow(tempTOT)
perpop[i,4]<-nrow(tempWithin)
perpop[i,5]<-nrow(tempOut)
perpop[i,6]<-nrow(tempTOT)/as.numeric(perpop[i,2])
perpop[i,7]<-nrow(tempWithin)/as.numeric(perpop[i,2])
perpop[i,8]<-nrow(tempOut)/as.numeric(perpop[i,2])
perpop[i,9]<-mean(tempTOT$length)
popvarie<-c(tempTOT$source1,tempTOT$source2)
perpop[i,10]<-sum(tempTOT$length)
perpop[i,11]<-length(unique(popvarie))
}
perpop2<-merge(perpop, infoID[,c(1,8,11,12)], by.x=perpop[,1], by.y=infoID[,1])

write.table(perpop,"popInfoIBDsharing.txt",sep="\t", row.names = F, quote=F)
perpop<-read.table("popInfoIBDsharing.txt",sep="\t",header=T, as.is=T,comment.char = "", fill=T, quote="")

#------------------------------------------------------------------
 ### SECTION 1 visualize exchange between populations as number of events
 #------------------------------------------------------------------

 # create a file to plot a symmetric matrix of exchange
 bestmirror<-ibd
 bestmirror$source1<-ibd$source2
 bestmirror$source2<-ibd$source1
 bestdouble<-rbind(ibd,bestmirror)
 
 # now plot as symmetric matrix 
 p<-ggplot(bestdouble, aes(x = source1, y = source2,size=length,color=LOD)) + 
   labs(x="source1", y="source2", title="IBD sharing") +
   geom_point(shape=21) +
   scale_color_gradient(low="lightblue", high="red") +  # the color code corresponds to the LOD score
   theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                      axis.text.y=element_text(size=9),
                      plot.title=element_text(size=11),
                      legend.text=element_text(size=7)) +
   theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
   scale_x_discrete(limits=poporder)+       #the poporder can be a different population order that you want to display in the plot
   scale_y_discrete(limits=poporder)
 p
 ggsave("IBDsharingblocksSymmetricPlay.pdf", useDingbats=FALSE)
 
# the plot above did not have the diagonal (within pop sharing)
 # this plot shows the sharing of segments within each population, color code for LOD, length on y axis
 ibd<-ibd[which(ibd$LOD>10),]
 
 diagonal<-ibd[which(ibd$source1==ibd$source2),]
pintern<-ggplot(diagonal, aes(x = source1,y=length,color=LOD)) + 
  labs(x="source1", title="IBD sharing within population LOD>10") +
  geom_point(shape=21) +
  scale_color_gradient(low="lightblue", high="red") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11),
                     legend.text=element_text(size=7)) +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
  scale_x_discrete(limits=poporder)
pintern
ggsave("IBDsharingblocksWithinPopPlay.pdf", useDingbats=FALSE)

#------------------------------------------------------------------
Matrix of exchange between populations
#------------------------------------------------------------------

# matrix with the total number of shared blocks
matrixIBD<-matrix(NA, length(poporder),length(poporder), dimnames=list(poporder, poporder))

for (i in 1:length(poporder)){
  for (k in 1:length(poporder)){
    pop.i=poporder[i]
    pop.k=poporder[k]
    temp<-ibd[union(which(ibd$source1==pop.i),which(ibd$source2==pop.i)),]
    if (pop.i==pop.k){
      matrixIBD[i,k]<- length(which(temp$source1==temp$source2))
    } else {
      tempp<-rbind(temp[which(temp$source1==pop.k),],temp[which(temp$source2==pop.k),])
      matrixIBD[i,k]<-nrow(tempp)
    }
  }
}
write.table(matrixIBD,"matrix_refinedIBD_merge_sharing.txt", sep="\t")

# make a matrix with the average length of blocks
matrixIBDAverageLength<-matrix(NA, length(poporder),length(poporder), dimnames=list(poporder, poporder))

for (i in 1:length(poporder)){
  for (k in 1:length(poporder)){
    pop.i=poporder[i]
    pop.k=poporder[k]
    temp<-ibd[union(which(ibd$source1==pop.i),which(ibd$source2==pop.i)),]
      if (pop.i==pop.k){
        temp2<-temp[(which(temp$source1==temp$source2)),]
      matrixIBDAverageLength[i,k]<- mean(temp2$length)
    } else {
      tempp<-rbind(temp[which(temp$source1==pop.k),],temp[which(temp$source2==pop.k),])
      matrixIBDAverageLength[i,k]<-mean(tempp$length)
    }
  }
}
write.table(matrixIBDAverageLength,"matrix_IBD_averageLength.txt", sep="\t")


# make a matrix with the TOTAL length of blocks
matrixIBDTotLength<-matrix(NA, length(poporder),length(poporder), dimnames=list(poporder, poporder))

for (i in 1:length(poporder)){
  for (k in 1:length(poporder)){
    pop.i=poporder[i]
    pop.k=poporder[k]
    temp<-ibd[union(which(ibd$source1==pop.i),which(ibd$source2==pop.i)),]
    if (pop.i==pop.k){
      temp2<-temp[(which(temp$source1==temp$source2)),]
      matrixIBDTotLength[i,k]<- sum(temp2$length)
    } else {
      tempp<-rbind(temp[which(temp$source1==pop.k),],temp[which(temp$source2==pop.k),])
      matrixIBDTotLength[i,k]<-sum(tempp$length)
    }
  }
}
write.table(matrixIBDTotLength,"matrix_IBD_totalLength.txt", sep="\t")

pops<-table(infoID$PopName)
#pops<-pops[which(pops>0)]
pops<-pops[rownames(matrixIBD)]

#adjust for population size
matrixIBDadjustpopsize<-matrixIBD
for (i in 1:length(poporder)){
  for (k in 1:length(poporder)){
    pop.i=poporder[i]
    pop.k=poporder[k]
      matrixIBDadjustpopsize[i,k]<- matrixIBD[i,k]/(pops[i]*pops[k])
  }
}
write.table(matrixIBDadjustpopsize,"matrix_IBDsharingAdjustPopSize.txt", sep="\t")

matrixIBDadjustpopsizelength<-matrixIBDTotLength
for (i in 1:length(poporder)){
  for (k in 1:length(poporder)){
    pop.i=poporder[i]
    pop.k=poporder[k]
    matrixIBDadjustpopsizelength[i,k]<- matrixIBDTotLength[i,k]/(pops[i]*pops[k])
  }
}
write.table(matrixIBDadjustpopsizelength,"matrix_IBDsharingAdjustPopSizelength.txt", sep="\t")

# ---------------
# Melt everything for ggplot

library(reshape)
library(ggplot2)
meltIBD<-melt(matrixIBD)
colnames(meltIBD)<-c("source1", "source2", "n_sharing")

meltIBDaverage<-melt(matrixIBDAverageLength)
meltIBD$averageLength<-meltIBDaverage$value
meltIBDlength<-melt(matrixIBDTotLength)
meltIBD$totalLength<-meltIBDlength$value
meltIBDadjuxt<-melt(matrixIBDadjustpopsize)
meltIBD$sharingadjust<-meltIBDadjuxt$value
meltIBDlengthadjuxt<-melt(matrixIBDadjustpopsizelength)
meltIBD$lengthadjust<-meltIBDlengthadjuxt$value


meltIBD2<-meltIBD[-(which(meltIBD$source1==meltIBD$source2)),] #exclude same pop sharing

# you can filter for more significant pairs, like pairs that share more than once, or more than the median
meltIBD22<-meltIBD2[which(meltIBD2$n_sharing!=0),]
meltIBD22<-meltIBD2[which(meltIBD2$n_sharing>1),]
meltIBD22<-meltIBD22[which(meltIBD22$sharingadjust>median(meltIBD22$sharingadjust)),]

gg<-ggplot(meltIBD22,aes(x=source1, y=source2, fill=sharingadjust, size=lengthadjust))+
  geom_point(shape=21)+
  # scale_fill_gradient(low="lightblue", high="red") +
  # scale_fill_brewer(palette ="Spectral") +
  scale_fill_gradientn(colours = terrain.colors(7))+
  # scale_fill_distiller(palette = "RdPu")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1,colour=as.character(info$color))) +
  scale_x_discrete(limits=poporder)+
  scale_y_discrete(limits=poporder)
gg

ggsave("matrixAdjustSize_length_sharing_IBDmerged.pdf", useDingbats=FALSE) # Figure 4A

# ------------------------------------------------------------------------
# geographic distances
# ------------------------------------------------------------------------

infoRed<-info[which(info$lat!="NO"),]

infoIDred<-infoID[which(infoID$population%in%poporder),]
pops<-table(infoIDred$population)
pops<-pops[poporder]


library(fields)
lista<-(cbind(as.numeric(infoRed$lon), as.numeric(infoRed$lat)))
MatrixGeo<-rdist.earth (lista, miles=FALSE)
rownames(MatrixGeo)<-infoRed$population
colnames(MatrixGeo)<-infoRed$population

geomelt<-melt(MatrixGeo)



# ------------------------------------------------------------------------
### SECTION 5
# map visualization with ggplot and with "connecting flight" mode.  - Figure 4B
## now the map with network connection like in Barbieri et al. 2017 Sci Rep
# https://www.nature.com/articles/s41598-017-17728-w/figures/5

meltIBD3<-meltIBD2[which(meltIBD2$source1%in%geomelt$X1),]
meltIBD3<-meltIBD3[which(meltIBD3$source2%in%geomelt$X1),]
#meltIBD3<-meltIBD3[which(meltIBD3$source2!=meltIBD3$source1),]

#geomelt<-geomelt[which(geomelt$X1!=geomelt$X2),]
meltIBD3$geodist<-as.character(geomelt$value)
meltIBD3<-meltIBD3[which(meltIBD3$n_sharing!=0),]

meltIBD3$lengthGeo<-meltIBD3$lengthadjust*as.numeric(meltIBD3$geodist)
meltIBD3$sharingGeo<-meltIBD3$sharingadjust*as.numeric(meltIBD3$geodist)

meltIBD33<-meltIBD3[which(meltIBD3$n_sharing>1),]
meltIBD33<-meltIBD33[which(meltIBD33$sharingadjust>median(meltIBD33$sharingadjust)),]

library(maps)
library('geosphere')

col.1 <- adjustcolor("darkorange", alpha=0.4)
col.2 <- adjustcolor("darkred", alpha=0.4)
edge.pal <- colorRampPalette(c(col.1, col.2), alpha = TRUE)
edge.col <- edge.pal(100)

pdf("mapSharingNetworkRefinedIBD_.pdf")
  
  map(database = "world", regions = ".",ylim=c(-30,30), xlim=c(-111,-50), col="grey90", fill=TRUE,  lwd=0.1)
  points(x=infoRed$lon, y=infoRed$lat, pch=19,  cex=0.5, col="darkorange")
  
  for(i in 1:nrow(meltIBD33))  {
    node1 <- infoRed[infoRed$population == as.character(meltIBD33[i,]$source1),]
    node2 <- infoRed[infoRed$population == as.character(meltIBD33[i,]$source2),]
    
    arc <- gcIntermediate(as.numeric(c(node1[1,]$lon, node1[1,]$lat)), 
                          as.numeric(c(node2[1,]$lon, node2[1,]$lat)), 
                          n=1, addStartEnd=TRUE )
    edge.ind <- round(100*meltIBD33[i,]$sharingadjust / max(meltIBD33$sharingadjust))
    
    lines(arc, col=edge.col[edge.ind], lwd=edge.ind/10)
  }
  text(x=as.numeric(infoRed$lon), y=as.numeric(infoRed$lat), labels=infoRed$shortPop,  cex=0.2, col="white")
  
  dev.off()

}

