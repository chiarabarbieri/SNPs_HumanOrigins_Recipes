
### R
#### commands for IBD block analysis
# in this script, there will be several examples of plots aimed to show the IBD sharing patterns between populations.
# information for each population will be used.
# First, create a file that for each ID associates a population, some grouping categories (continent, macro greographic region, language family, etc), and geographic coordinate to plot on a map.
library(ggplot2)


info<-read.table("infosetIBD.txt", header=T) # this file includes a list of all the sample ID and relative population
poporder<-read.table("poplist.txt")[,1]  # this file includes one column with all the populations i want to consider in the analysis, one after the other, in the order i prefer
poporder<-as.vector(poporder)
infored<-info[which(info$population%in%poporder),] # this command will optionally subset a list of individuals in the case in the original files i have some populations i do not want to consider in the analysis
pops<-table(infored$populationOriginal)
pops<-pops[which(pops!=0)] 

# now workout the outputs from the beagle run, for each chromosome file

for (i in 1:22){
  ibd <- read.table(gzfile(paste("BeaglePhased",i,".ibd.gz",sep="", collapse="")))
ibd$source1<- info$population [match(ibd$V1, info$sample_ID)] # add a column with the population corresponding to the first ID of the pair
ibd$source2<- info$population [match(ibd$V3, info$sample_ID)] # add a column with the population corresponding to the second ID of the pair
ibd$length<-ibd$V7-ibd$V6 # add a column with the length of the segment in cM
colnames(ibd)[8]<-"LOD" # this is the column with the LOD score, which indicates how robust is the segment, proportional to how many times it is found in the beagle run
rownames(ibd) <- NULL
write.table(ibd,paste("ibd_adjusted_centimorgan_ch",i,".txt",sep="", collapse=""), sep="\t", row.names=FALSE, col.names=FALSE)
}


### now outside or R, in the shell, run the command
#  cat ibd_adjust* > IBD_adjustAllChromosomes.txt
#------------ back to R: 

ibd<-read.table("IBD_adjustAllChromosomes.txt",sep="\t", header=F)
colnames(ibd)<-c("firstID","firstHapIndex","secondID","secondHapIndex", "chromosome", "start","end","LOD","source1","source2","length")

# this little loops sorts the two population sources, so that an exchange between popA + popB will appear only in this order and not as popB + popA
coso<-cbind(as.character(ibd$source1),as.character(ibd$source2))
for(i in 1:nrow(ibd)){
  cosotemp<-coso[i,]
  cososort<-sort(cosotemp)
  ibd$source1[i]<-cososort[1]
  ibd$source2[i]<-cososort[2]
}
write.table(ibd,"IBD_adjustCentimorganAllChromosomesSort.txt",sep="\t")
ibd<-read.table("IBD_adjustCentimorganAllChromosomesSort.txt", header=T)

diffibd<-ibd[which(ibd$source1!=ibd$source2),] # ibd blocks shared by different pops (exclude within pop sharing)
 best<-diffibd[which(diffibd$LOD>10),] # use blocks with a LOD score higher than 10 to keep only the most significant ones

 ### SECTION 1
 # create a file to plot a symmetric matrix of exchange
 bestmirror<-best
 bestmirror$source1<-best$source2
 bestmirror$source2<-best$source1
 bestdouble<-rbind(best,bestmirror)
 
 # now plot as symmetric matrix 
 p<-ggplot(bestdouble, aes(x = source1, y = source2,size=length,color=LOD)) + 
   labs(x="source1", y="source2", title="IBD sharing LOD>10") +
   geom_point(shape=21) +
   scale_color_gradient(low="lightblue", high="red") +  # the color code corresponds to the LOD score
   theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                      axis.text.y=element_text(size=9),
                      plot.title=element_text(size=11),
                      legend.text=element_text(size=7)) +
   theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
   scale_x_discrete(limits=poporder)+
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

#______________________________________________________________________________________________________________

### SECTION 2
### IBD sharing patterns for each population
 
perpop<-matrix(NA,length(poporder),6)
colnames(perpop)<-c("population","samplesize","numbersharing","lengthsharing","totallenghtsharing","howmanypops")
perpop[,1]<-poporder
perpop[,2]<-pops[poporder]


for (i in 1:nrow(perpop)){
popp<-poporder[i]
temp<-best[c(which(best$source1%in%popp),which(best$source2%in%popp)),]
perpop[i,3]<-nrow(temp)
perpop[i,4]<-mean(temp$length)
popvarie<-c(temp$source1,temp$source2)
perpop[i,6]<-length(unique(popvarie))
perpop[i,5]<-sum(temp$length)
}

write.table(perpop,"popInfoIBDsharingPlay.txt",sep="\t")
perpop<-read.table("popInfoIBDsharingPlay.txt",header=T,sep="\t")

# now plots to display four values of interest

p1<-ggplot(data=perpop, aes(x=population, y=lengthsharing)) +
  geom_bar(stat="identity")+
   labs(x="populations", y="IBDblocks", title="average length IBD blocks") +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
  scale_x_discrete(limits=poporder)

p1
dev2bitmap("AveragelengthIBDblocksLOD10.pdf",type="pdfwrite")

p2<-ggplot(data=perpop, aes(x=population, y=howmanypops)) +
  geom_bar(stat="identity")+
   labs(x="populations", y="number pops", title="number of pops to share with") +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
  scale_x_discrete(limits=poporder)
p2

dev2bitmap("numberpopstosharewithLOD10.pdf",type="pdfwrite")

perpop$freq<-perpop$numbersharing/perpop$samplesize

p3<-ggplot(data=perpop, aes(x=population, y=freq)) +
  geom_bar(stat="identity")+
   labs(x="populations", y="number of sharing / pop size", title="freq of sharing") +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
  scale_x_discrete(limits=poporder)
p3

dev2bitmap("freqofSharingLOD10.pdf",type="pdfwrite")

p4<-ggplot(data=perpop, aes(x=population, y=totallenghtsharing)) +
  geom_bar(stat="identity")+
  labs(x="populations", y="IBDblocks", title="total length IBD blocks") +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
  scale_x_discrete(limits=poporder)
p4
dev2bitmap("totallengthSharingLOD10.pdf",type="pdfwrite")

# four plots together with multiplot 
# copy the function from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
pdf("combineInfoIBD_popsPlay.pdf")
multiplot(p1, p4,p2, p3, cols=2)  # function multiplot at the end of this script
dev.off()

#______________________________________________________________________________________________________________
### SECTION 3
#### plot the admixture with two target: in this case, Yoruba (for Africa) and Spain (for Europe)

target<-c("Yoruba", "Spanish")
AfricaEurope<-ibd[c(which(ibd$source1%in%target),which(ibd$source2%in%target)),]
AfricaEurope$SourceLocal<-NA
AfricaEurope$SourceOutside<-NA
write.table(AfricaEurope,"IBDAFricaEurope.txt", sep="\t")
AfricaEurope<-read.table("IBDAFricaEurope.txt", sep="\t", header=T)

Africa<-AfricaEurope[c(which(AfricaEurope$source1=="Yoruba"),which(AfricaEurope$source2=="Yoruba")),]
Africa$SourceLocal<-"Youruba"
Africa$SourceOutside<-"Youruba"
AfricaLess<-Africa[-which(as.character(Africa$source1)==as.character(Africa$source2)),]
for(i in 1:nrow(AfricaLess)){
  temp<-c(as.character(AfricaLess[i,9]),as.character(AfricaLess[i,10]))
  AfricaLess$SourceLocal[i]<-temp[-which(temp=="Yoruba")]
  AfricaLess$SourceOutside[i]<-temp[which(temp=="Yoruba")]
}
Africa<-rbind(AfricaLess, Africa[which(as.character(Africa$source1)==as.character(Africa$source2)),])

targetss<-"Spanish"
Spain<-AfricaEurope[c(which(AfricaEurope$source1==targetss),which(AfricaEurope$source2==targetss)),]
Spain$SourceLocal<-targetss
Spain$SourceOutside<-targetss
SpainLess<-Spain[-which(as.character(Spain$source1)==as.character(Spain$source2)),]
for(i in 1:nrow(SpainLess)){
  temp<-c(as.character(SpainLess[i,9]),as.character(SpainLess[i,10]))
  SpainLess$SourceLocal[i]<-temp[-which(temp==targetss)]
  SpainLess$SourceOutside[i]<-temp[which(temp==targetss)]
}
Spain<-rbind(SpainLess, Spain[which(as.character(Spain$source1)==as.character(Spain$source2)),])


mergeAfricaEuropeIBD<-rbind(Africa, Spain)
mergeAfricaEuropeIBDLOD10<-mergeAfricaEuropeIBD[which(mergeAfricaEuropeIBD$LOD>10),]

pdf("DistributionLengthIBDAfricaEurope_Play.pdf")
gg2<-ggplot(mergeAfricaEuropeIBDLOD10,aes(x=SourceLocal,y=length, fill=SourceOutside)) +
  geom_boxplot(alpha=0.7) +
  ggtitle("distribution length IBD blocks from Africa and Europe, LOD>10") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
  scale_x_discrete(limits=poporder)+
  scale_fill_brewer(palette = "Accent")
gg2
dev.off()


# ------------------------------------------------------------------------
### SECTION 4
### make a matrix of distance with the inverse of the number of haplotypes shared between populations
ibd<-ibd[which(ibd$LOD>10),]

matrixIBD<-matrix(NA, length(poporder),length(poporder), dimnames=list(poporder, poporder))
listpopInfo<-read.table("listpopInfo.txt", header=T)

for (i in 1:length(poporder)){
  for (k in 1:length(poporder)){
    pop.i=poporder[i]
    pop.k=poporder[k]
    temp<-rbind(ibd[which(ibd$source1==pop.i),],ibd[which(ibd$source2==pop.i),])
    if (pop.i==pop.k){
      matrixIBD[i,k]<- length(which(temp$source1==temp$source2))
    } else {
      tempp<-rbind(temp[which(temp$source1==pop.k),],temp[which(temp$source2==pop.k),])
      matrixIBD[i,k]<-nrow(tempp)
    }
  }
}

# make a matrix with the average length of blocks (it will come in use a bit later)
matrixIBDAverageLength<-matrix(NA, length(poporder),length(poporder), dimnames=list(poporder, poporder))

for (i in 1:length(poporder)){
  for (k in 1:length(poporder)){
    pop.i=poporder[i]
    pop.k=poporder[k]
    temp<-rbind(ibd[which(ibd$source1==pop.i),],ibd[which(ibd$source2==pop.i),])
    if (pop.i==pop.k){
      matrixIBDAverageLength[i,k]<- length(which(temp$source1==temp$source2))
    } else {
      tempp<-rbind(temp[which(temp$source1==pop.k),],temp[which(temp$source2==pop.k),])
      matrixIBDAverageLength[i,k]<-nrow(tempp)
    }
  }
}


inv.matrixIBD<-1/matrixIBD # the inverse of the matrix: more haplotypes shared, less distance between two pops
inv.matrixIBD[inv.matrixIBD==Inf] <- 2 # i set an arbitrary max distance of 2 for couple who do not share any IBD block
# now some population analysis from the distance matrix. plot population distance with MuldiDimensional Scaling (MDS)
# i load some extra info on each population from a file i elaborated

listpopInfo<-read.table("listpopInfoPlay.txy", header = T, sep = "\t")

library("MASS")

m5<- as.dist(inv.matrixIBD, diag=F, upper=F)

isoMDS(m5)-> distanti
pdf("MDS_2D_IBD_dissimilarity_play.pdf")
plot(distanti$points ,ylab="",xlab="", sub=paste("stress=",round(distanti$stress), collapse = ""), main="non metric MDS", pch=listpopInfo$pch, col=as.character(listpopInfo$color))
text(distanti$points,rownames(distanti$points),cex=0.3)
#dev2bitmap("MDS_2D_IBD_dissimilarity.pdf",type="pdfwrite")
dev.off()

#plot 3 dimension MDS
library(scatterplot3d)
scatterplot3d(wt,disp,mpg, pch=16, highlight.3d=TRUE,
              type="h", main="3D Scatterplot") 
isoMDS(m5,k=3)-> distanti
scatterplot3d(distanti$points, ylab="",xlab="", zlab="",type="h", sub=paste("stress=",round(distanti$stress)), main="non metric MDS", pch=listpopInfo$pch, color=as.character(listpopInfo$color))
s3d <- scatterplot3d(distanti$points, ylab="",xlab="", zlab="",type="h", sub=paste("stress=",round(distanti$stress)), main="non metric MDS", pch=listpopInfo$pch, color=as.character(listpopInfo$color))
text(s3d$xyz.convert(distanti$points),rownames(distanti$points),cex=0.3)


# a neighbor joining tree
library(ape, pegas)
treeNJ<-nj(m5)
treeNJ$edge.length[treeNJ$edge.length < 0] = 0.002

plot.phylo(treeNJ, type="u", tip.col=as.character(listpopInfo$color), cex=0.3 )


# ------------------------------------------------------------------------
### SECTION 5
# map visualization with ggplot and with "connecting flight" mode

# make working file for ggplot
library(reshape)
library(ggplot2)
meltIBD<-melt(matrixIBD)
colnames(meltIBD)<-c("source1", "source2", "n_sharing")

meltIBDaverage<-melt(matrixIBDAverageLength)
meltIBD$averageLength<-meltIBDaverage$value
write.table(meltIBD,"meltIBD_length_nofSharing.txt", sep="\t")

##### go on a MAP

infoRed<-listpopInfo[which(listpopInfo$lat!="NO"),] # i have some populations as reference for IBD sharing that i don't want to plot because they don't fit in the map
write.table(infoRed,"infoRed.toPlotPlay.txt", sep = "\t")
infoRed<-read.table("infoRed.toPlotPlay.txt", header=T)
library(rworldmap)
library(mapproj)
library(ggmap)
library(ggplot2)


map.world <- map_data(map="world")

library(ggrepel)

gg <- ggplot()
gg <- gg + theme()
gg <- gg + geom_map(data=map.world, map=map.world, aes(map_id=region, x=long, y=lat), fill="white", colour="black", size=0.15)
gg<- gg+coord_quickmap(ylim=c(-30,30), xlim=c(-111,-50)) # fix the limits with geo lon lat coordinates, according to the space you want to show
gg <- gg + geom_point(data=infoRed, aes(x=lon, y=lat, color=group)) +  # i color them for my column called "group"
  geom_label_repel(data=infoRed, aes(x=lon, y=lat,label=population), size=2, point.padding = 0.5,)
gg
ggsave("MapDotscolorelliPlay.pdf",gg, useDingbats=FALSE) # this is just a simple example to see if the population dots fall in the map area

testmaxminlimits<-meltIBD[-c(which(meltIBD$source1==meltIBD$source2),which(is.nan(meltIBD$averageLength))),]
limitslength<-c(min(testmaxminlimits$averageLength),max(testmaxminlimits$averageLength))
limitsSharing<-c(min(testmaxminlimits$n_sharing),max(testmaxminlimits$n_sharing))

# in this map i visualize all the populations that share with a given Target.
# the larger the dot, the more events of sharing
# color coded for average length of blocks (that's why i created the matrix of lenght of sharing just before)
# the loop will create and save a map for each population in my infoRed table.

poptarget<-infoRed$population
for (i in 1:length(poptarget)){
  target<-as.character(poptarget[i])
IBDtarget<-meltIBD[which(meltIBD$source2==target),]
IBDtargetRed<-IBDtarget[which(IBDtarget$source1%in%infoRed$population),]
IBDtargetRed$lon<-infoRed$lon
IBDtargetRed$lat<-infoRed$lat
IBDtargetRed<-IBDtargetRed[-which(IBDtargetRed$source1==target),]
targetCoordinate<-infoRed[which(infoRed$population==target),]
gg <- ggplot()
gg <- gg + theme()
gg <- gg + geom_map(data=map.world, map=map.world, aes(map_id=region, x=long, y=lat), fill="white", colour="black", size=0.15)
gg<- gg+coord_quickmap(ylim=c(-30,30), xlim=c(-111,-50))
gg <- gg + geom_point(data=IBDtargetRed, aes(x=lon, y=lat, color=averageLength, size=n_sharing)) +
scale_size_continuous(name="Number of sharing", range = c(0.5,7), limits=limitsSharing)+
  scale_color_gradient(low = "blue", high = "red",limits=limitslength) +
  ggtitle(paste("IBD sharing with ", target,collapse = ""))
gg<-gg + geom_point(data=targetCoordinate,aes(x=lon, y=lat, fill="black",shape=17, size=20)) +
                  scale_shape_identity()+ labs(fill="Target location")
ggsave(paste("MapDot_IBDplay",target,".pdf", sep = "", collapse = ""),gg, useDingbats=FALSE)
}


#_______________________________________________________
## now the map with network connection like in Barbieri et al. 2017 Sci Rep
# https://www.nature.com/articles/s41598-017-17728-w/figures/5

poptarget<-as.character(infoRed$population)
meltIBDnoIdenticalPairs<-meltIBD[-which(meltIBD$source1==meltIBD$source2),]
meltIBDForNetwork<-meltIBDnoIdenticalPairs[which(meltIBDnoIdenticalPairs$source1%in%poptarget),]
meltIBDForNetwork<-meltIBDForNetwork[which(meltIBDForNetwork$source2%in%poptarget),]
meltIBDForNetwork<-meltIBDForNetwork[which(meltIBDForNetwork$n_sharing!=0),]

library(maps)
library('geosphere')
library(rworldmap)
col.1 <- adjustcolor("darkorange", alpha=0.4)
col.2 <- adjustcolor("darkred", alpha=0.4)
edge.pal <- colorRampPalette(c(col.1, col.2), alpha = TRUE)
edge.col <- edge.pal(100)


pdf("mapSharingNetworkIBD_all_play.pdf")
map(database = "world", regions = ".",ylim=c(-30,30), xlim=c(-111,-50), col="grey90", fill=TRUE,  lwd=0.1)
points(x=infoRed$lon, y=infoRed$lat, pch=19,  cex=0.5, col="darkorange")

for(i in 1:nrow(meltIBDForNetwork))  {
  node1 <- infoRed[infoRed$population == as.character(meltIBDForNetwork[i,]$source1),]
  node2 <- infoRed[infoRed$population == as.character(meltIBDForNetwork[i,]$source2),]
  
  arc <- gcIntermediate( c(node1[1,]$lon, node1[1,]$lat), 
                         c(node2[1,]$lon, node2[1,]$lat), 
                         n=1, addStartEnd=TRUE )
  edge.ind <- round(100*meltIBDForNetwork[i,]$n_sharing / max(meltIBDForNetwork$n_sharing))
  
  lines(arc, col=edge.col[edge.ind], lwd=edge.ind/10)
}
text(x=infoRed$lon, y=infoRed$lat, labels=infoRed$population, pch=19,  cex=0.1, col="white") # i add tiny population text labels to cross check my results.
dev.off()


