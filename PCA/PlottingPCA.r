###  R
# this script elaborates results from smartpca to creat PCA plots (Figure 2)
# first step is to modify the .evec file and add relevant information for each sample
# for example the population, a color and a symbol (pch) for plotting

# i import my eigenvectors modifies with relevant information
infoEigenMY<-read.table("EIGENVECTOROUT_modified.txt", sep = "\t", header=T)

pops<-table(infoEigenMY$population)
popolazionilist<-unlist(labels(pops))

# create a file that i will use as legend
aa<-c(1:length(popolazionilist))
coso<-matrix(NA,length(pops),4)
colnames(coso)<-c("population","group","color","pch")
for (i in 1:length(pops)){
  coso[i,1]<-popolazionilist[i]
  coso[i,2]<-as.character(infoEigenMY$group[which(infoEigenMY$population==popolazionilist[i] )[1]])
  coso[i,3]<-as.character(infoEigenMY$color[which(infoEigenMY$population==popolazionilist[i] )[1]])
  coso[i,4]<-as.character(infoEigenMY$pch[which(infoEigenMY$population==popolazionilist[i] )[1]])
}
coso2<-as.data.frame(coso)
cososort<-coso2[order(coso2$group, coso2$population),]
write.table(cososort,"listpopInfo.txt",sep = "\t")

poplistinfo<-read.table("listpopInfo.txt",sep = "\t", header=T) #this is a list of populations with respective group, color code and symbol. it will be useful for other analysis.

# now plot the PCA

layout(matrix(c(1,2), ncol=1), heights=c(5, 1)) # create a layout with two spaces, one is the PCA plot, the other is for the legend
par(xpd = T, mar = par()$mar + c(5,3,0.3,0.3))
plot(infoEigenMY$FIRST,infoEigenMY$SECOND, col=as.character(infoEigenMY$color), pch=infoEigenMY$pch, cex=1, cex.axis=0.6, cex.lab=0.6,xlab = "",ylab = "", main="PCA on my dataset")
title(ylab = "SECOND", line = 1.5, cex.lab = 0.5)
title(xlab = "FIRST", line = 1.5, cex.lab = 0.5)
text(infoEigenMY$FIRST,infoEigenMY$SECOND,labels = infoEigenMY$sample_ID, cex=0.2) # to plot minuscule labels and identify all the samples on the plot

plot.new()
par(mar=rep(0, 4))
legend("center", legend=poplistinfo$population, col=as.character(poplistinfo$color), pch=poplistinfo$pch, ncol=5, cex=0.8, y.intersp=1.4, x.intersp=1)

dev2bitmap("PCA_1and2_colorPOP_manual_Myset.pdf", type="pdfwrite") # it is a bit ugly. play with the layout and par margins to make it fit better.


# another simple example  of visualization in ggplot

library(ggplot2)
theme_set(theme_bw())  # pre-set the bw theme.
#infoEigen<-read.table("Eigen_smartPCA_withInfo.txt", sep = "\t", header=T)

pdf("PCA_1and2_colormyDatavsPublished.pdf")

gg<-ggplot(infoEigenMY, aes(x=FIRST, y=SECOND) ) + 
  geom_point(aes(col=set)) + 
  geom_text(aes(label=sample_ID),hjust=0, vjust=0, size=1.2)+
  labs(subtitle="PCA first two coordinates test ggplot", 
       y="SECOND", 
       x="FIRST", 
       title="Scatterplot", 
       caption = "populations")
#gg + scale_color_manual(values=groupcol)
gg
dev.off()



# ---------------------------------------------
### PLOT THE NEIGHBOR JOINING TREE FROM FST POPULATION DISTANCES (Figure S3)
# ---------------------------------------------

# specify fstonly=YES if you want only the fst in the parfile of smartpca.

FST<-read.table("FST_all") # your file from smartpca
FST<-FST[,-1]
rownames(FST)<-infoFST$population
colnames(FST)<-infoFST$population

m5<- as.dist(FST, diag=F, upper=F)

library(ape, pegas)
treeNJ<-nj(m5)
treeNJ$edge.length[treeNJ$edge.length < 0] = 0.002 # a little trick to adjust eventual "negative" branches

pdf("NJ_albero_allSNPS.pdf", useDingbats=FALSE)
plot.phylo(treeNJ, type="u", tip.col=as.character(infoFST$color), cex=0.4 )  # when you have a infoFST file with same order of the FST matrix and a color code assigned to each population
dev.off()


