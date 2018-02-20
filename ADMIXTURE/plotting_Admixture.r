R --
  # commands to elaborate admixture results and likelihood files, and plot

# cross validation error

read.table("CV.txt")->cv  

library(ggplot2)

p <- ggplot(cv, aes(x=V3, y=V4)) + 
  geom_boxplot()+ 
  ggtitle("values associated to each K")
p
pdf("distributionLikelihoodK_100run.pdf")
p
dev.off()


### control the likelihood of the 100 runs for each K

K_range<-c(4:10)  # select the number of K you run in Admixture

catll<-c()
for (i in K_range){
 ll<-read.table(paste("LL.K",i,".txt", sep="", collapse = ""))
 ll$K<-paste("ll",i,sep="", collapse = "")
 catll<-rbind(catll,ll)
}


# plot the distribution of likelihood associated to each K, to see if it looks regular or if it's bimodal or other problems
pdf("distributionLikelihoodK_perKSet1.pdf")
ggplot(catll,aes(x=V2,fill=K)) + geom_histogram(position="dodge") + ggtitle("likelihood associated to each K")
dev.off()

# find the run with the highest likelihood for each K
listina<-names(table(catll$K))
bestruns<-matrix(NA,length(K_range),2)
bestruns[,1]<-listina
for (i in 1:length(K_range)){
  bloc<-catll[which(catll$K==listina[i]),]
  massimo<-max(bloc$V2)
  bestruns[i,2]<-which(bloc$V2==massimo)
}

rownames(bestruns)<-gsub("ll","",bestruns[,1])

------------------------------
#it should look like  this
> bestruns
[,1]   [,2]
10 "ll10" "63"
4 "ll4"  "92"
5 "ll5"  "21"
6 "ll6"  "79"
7 "ll7"  "1" 
8 "ll8"  "54"
9 "ll9"  "71"
------------------------------

# prepare to plot: information to plot on the admixture bars, from your external info population file

  
MYINFO<-read.table("infoAdmixture.txt",header=T, as.is=T, sep="\t")
# this table must contain a row for each sample and the following columns: 
# oderAdmix	: is the order of the admixture file .Q
# sample_ID : the ID of each sample
# population	: the name of the population
# orderplot : the order you want to plot the samples (e.g. for continent, geography...)
------------------------------
  >head(MYINFO)
  oderAdmix	sample_ID	population	orderplot
  40	HGDP00920	Yoruba	1
  41	HGDP00924	Yoruba	2
  42	HGDP00925	Yoruba	3
  43	HGDP00926	Yoruba	4
------------------------------

table(MYINFO$population)->pops
namespop<-unique(MYINFO$population)

my.labels <- vector()   ## plotting pop labels instead of sample ids
for (k in 1:length(namespop)){
  paste("^",namespop[k],"$",sep="")->a
  length(grep(a, MYINFO$population)) -> my.labels[k]
}

labels.coords<-vector()  ### where to plot pop labels
labels.coords[1]<-my.labels[1]/2
for (i in 1:(length(my.labels)-1)) {
  labels.coords[i]+(my.labels[i]/2+my.labels[i+1]/2)->labels.coords[i+1]
}
z<-vector()
z[1]<-my.labels[1]
for (i in 1:(length(my.labels)-1)){
  z[i]+my.labels[i+1]->z[i+1]
}

# now plot for each K
K<-4 # chose the K to plot. Start with the lowest.

choiceRun<-bestruns[which(rownames(bestruns)==K),]
valuesToplot<-read.table(paste("K",K,".Run", choiceRun[2], ".Q", sep="", collapse = ""))

tblSort<-tbl[MYINFO$oderAdmix,]

pdf(paste("AdmixtureForK",K,".pdf", sep="", collapse = ""),pointsize=8, height=3.5)

barplot(t(as.matrix(tblSort)), col=c("brown","green","purple","orange"), axisnames=F, axes=F, space=0,border=NA)
axis(3, labels = FALSE, tick=F)
for (i in c(0,z)){
  lines(x=c(i,i),y=c(0,3), lwd=0.7, col="white")
}
text(labels.coords, par("usr")[3] + 1.03 , srt = 45, adj = 0, labels = namespop, cex=0.7, xpd = TRUE)
dev.off()

# and repeat for K=5, K=6 etc.
# This plotting system has a problem: for each new K i have to add a color manually, and the order of the colors for ancestry block is shuffled in every K
# so i have to manually sort the colors to keep the color scheme in all K.
# other alternatives ways to plot i tried did not satisfy me yet. 

# use colorbrewer or similar to find a  combination of colors graphically satisfying. 