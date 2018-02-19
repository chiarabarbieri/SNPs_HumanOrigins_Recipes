

### prepares the outfile to manipulate in R

maldersuccess<-read.table("outputSuccessMalder.cat.txt")
colnames(maldersuccess)<-c("test", "status",     "p-value", "test_pop","ref_A","ref_B", "Zscore",  "divergenceTime", "null", "interval")

library(ggplot2)
ThisYear<-2018  # to have the admixture times in calendar years
generationyeras<-30  # chose your generation time

maldersuccess$AgeCalendarYear<-ThisYear-(maldersuccess$divergenceTime*generationyeras)
maldersuccess$Min<-ThisYear-((maldersuccess$divergenceTime+maldersuccess$interval)*generationyeras)
maldersuccess$Max<-ThisYear-((maldersuccess$divergenceTime-maldersuccess$interval)*generationyeras)
maldersuccess$couple<-paste(maldersuccess$ref_A,maldersuccess$ref_B)

poporder<-c()
for (i in 1:length(pops)){
  poporder[i]<-paste("POP",i,collapse = "", sep="")
}
# or just list the populations in any order you want to plot them

pdf("divergenceAge_EuropeAfrica_Malder.pdf")
gg<-ggplot(maldersuccess,aes(x=AgeCalendarYear,y=test_pop, group=couple,color= couple))+
            geom_errorbarh(aes(xmin=Min,xmax=Max), alpha=0.2) +
            geom_point(na.rm=TRUE)+
            ggtitle("Admixture Times from Africa and Europe") +
            theme_bw() +
            #scale_fill_brewer(palette = "Accent")+
            scale_color_manual(labels = c("Europe", "Africa"), values = c("darkgreen","purple"))+
            scale_y_discrete(limits=poporder)     # this will plot all the Target admixed in the right order on the Y axis       
gg
dev.off()

# i want to exclude the admixture not supported by p value <0.01
maldersignificant<-maldersuccess[which(maldersuccess[,3]<0.001),]

pdf("divergenceAge_EuropeAfrica_Xavante_MalderSignificant0001.pdf")
gg<-ggplot(maldersignificant,aes(x=AgeCalendarYear,y=test_pop,color= couple))+
            geom_errorbarh(aes(xmin=Min,xmax=Max), alpha=0.2) +
            geom_point(na.rm=TRUE)+
            ggtitle("Admixture Times from Africa and Europe") +
            theme_bw() +
            scale_fill_brewer(palette = "Accent")+
              scale_y_discrete(limits=poporder)            
gg
dev.off()

# now color for the Z-score 
pdf("divergenceAge_EuropeAfrica_Xavante_MalderSignificant0001_zscore.pdf")
gg<-ggplot(maldersignificant,aes(x=AgeCalendarYear,y=test_pop,color= Zscore))+
            geom_errorbarh(aes(xmin=Min,xmax=Max), alpha=0.2) +
            geom_point(na.rm=TRUE)+
            scale_colour_gradientn(colours = c("darkred","red","orange"))+
             ggtitle("Admixture Times from Africa and Europe") +
            theme_bw() +
              scale_y_discrete(limits=poporder)            
gg
dev.off()




# _________ plot the admixture LD decay

LD<-read.table("LD_decaycurveALLfilter.txt",header=T)
LD<-LD[-which(LD$d=="Inf"),] # exclude values of d that are marked as Infinite

LD<-LD[which(LD$Target%in%c("Bolivian", "Zapotec", "Quechua")),]

library(ggplot2)
poporder1<-c("Bolivian", "Zapotec", "Quechua")

gg<-ggplot(LD,aes(x=d, y=value, color=variable))+
geom_point(aes(x=d,y = Spanish.Ecuador, col="Spanish.Ecuador"), shape=3, size=0.02)+
geom_point(aes(x=d,y = Yoruba.Ecuador, col="Yoruba.Ecuador"), shape=3, size=0.02)+
  scale_color_manual(labels = c("Europe", "Africa"), values = c("darkgreen","purple"))+
 ggtitle("Admixture LD decay from Africa and Europe") +
  theme_bw() +
facet_wrap(~Target, ncol=1)
gg

ggsave("LDdecay_2sources.pdf",useDingbats=FALSE)

