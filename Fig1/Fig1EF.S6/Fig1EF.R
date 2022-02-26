#################################################################################################
#Description: generate Fig1E-F
#Author: Sheng Guo PhD  guosheng@crownbio.com
#Date:   Aug-2021
#Apache License
#Version 2.0, January 2004
#https://opensource.org/licenses/Apache-2.0
#################################################################################################
library(ggpubr)

MINSAMPLES = 20 		#mimimal number of samples for a cancer to be included
MINPDX = 10				#mimimal number of PDX for a cancer to be included
MAX_PASSAGE= 10 		#max passage for a samples to be inclued 
MIN_SAMPLE2PLOT=30		#for individual cancer, minimal number of samples to be plotted
MIN_PASSAGE2PLOT=4		#for individual cancer, minimal number of passage difference to be plotted
MAX_SAMPLE_PER_PDX = 20


nod<-read.csv('nod.csv')

#--remove Cancer with few samples--#
a=as.data.frame(table(nod$Cancer)<MINSAMPLES)
a=cbind(a,a)
colnames(a) = c('t1','t2')
toRemove=rownames(a[a$t1==TRUE,])
nod=nod[!(nod$Cancer %in% toRemove),]

#--remove Cancer with few unique PDX models--#
Cancers = unique( nod$Cancer)
CancerKeep = c()
for(i in 1:length(Cancers)){
	nod2 = nod[nod$Cancer== Cancers[i],]
	if(length(unique(nod2$DuplicateModel))>MINPDX){
		CancerKeep = c(CancerKeep, Cancers[i])
	}
}
length(CancerKeep)



nod=nod[(nod$Cancer %in% CancerKeep),]

#--for each PDX, randomly keep at most 5 samples--#
Cancers = unique( nod$Cancer)
nod_new <-c()
for(i in 1:length(Cancers)){
	nod2 = nod[nod$Cancer== Cancers[i],]
	u<-unique(nod2$DuplicateModel)

	for(j in 1:length(u)){
		nod3 = nod2[nod2$DuplicateModel== u[j],]
		if(dim(nod3)[1]>MAX_SAMPLE_PER_PDX){
			nod4 = nod3[sample(nrow(nod3), MAX_SAMPLE_PER_PDX), ]		
			nod_new<-rbind(nod_new,nod4)
		}else{
			nod_new<-rbind(nod_new,nod3)
		}		
	}
}

nod<-nod_new

#--for each PDX, compute Tumor Purity diffrence by Late - Early passage within a PDX--#
TP_vector_all<-c()
P_vector_all<-c()
pdf(file='Fig1EF.S6.pdf')
Cancers = unique( nod$Cancer)
for(i in 1:length(Cancers)){
	nod2 = nod[nod$Cancer== Cancers[i],]
	u<-unique(nod2$DuplicateModel)

	TP_vector<-c()
	P_vector<-c()
	
	tm<- sapply(strsplit(nod2$TestModel, "-"), `[`, 2)
	tm<-gsub("R[0-9]+","",tm)		
	tm<-gsub("^P","",tm)
	tm<-as.numeric(tm)
	
	TP_vector<-nod2$TumorPurity
	P_vector<-tm
	
	TP_vector_all<-c(TP_vector_all,TP_vector)
	P_vector_all<-c(P_vector_all,P_vector)

	TP<-cbind(P_vector,TP_vector)
	colnames(TP)=  c('Passage','TumorPurity')
	TP<-as.data.frame(TP)
	
	TP = TP[TP[,1]<=MAX_PASSAGE,]
	TP<-na.omit(TP)
	
	if(dim(TP)[1]<MIN_SAMPLE2PLOT){next;}
	if(unique(TP[,1])<MIN_PASSAGE2PLOT){next;}
	gg<-ggscatter(TP, x="Passage", y="TumorPurity", 
			  xlab="Passage",
			  ylab="Tumor purity (%)",
			  main=Cancers[i],
			  xlim=c(0,max(TP[,1])+1),
			  shape = 21, size = 1,
			  add="reg.line",
			  conf.int=FALSE			  
			   ) +
			  stat_cor(size=6,label.x.npc = "center", label.y.npc = "bottom")+
			  theme_pubr(base_size = 20)
	
	print(gg)	
}
dev.off()


