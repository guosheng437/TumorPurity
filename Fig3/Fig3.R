########################################################################################################################
#Description: generate Fig3F-G for tumor purity by cancer
#Author: Sheng Guo PhD  guosheng@crownbio.com
#Date:   Aug-2021
#Apache License
#Version 2.0, January 2004
#https://opensource.org/licenses/Apache-2.0
#
#Input files:
#	1. Locallo2019.txt  Tumor purity estimated by the ESTIMATE algorithm, taken from
#	   	Locallo, A., Prandi, D., Fedrizzi, T. and Demichelis, F. (2019) TPES: tumor purity estimation from SNVs. 
#		Bioinformatics, 35, 4433-4435.
#
#		Note: Tumor purity estimated by the ESTIMATE algorithm is from the paper below via TCGAbiolinks
#		Aran D, Sirota M, Butte AJ. Systematic pan-cancer analysis of tumour purity. 
#		Nat Commun. Nature Publishing Group; 2015;6: 8971. pmid:26634437
#	2. NOD_medianTumorPurityByCancer.csv: median tumor purity by cancers for NOD/SCID mice estimated by Crown's 
#	   deep NGS assay with estimation error <1%, method is described in
#		Chen, X., Qian, W., Song, Z., Li, Q.X. and Guo, S. (2020) Authentication, characterization and contamination 
#		detection of cell lines, xenografts and organoids by barcode deep NGS sequencing. NAR Genom Bioinform, 2, lqaa060
########################################################################################################################
library(TCGAbiolinks) #Bioconductor version: Release (3.14)
library(ggpubr)
library(tidyverse)

MIN_TCGA_SAMPLE=20

#--read in NOD/SCID tumor purity
nod<-read.csv(file='NOD_medianTumorPurityByCancer.csv');

#--create cancer mapping between PDX and TCGA
nod_cancer = c("KI","HN","LU","BL","CV","OV","LI","UT","CR","BR","BN","ME")
tcga_cancer = c("KIRC","HNSC","LUAD+LUSC","BLCA","CESC","OV","LIHC","UCEC+UCS","COAD","BRCA","LGG+GBM","SKCM")

#--there are 5 methods to calculate tumor purity in TCGAbiolinks
#ESTIMATE	ABSOLUTE	LUMP	IHC	CPE
#METHOD=c('ABSOLUTE','ESTIMATE','LUMP','IHC','CPE');
#ABSOLUTE: too few data points to use, so only ESIMATE based tumor purity is used
#METHOD=c('ABSOLUTE','ESTIMATE','LUMP','IHC','CPE');
METHOD=c('ESTIMATE');

for(m in 1:length(METHOD)){
	tcga_tp<-c()
	for(i in 1:length(nod_cancer)){

		a<-c()
		if(i==3){
			a1<-Tumor.purity[ Tumor.purity$Cancer.type=='LUAD', METHOD[m]]
			a2<-Tumor.purity[ Tumor.purity$Cancer.type=='LUSC', METHOD[m]]
			a1<-as.vector(a1)
			a1<-a1[a1!='NaN']
			a1<-gsub(',','.',a1)
			a1<-as.numeric(a1)
			a2<-as.vector(a2)
			a2<-a2[a2!='NaN']
			a2<-gsub(',','.',a2)
			a2<-as.numeric(a2)
			a<-c(a1,a2)
		}else if(i==8){
			a1<-Tumor.purity[ Tumor.purity$Cancer.type=='UCS', METHOD[m]]
			a2<-Tumor.purity[ Tumor.purity$Cancer.type=='UCEC', METHOD[m]]
			a1<-as.vector(a1)
			a1<-a1[a1!='NaN']
			a1<-gsub(',','.',a1)
			a1<-as.numeric(a1)
			a2<-as.vector(a2)
			a2<-a2[a2!='NaN']
			a2<-gsub(',','.',a2)
			a2<-as.numeric(a2)
			a<-c(a1,a2)
			a<-c(a1,a2)
		}else if(i==11){
			a1<-Tumor.purity[ Tumor.purity$Cancer.type=='LGG', METHOD[m]]
			a2<-Tumor.purity[ Tumor.purity$Cancer.type=='GBM', METHOD[m]]
			a1<-as.vector(a1)
			a1<-a1[a1!='NaN']
			a1<-gsub(',','.',a1)
			a1<-as.numeric(a1)
			a2<-as.vector(a2)
			a2<-a2[a2!='NaN']
			a2<-gsub(',','.',a2)
			a2<-as.numeric(a2)
			a<-c(a1,a2)
			a<-c(a1,a2)
		}else{
			a<-Tumor.purity[ Tumor.purity$Cancer.type==tcga_cancer[i], METHOD[m]]
			a<-as.vector(a)
			a<-a[a!='NaN']
			a<-gsub(',','.',a)
			a<-as.numeric(a)
		}
		
		if(length(a)<MIN_TCGA_SAMPLE){
			print(paste("warning: too few samples for ", nod_cancer[i]))	
		}
		
		tcga_tp<-c(tcga_tp,median(a))
	}

	nod_tp <- nod[nod$Cancer %in% nod_cancer,3]
	tcga_tp = tcga_tp *100

	median(tcga_tp)
	median(nod_tp)

	TP<-cbind(tcga_tp,nod_tp)
	rownames(TP)<-nod_cancer
	TP<-as.data.frame(TP)

	#pdf(file=paste('TCGA_NOD_TP_Aran2015_', METHOD[m], '.pdf', sep=''))
	pdf(file='Fig3G.pdf')
	gg<-ggscatter(TP, x="tcga_tp", y="nod_tp", 
			  xlab=paste(METHOD[m], "-derived TCGA tumor purity (%)", sep=""),
			  #ylab="NOD/SCID PDX tumor purity (%)",
			  ylab="PDX tumor purity (%)",
			  #xlim=c(65,92),
			  #ylim=c(65,92),
			  add="reg.line",
			  conf.int=FALSE,
			  label=rownames(TP),
			  repel=T
			   ) +
			  stat_cor(size=6)+
			  theme_pubr(base_size = 20)
	#gg<-gg +  grids(linetype = "dashed")
	#gg<-gg +geom_abline(slope=1,intercept=0,col='red',lwd=1,lty=2)
	print(gg)
	dev.off()
}

#--read in Locallo2019 ABSOLUTE data--#
tcga_locallo<-read.table("Locallo2019.txt",head=TRUE)
table(tcga_locallo$TCGA.project)
tcga_locallo_tp<-aggregate(data=tcga_locallo, ABSOLUTE_Locallo2019~TCGA.project, FUN=mean)

nod_cancer = c( "KI","PA","HN","LU","GA","BL","CV","OV","LI","UT","CR","BR","BN","ME")
tcga_cancer = c("KIRC","PAAD","HNSC","LUAD+LUSC","STAD", "BLCA","CESC","OV","LIHC","UCEC+UCS","COAD","BRCA","LGG+GBM","SKCM")

METHOD='ABSOLUTE'

tcga_tp<-c()
for(i in 1:length(nod_cancer)){

	a<-c()
	if(i==4){
		a1<-tcga_locallo_tp[tcga_locallo_tp$TCGA.project=='LUAD',2]
		a2<-tcga_locallo_tp[tcga_locallo_tp$TCGA.project=='LUSC',2]
		a1<-as.vector(a1)
		a1<-a1[a1!='NaN']
		a1<-gsub(',','.',a1)
		a1<-as.numeric(a1)
		a2<-as.vector(a2)
		a2<-a2[a2!='NaN']
		a2<-gsub(',','.',a2)
		a2<-as.numeric(a2)
		a<-c(a1,a2)
	}else if(i==10){
		a1<-tcga_locallo_tp[tcga_locallo_tp$TCGA.project=='UCS',2]
		a2<-tcga_locallo_tp[tcga_locallo_tp$TCGA.project=='UCEC',2]
		a1<-as.vector(a1)
		a1<-a1[a1!='NaN']
		a1<-gsub(',','.',a1)
		a1<-as.numeric(a1)
		a2<-as.vector(a2)
		a2<-a2[a2!='NaN']
		a2<-gsub(',','.',a2)
		a2<-as.numeric(a2)
		a<-c(a1,a2)
		a<-c(a1,a2)
	}else if(i==13){
		a1<-tcga_locallo_tp[tcga_locallo_tp$TCGA.project=='LGG',2]
		a2<-tcga_locallo_tp[tcga_locallo_tp$TCGA.project=='LGG',2]
		a1<-as.vector(a1)
		a1<-a1[a1!='NaN']
		a1<-gsub(',','.',a1)
		a1<-as.numeric(a1)
		a2<-as.vector(a2)
		a2<-a2[a2!='NaN']
		a2<-gsub(',','.',a2)
		a2<-as.numeric(a2)
		a<-c(a1,a2)
		a<-c(a1,a2)
	}else{
		a<-tcga_locallo_tp[tcga_locallo_tp$TCGA.project==tcga_cancer[i],2]
		a<-as.vector(a)
		a<-a[a!='NaN']
		a<-gsub(',','.',a)
		a<-as.numeric(a)
	}
	
	
	tcga_tp<-c(tcga_tp,median(a))
}

nod_tp <- nod[nod$Cancer %in% nod_cancer,3]
tcga_tp = tcga_tp *100

median(tcga_tp)
median(nod_tp)

TP<-cbind(tcga_tp,nod_tp)
rownames(TP)<-nod_cancer
TP<-as.data.frame(TP)

#pdf(file=paste('TCGA_NOD_TP_Locallo2019_', METHOD, '.pdf', sep=''))
pdf(file='Fig3F.pdf')
gg <- ggscatter(TP, x="tcga_tp", y="nod_tp", 
		  xlab="ABSOLUTE-derived TCGA tumor purity (%)",
		  ylab="PDX tumor purity (%)",
		  add="reg.line",
		  conf.int=FALSE,
		  label=rownames(TP),
		  repel=T
		   ) +
		  stat_cor(size=6)+
		  theme_pubr(base_size = 20)
gg
dev.off()

 
########################################################################################################################
#--get an overall comparison between different methods
########################################################################################################################
b1<-as.tibble(read.table("Locallo2019.txt",head=TRUE)) #-read Tumor purity from Locallo2019 estimated by ABSOLUTE 
b1 = b1 %>% select(sample,ABSOLUTE_Locallo2019)
xx=table(b1$sample) == 1 #remove samples with more than one occurrence
xx=names(xx[xx==TRUE])
b1 = b1 %>% filter(sample %in% xx)

a1<-as.tibble(Tumor.purity)
a1<-a1 %>% filter(ABSOLUTE != 'NaN') %>% filter(ESTIMATE != 'NaN') %>% filter(LUMP != 'NaN') %>% filter(IHC != 'NaN') %>% filter(CPE != 'NaN')
a1$ESTIMATE<-as.numeric(gsub(',','.',a1$ESTIMATE))
a1$ABSOLUTE<-as.numeric(gsub(',','.',a1$ABSOLUTE))
a1$LUMP<-as.numeric(gsub(',','.',a1$LUMP))
a1$IHC<-as.numeric(gsub(',','.',a1$IHC))
a1$CPE<-as.numeric(gsub(',','.',a1$CPE))
a1 = a1 %>% select(Sample.ID, ESTIMATE,ABSOLUTE,LUMP,IHC) 
names(a1)[1] = 'sample'
a1 = a1 %>% mutate(sample = gsub('.{1}$', '',sample))


commonsamples = intersect(a1$sample,b1$sample)
a1=a1 %>% filter(sample %in% commonsamples)
b1=b1 %>% filter(sample %in% commonsamples)

a1=a1 %>%left_join(b1)

a1 = a1 %>% pivot_longer(!sample, names_to = "Method", values_to = "TumorPurity")
a1 = a1 %>% mutate(TumorPurity = TumorPurity*100)


pdf(file='TCGA_TumorPurity_5Methods.pdf', width=8, height=4)
g<-ggboxplot(a1, x = "Method", y = "TumorPurity", color = "Method", 
		  xlab='Estimation method', ylab='Tumor purity (%)',notch =FALSE,
          add = "jitter", legend = "none" , size=0.25,	
		  order = c('ABSOLUTE_Locallo2019','ABSOLUTE','LUMP','ESTIMATE','IHC'),		  
		  add.params = list(size = 0.2, alpha = 0.25)
		  )+
		font('x.text',size=6)+
		rotate_x_text(angle = 45)
		#theme_pubr(base_size = 20)
ggpar(g, palette = "p1_aaas")
dev.off()

#--compare ABSOLUTE tumor purity between Locallo2019 and Aran2015
b1<-as.tibble(read.table("Locallo2019.txt",head=TRUE)) #-read Tumor purity from Locallo2019 estimated by ABSOLUTE 
b1= b1 %>% select(sample,ABSOLUTE_Locallo2019)
xx=table(b1$sample) == 1 #remove samples with more than one occurrence
xx=names(xx[xx==TRUE])
b1 = b1 %>% filter(sample %in% xx)


a1<-as.tibble(Tumor.purity)
a1 = a1 %>% select(Sample.ID, ABSOLUTE)
a1<-a1 %>% filter(ABSOLUTE != 'NaN') 
a1$ABSOLUTE<-as.numeric(gsub(',','.',a1$ABSOLUTE))
names(a1)[1] = 'sample'
a1 = a1[!grepl("B$", a1$sample, perl=TRUE),] 

a1 = a1 %>% mutate(sample = gsub('.{1}$', '',sample))

commonsamples = intersect(a1$sample,b1$sample)
a1=a1 %>% filter(sample %in% commonsamples)
b1=b1 %>% filter(sample %in% commonsamples)

a1=a1 %>%left_join(b1)

a1 = a1 %>% mutate(ABSOLUTE = ABSOLUTE*100,ABSOLUTE_Locallo2019=ABSOLUTE_Locallo2019*100)

#pdf(file="ABSOLUTE_Aran2015_vs_Locallo2019_.pdf")
gg <- ggscatter(a1, x="ABSOLUTE", y="ABSOLUTE_Locallo2019", 
		  xlab="Tumor purity (%) by ABSOLUTE in Aran2015",
		  ylab="Tumor purity (%) by ABSOLUTE in Locallo2019",
		  add="reg.line",
		  size=1,
		  conf.int=FALSE
		   ) +
		  stat_cor(size=6,method='spearman')+
		  theme_pubr(base_size = 20)
#gg<-gg + grids(linetype = "dashed")
gg<-gg + geom_abline(slope=1,intercept=0,col='red',lwd=1,lty=2)
gg
#dev.off()
