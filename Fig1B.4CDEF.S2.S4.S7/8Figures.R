##############################################################################################################################
#Description: (1)use TCGA patient tumor to establish a non-linear non-parametric loess relationship between ESTIMATE score and 
#                Tumor Purity, so to extrapolate Tumor Purity for human cells of PDX tumors
#             (2)Use the same loess curve to infer Tumor Purity for 157 PDX with RNASeq and NGS-QC data, then plot correlation 
#				 between ESTIMATE-derived Tumor Purity and NGS-QC panel derived REAL Tumor purity to show that ESTIMATE algorithm
#				 is not very accurate for PDX RNAseq data.
#Author: Sheng Guo PhD  guosheng@crownbio.com
#Date:   Feb-2022
#Apache License
#Version 2.0, January 2004
#https://opensource.org/licenses/Apache-2.0
###############################################################################################################################
library(tidyverse)
library(data.table)
library(ggpubr)
library(TCGAbiolinks) #--Bioconductor version: Release (3.14) 
library(estimate)	  #--version 1.0.13. Last updated on 2016-09-26


#--get TCGA ESTIMATE SCORE AND PURITY
tcga<-fread('TCGA.txt')  	 #--source: https://bioinformatics.mdanderson.org/public-software/estimate/, version 1.0.3 downloaded on 2022-01
TP = Tumor.purity[,c(1,2,3)] #--from TCGAbiolinks
TP<-TP[TP[,3]!='NaN',]
TP[,3] = as.numeric(gsub(',','.',TP[,3]))
TP[,1] = as.character(TP[,1])
TP[,2] = as.character(TP[,2])
TP = as_tibble(TP)
names(TP)[1] = 'ID'

tcga$ID = paste(tcga$ID,'A',sep='')
tcga = left_join(TP,tcga)
tcga = na.omit(tcga)

########################################################################
#--nonparametric regression to fit the TUMOR_PURITY~ESTIMATE_SCORE curve
#--then use the fitted smooth line to get tumor purity for PDX samples 
########################################################################
#--to get extrapolation for loess fitting
tcga.lo <- loess(ESTIMATE ~ ESTIMATE_score, tcga, control = loess.control(surface = "direct"))
j <- order(tcga$ESTIMATE_score)
lines(tcga$ESTIMATE_score[j],tcga.lo$fitted[j],col='red',lwd=1.5)

#pdf(file='TCGA_ESTIMATEscore_vs_TumorPurity.pdf')	
pdf(file='FigS2.pdf')	
gg<-ggscatter(data=tcga, x="ESTIMATE_score", y="ESTIMATE", 
			  xlab="ESTIMATE score",
			  ylab="ESTIMATE-derived tumor purity (%)",
			  title='ESTIMATE score and tumor purity in 7098 TCGA samples'
			   ) +
   		      rotate_x_text(angle = 45) + 
			  theme(panel.grid.minor = element_line(colour="red", size=0.5), panel.grid.major = element_line(colour = "grey")) +
			  scale_y_continuous(minor_breaks = seq(0 , 1, 0.025), breaks = seq(0, 1, 0.025)) +
			  scale_x_continuous(minor_breaks = seq(-4000 , 5000, 500), breaks = seq(-4000, 5000, 500))
gg<- gg + geom_line(aes(x=tcga$ESTIMATE_score[j],y=tcga.lo$fitted[j]), size=1,color = "red", linetype = "solid") 
gg
dev.off()

#--save tcga data for reproducing the results in the future
write.csv(file='7098TCGAsamples_ESTIMATEscore_TumorPurity.csv', tcga)

#--read in PDX ESTIMATE scores calculated by R package 'estimate (version 1.0.13. Last updated on 2016-09-26)'
pdx1<-fread('PDX_Human.txt')    #--ESTIMATE score for PDX tumors using only human reads mapped to ESIMATE signature genes
pdx2<-fread('PDX_Mouse.txt')    #--ESTIMATE score for PDX tumors using only mouse reads that are mapped to human orthologs in ESIMATE signature genes
pdx3<-fread('PDX_Hybrid.txt') 	#--ESTIMATE score for PDX tumors using both human and mouse reads 

pdx<-left_join(pdx1,pdx2,by='NAME')
pdx<-left_join(pdx,pdx3,by='NAME')
pdx = pdx %>% arrange(by=NAME)

#--exclude leukemia, sarcoma etc due to heavey immune/stroma influence
pdx = pdx %>% filter(!grepl('AL', NAME))
pdx = pdx %>% filter(!grepl('AM', NAME))
pdx = pdx %>% filter(!grepl('XX', NAME))
pdx = pdx %>% filter(!grepl('LY', NAME))
pdx = pdx %>% filter(!grepl('GS', NAME))
pdx = pdx %>% filter(!grepl('SA', NAME))

pdx_ESTIMATE = pdx %>% select(NAME,ESTIMATEScore.x, ESTIMATEScore.y,ESTIMATEScore)
names(pdx_ESTIMATE) = c('PDX','Human','Mouse','Combined')

pdx_ESTIMATE = pdx_ESTIMATE %>%arrange(by=Human)
pdx_ESTIMATE = pdx_ESTIMATE %>% filter(Human <0)

pdx_ESTIMATE$Human.purity=predict(tcga.lo,pdx_ESTIMATE$Human)
pdx_ESTIMATE$Mouse.purity=predict(tcga.lo,pdx_ESTIMATE$Mouse)
pdx_ESTIMATE$Combined.purity=predict(tcga.lo,pdx_ESTIMATE$Combined)

#pdf(file='PDX_ESTIMATE_TumorPurity.pdf')
pdf(file='FigS4S7.pdf')
gghistogram(pdx_ESTIMATE$Human.purity*100, fill = "#00AFBB", xlab='Tumor purity (%)', ylab='Count',
			rug = TRUE) +
	theme_pubr(base_size = 20)

gghistogram(pdx_ESTIMATE$Mouse.purity*100, fill = '#00AFBB', xlab='Tumor purity (%)', ylab='Count',
			rug = TRUE) +
	theme_pubr(base_size = 20)

gghistogram(pdx_ESTIMATE$Combined.purity*100, fill = "#00AFBB", xlab='Tumor purity (%)', ylab='Count',
			rug = TRUE) +
	theme_pubr(base_size = 20)
dev.off()

########################################################################
#--check if there is any passage effect using first 10 passages due to 
#--sample size requirement
########################################################################
pdx.purity= cbind(pdx_ESTIMATE %>% select(PDX,Human,Mouse,Combined), pdx_ESTIMATE$Human.purity)
names(pdx.purity)[5]='TumorPurity'
pdx.purity = pdx.purity %>% mutate(Passage = str_extract(string =  PDX, pattern = "P[0-9]+"))
pdx.purity = pdx.purity %>% mutate(Passage = str_extract(string =  Passage, pattern = "[0-9]+"))
pdx.purity = pdx.purity %>% mutate(Passage = parse_number(Passage))

pdx.purity = pdx.purity %>% filter(Passage<11)
pdx.purity = pdx.purity %>% mutate(TumorPurity = TumorPurity*100)
pdx.purity = pdx.purity %>% mutate(Passage = paste('P', Passage,sep=''))

#pdf(file='ESTIMATE_TumorPurityByPassage_Human.pdf', width=8, height=4)
pdf(file='Fig1B.pdf', width=8, height=4)
g<-ggboxplot(pdx.purity, x = "Passage", y = "TumorPurity", color = "Passage", 
		  xlab='Passage', ylab="ESTIMATE-derived tumor purity (%)\n in human cells of PDX tumor",notch =FALSE,
          add = "jitter", legend = "none" , size=0.25,
		  order = paste('P',0:10,sep=''),
		  add.params = list(size = 0.2, alpha = 0.25)
		  )+
		rotate_x_text(angle = 45)
		#theme_pubr(base_size = 20)
ggpar(g, palette = "p1_aaas") # nature
dev.off()

#################################################################################################################
#Job 2: Use the same loess curve to infer Tumor Purity for 157 NOD/SCID PDX with RNASeq and NGS-QC data, 
#then plot correlation between ESTIMATE-derived Tumor Purity and NGS-QC panel derived REAL Tumor purity, 
#so to show that ESTIMATE algorithm is not very accurate for PDX RNAseq data.
#################################################################################################################
#--157 PDX with RNAseq data, from which ESTIMATE scores were calculated using R package 'estimate' version 1.0.13
#--meanwhile, same samples for RNAseq were also assasy by our NGS-QC panel to get accurate (error <1%) estimation
#--of tumor purity, for more information, see https://academic.oup.com/nargab/article/2/3/lqaa060/5893932
PDX157Data = fread(file='PDX_mouse_ratio.txt') 

PDX157Data = PDX157Data %>% select(TestSample,MouseRatio_Chip,ESTIMATEScore, ImmuneScore, StromalScore)
PDX157Data = PDX157Data %>% mutate(ESTIMATE_TumorPurity=predict(tcga.lo,PDX157Data$ESTIMATEScore)*100)
PDX157Data = PDX157Data %>% mutate(REAL_TumorPurity=100-MouseRatio_Chip)

#--spearman is more suitable due to outliers--#
cor.test(PDX157Data$REAL_TumorPurity, PDX157Data$ESTIMATE_TumorPurity, method='spearman')
cor.test(PDX157Data$REAL_TumorPurity, PDX157Data$ImmuneScore, method='spearman')
cor.test(PDX157Data$REAL_TumorPurity, PDX157Data$StromalScore, method='spearman')
cor.test(PDX157Data$REAL_TumorPurity, PDX157Data$ESTIMATEScore, method='spearman')

#pdf(file='ESTIMATE_vs_NGSassay.pdf')
pdf(file='Fig4C.pdf')
p<-ggscatter(PDX157Data, x="ESTIMATE_TumorPurity", y="REAL_TumorPurity", 		  
		  ylab="NGS assay-derived tumor purity (%)",
		  xlab="ESTIMATE-derived tumor purity (%)",		  
		  xlim=c(20,100),
		  ylim=c(20,100),
		  size=1.5,
		  add="reg.line",
		  add.params = list(size = 1.0, alpha = 0.25),
		  conf.int=FALSE)+
		  stat_cor(size=6, method='spearman', label.x.npc = "left", label.y.npc = "bottom")+
		  annotate("text", x = 60, y = 100, label = "",size=7, parse = TRUE)+
		  theme_pubr(base_size = 20)
p
dev.off()

#pdf(file='ImmuneScore_vs_NGSassay.pdf')
pdf(file='Fig4D.pdf')
p<-ggscatter(PDX157Data, x="ImmuneScore", y="REAL_TumorPurity", 		  
		  ylab="NGS assay-derived tumor purity (%)",
		  xlab="Immune score",		  
		  #xlim=c(20,100),
		  #ylim=c(20,100),
		  color='red',
		  size=1.5,
		  #add="reg.line",
		  add.params = list(size = 1.0, alpha = 0.25),
		  conf.int=FALSE)+
		  stat_cor(size=6, method='spearman', label.x.npc = "left", label.y.npc = "bottom")+
		  annotate("text", x = 60, y = 100, label = "",size=7, parse = TRUE)+
		  theme_pubr(base_size = 20)
p
dev.off()

#pdf(file='StromalScore_vs_NGSassay.pdf')
pdf(file='Fig4E.pdf')
p<-ggscatter(PDX157Data, x="StromalScore", y="REAL_TumorPurity", 		  
		  ylab="NGS assay-derived tumor purity (%)",
		  xlab="Stroma score",		  
		  #xlim=c(20,100),
		  #ylim=c(20,100),
		  color='green',
		  size=1.5,
		  #add="reg.line",
		  add.params = list(size = 1.0, alpha = 0.25),
		  conf.int=FALSE)+
		  stat_cor(size=6, method='spearman', label.x.npc = "left", label.y.npc = "bottom")+
		  annotate("text", x = 60, y = 100, label = "",size=7, parse = TRUE)+
		  theme_pubr(base_size = 20)
p
dev.off()


#pdf(file='ESTIMATEScore_vs_NGSassay.pdf')
pdf(file='Fig4F.pdf')
p<-ggscatter(PDX157Data, x="ESTIMATEScore", y="REAL_TumorPurity", 		  
		  ylab="NGS assay-derived tumor purity (%)",
		  xlab="ESTIMATE score",		  
		  #xlim=c(20,100),
		  #ylim=c(20,100),
		  color='blue',
		  size=1.5,
		  #add="reg.line",
		  add.params = list(size = 1.0, alpha = 0.25),
		  conf.int=FALSE)+
		  stat_cor(size=6, method='spearman', label.x.npc = "left", label.y.npc = "bottom")+
		  annotate("text", x = 60, y = 100, label = "",size=7, parse = TRUE)+
		  theme_pubr(base_size = 20)
p
dev.off()