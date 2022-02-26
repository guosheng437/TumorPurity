#################################################################################################
#Description: generate Fig4A-B
#Author: Sheng Guo PhD  guosheng@crownbio.com
#Date: Aug-2021
#Apache License
#Version 2.0, January 2004
#https://opensource.org/licenses/Apache-2.0
#################################################################################################
library(ggpubr)

a<-read.table('PDX_mouse_ratio.txt', sep='\t',header=TRUE, row.names=1)

a$TP_WES = 100-a$MouseRatio_WES
a$TP_RNAseq = 100-a$MouseRatio_RNAseq
a$TP_Chip = 100-a$MouseRatio_Chip

cor(a$TP_Chip, a$TP_WES)
cor(a$TP_Chip, a$TP_RNAseq)

pdf(file='Fig4A.pdf')
ggscatter(a, x="TP_WES", y="TP_Chip", 
		  #ylab="Ground truth tumor purity (%)",
		  ylab="NGS assay-derived tumor purity (%)",
		  xlab="WES-derived tumor purity (%)",		  
		  #xlim=c(0,16.5),
		  #ylim=c(0,16.5),
		  size=1.5,
		  add="loess",
		  conf.int=FALSE) +
		  #stat_cor(size=6)+
		  annotate("text", x = 70, y = 100, label = "italic(R) == 0.93",size=7, parse = TRUE)+
		  theme_pubr(base_size = 20)
dev.off()

pdf(file='Fig4B.pdf')
ggscatter(a, x="TP_RNAseq", y="TP_Chip", 
		  #ylab="Ground truth tumor purity (%)",
		  ylab="NGS assay-derived tumor purity (%)",
		  xlab="RNAseq-derived tumor purity (%)",		  
		  xlim=c(70,100),
		  #ylim=c(0,16.5),
		  size=1.5,
		  add="reg.line",
		  #add.params = list(size = 1.0, alpha = 0.25),
		  conf.int=FALSE) +
		  #stat_cor(size=6)+
		  annotate("text", x = 75, y = 100, label = "italic(R) == 0.61",size=7, parse = TRUE)+
		  theme_pubr(base_size = 20)
dev.off()
