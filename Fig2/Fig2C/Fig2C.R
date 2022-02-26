#################################################################################################
#Description: generate Fig2C
#Author: Sheng Guo PhD  guosheng@crownbio.com
#Date: Aug-2021
#Note: only use NOD and BALB/c mice with passage information
#Apache License
#Version 2.0, January 2004
#https://opensource.org/licenses/Apache-2.0
#################################################################################################
rm(list=ls())
library(ggpubr)

TP<-read.table('TumorPurity.txt', header=TRUE, row.names=1)


pdf(file='Fig2C.pdf')
ggscatter(TP, x="TumorPurity_Balbc", y="TumorPurity_NOD", 
		  xlab="BALB/c nude tumor purity (%)",
		  ylab="NOD/SCID tumor purity (%)",
		  add="reg.line",
		  conf.int=FALSE,
		  label=rownames(TP)
		   ) +
		  stat_cor(size=6)+
		  theme_pubr(base_size = 20)
dev.off()
