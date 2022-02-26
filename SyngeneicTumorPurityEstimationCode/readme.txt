###################################################################################
#Description: README for computer code accompanying manuscript 
# 		"Tumor Purity in Preclinical Mouse Tumor Models"
#Author: Sheng Guo PhD  guosheng@crownbio.com
#Date:   Aug-2021
#
#
#Apache License
#Version 2.0, January 2004
#https://opensource.org/licenses/Apache-2.0
###################################################################################
This fold contains perl code and demo data for estimating tumor purity in syngeneic tumors


1. Running environment:
	Windows/Linus OS with perl installation (tested version: 5.28.1)


2. Description
This script implements the maximal likelihood algorithm described in the section 
"Tumor purity estimation for syngeneic models" of MATERIALS AND METHODS

It requireds three input files:
	--mouse_strain_mutfre: record the SNP identity and frequency (always one for a nucleotide) in a mouse strain
	--syngeneic_mutfreq:   record the SNP identity and frequency of somatic mutations in a syngeneic model

	the above two files contains same SNPs
	
	--tumor_depth_file:    this is the read depth file for a syngeneic tumor profiled by a deep NGS assay reported in
	Chen, X., Qian, W., Song, Z., Li, Q.X. & Guo, S. Authentication, characterization and contamination detection of 
	cell lines, xenografts and organoids by barcode deep NGS sequencing. NAR Genom Bioinform 2, lqaa060 (2020).	
	For each SNP (e.g. 1_39985259), it records the number of reads on its four nucleotides in the following format:
	
	pos	refbase	Adepth	Tdepth	Cdepth	Gdepth
	1_39985259	T	32	223	0	0		

3. Example
	perl TumorPurityEstimation.pl EMT6.mutfre EMT6.somatic.mutfre EMT6-1.depth
	
	output: 
	EMT6-1.depth	0.626	224
	
	output explanation:
	EMT6-1.depth: syngeneic tumor depth file for its estimation of tumor purity
	0.626: tumor purity
	224:   number of informative SNPs used in the maximal-likelihood calculation
	