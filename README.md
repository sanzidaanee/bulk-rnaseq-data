Introduction

 This project is for the  course of  Bioinformatics: Learn Bulk RNA-Seq Data Analysis From Scratch offered by udemy.

Description

The aim of this project is to perform a bioinformatics analysis on bulk RNA-seq samples for finding and characterizing differentially expressed genes that are important 
for a particular disease.The RNA-seq data are retrieved from the references 1, where the total RNA was extracted from control and dexamethasone treated ASM (Airway smooth 
muscle) cells. Dysregulation  of ASM function  creates various respiratory conditions such as asthma and chronic pulmonary disease. Analysis is done on eight cell lines. 


The raw FASTQ files were analyzed in Linux first for differential expression (DE) analysis. Assess the quality of the count data, and detect major sources of variation in the data 
and then performed differential analysis of genes using the DESeq2 package in R and visualized data like PCA, MS, HeatMap and Volcano Plots to find the differential expression of 
genes.

	 	 	 		
The DESeq2 R package is used in R to model the count data using a negative binomial model and test for differentially expressed genes. 
Visualization of the results with PCA,MA, heatmaps and volcano plots will be performed and the significantly differentially expressed genes 
are identified for further use in GO and pathway analysis to investigate the relationship with major biological processes and diseases. 


References

Himes, B. E., Jiang, X., Wagner, P., Hu, R., Wang, Q., Klanderman, B., ... & Lu, Q. (2014). RNA-Seq transcriptome profiling identifies CRISPLD2 as a glucocorticoid r
esponsive gene that modulates cytokine function in airway smooth muscle cells. PloS one, 9(6), e99625.
