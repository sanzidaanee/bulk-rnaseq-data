# Install Bioconductor and all other packages
`
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install ("GenomicFeatures")
BiocManager::install("DESeq2")   


# Import data file with creating variable for count and metadate

cnt <- read.csv("/Users/sanzidaakhteranee/Documents/Online_Course_24/Udemy_bulk rnaseq/Files 2/counts.csv")  #Feature count matrix that have obtained after processing the FASTQ files in Linux




met <- read.csv("/Users/sanzidaakhteranee/Documents/Online_Course_24/Udemy_bulk rnaseq/Files 2/metadata2.csv")  # metadata describe the sample that the RNA sequence was obtained from cell line





# Read the count and metadata file

str(cnt)

str(met)


# Matches the row names of metadata to column names in counts_data

all(colnames(cnt) %in% rownames(met))

#Checking order of row and column names

all(colnames(cnt)==rownames(met))


# Installation of DESeq2

BiocManager::install("DESeq2")





#Calling of DESeq2 Library

library(DESeq2)


#Building DESeq2 Dataset

dds <- DESeqDataSetFromMatrix(countData = cnt, 
                              colData = met, 
                              design = ~dexamethasone) # study was designed with control and  dexamethasone treatment

print(dds)



# Removal of low count reads 

keep <- rowSums(counts(dds)) >= 10




#Subset dataset

dds <- dds[keep, ] # all the genes of count reads less than 10 will remove

print(dds)  # thousands of rows have been dropped having count less then 10




# Setting reference for DEG  analysis

dds$dexamethasone <- relevel(dds$dexamethasone, ref= "untreated")
deg <- DESeq(dds) # set a variable for differential expression analysis
res <- results(deg) # result


#Save the results at the local folder in CSV file.

write.csv(res, "test_udemy.csv") # save result in new csv files

#Summary Statistics of results

summary(res)


# Result by changing p value consider alpha 0.05

new_result <- results(deg, alpha = 0.05)
summary(new_result)



# Install packages for converting gene id into gene name

BiocManager::install("org.Hs.eg.db") #organism-specific annotation package for Homo sapiens (human) and map Ensembl IDs to gene symbols:

library(org.Hs.eg.db)


# Create a dataframe

new_result.df <- as.data.frame(new_result)
str(new_result.df)



#Transfer gene id to gene name (create a new column named Symbol in this data frame)

new_result.df$Symbol <- mapIds (org.Hs.eg.db, rownames (new_result.df), 
                                keytype = "ENSEMBL", column = "SYMBOL")
new_result.df
str(new_result.df)


#save result into a new file

write.csv(new_result.df, "final_test_udemy.csv")

`



#Install packages for building of PCA Plot

install.packages("gtable", dependencies = TRUE)



if (!requireNamespace("BiocManager", force = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

library(DESeq2)




# PCA plot

vsd <- vst(deg, blind= FALSE) # vst function is used to apply a variance stabilizing transformation to count data and transformation will take the experimental design into account

plotPCA(vsd, intgroup= "dexamethasone")



#Size factor estimation

sizeFactors(deg) # normalization factors that account for differences in sequencing depth or library size across samples



#Estimating the dispersion

plotDispEsts(deg)



#Install Tidyverse and ggplot packages

install.packages("dplyr")
install.packages("ggplot2")



#Calling library

library(dplyr)
library(ggplot2)


#Building MA plot

plotMA(new_result)







