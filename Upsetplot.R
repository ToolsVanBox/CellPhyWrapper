#!/usr/bin/env Rscript

### Empty the variables
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)

### Load libraries
library(VariantAnnotation)
library(UpSetR)

### get command arguments
args = commandArgs(trailingOnly = TRUE)
vcf_file = args[1]
outputdir = args[2]
prefix = args[3]

### Read file 
myvcf <- readVcf(vcf_file)
#rownames(myvcf) <- paste(as.character(seqnames(myvcf)),start(myvcf),sep=":")

### Fix the haplotypes in binairy fashion 
myvcf_gt_df <- as.data.frame(geno(myvcf)$GT)  
myvcf_gt_df[is.na(myvcf_gt_df)] <- 0
myvcf_gt_df[myvcf_gt_df == '0/0'] <- 0
myvcf_gt_df[myvcf_gt_df == '0|0'] <- 0
myvcf_gt_df[myvcf_gt_df == './.'] <- 0
myvcf_gt_df[myvcf_gt_df == '.|.'] <- 0
myvcf_gt_df[myvcf_gt_df == '0/1'] <- 1
myvcf_gt_df[myvcf_gt_df == '0|1'] <- 1
myvcf_gt_df[myvcf_gt_df == '1/1'] <- 1
myvcf_gt_df[myvcf_gt_df == '1|1'] <- 1

### Get upset matrix 
upset_matrix_raw <- matrix(as.numeric(unlist(myvcf_gt_df)),ncol=ncol(myvcf_gt_df))
rownames(upset_matrix_raw) <- rownames(myvcf_gt_df)
colnames(upset_matrix_raw) <- colnames(myvcf_gt_df)
upset_matrix <- upset_matrix_raw[which(rowSums(upset_matrix_raw) > 1),]

### Plot figures 
pdf(paste0(outputdir, "/", prefix, "_all_upset.pdf"), onefile=FALSE)
upset(as.data.frame(upset_matrix_raw),nsets=ncol(upset_matrix_raw), order.by = "freq")
dev.off()
pdf(paste0(outputdir, "/", prefix, "_shared_upset.pdf"), onefile=FALSE)
upset(as.data.frame(upset_matrix),nsets=ncol(upset_matrix), order.by = "freq")
dev.off()








