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
outgroup = args[4]

### Read file and add VAF if needed
myvcf <- readVcf(vcf_file)
# If no VAF column is present create one 
if (is.null(geno(vcf)$VAF)){
  mydf <- data.frame(matrix(NA, nrow = nrow(vcf)))
  # Calculate VAF per sample 
  for (sname in samples(header(vcf))) {
    ad <- sapply(geno(vcf)$AD[,sname],"[[",2)
    dp <- geno(vcf)$DP[,sname]
    VAF <- ad/dp
    mydf[[sname]] <- as.numeric(VAF)
  }
  # Add VAF to vcf
  mydf <- mydf[-1]
  geno(header(vcf))["VAF",] = list("1","Float","Varient Allele Frequency")
  geno(vcf)$VAF <- mydf
}

### get mutation VAF
vcf_names = gsub("_.*", "", names(vcf))
vaf_all = vcf@assays@data$VAF %>% apply(2, unlist)
vaf_all[is.na(vaf_all)] = 0
colnames(vaf_all) = samples(header(vcf))
### Create an upset matrix based upon vaf. 
# Vaf > 0 = present
upset_matrix_raw = lapply(colnames(vaf_all), function(cn) {
  ifelse(vaf_all[ ,cn] > 0, 1, 0)
}) %>% do.call(cbind, .) %>% `colnames<-`(colnames(vaf_all))
rownames(upset_matrix_raw) <- rownames(vcf)
upset_matrix <- upset_matrix_raw[which(rowSums(upset_matrix_raw) > 1),]

### Plot figures 
pdf(paste0(outputdir, "/", prefix, "_all_upset.pdf"), onefile=FALSE)
upset(as.data.frame(upset_matrix_raw),nsets=ncol(upset_matrix_raw), order.by = "freq", keep.order = T, 
      sets = colnames(upset_matrix_raw[,!colnames(upset_matrix_raw) %in% c(outgroup)]))
dev.off()
pdf(paste0(outputdir, "/", prefix, "_shared_upset.pdf"), onefile=FALSE)
upset(as.data.frame(upset_matrix_raw),nsets=ncol(upset_matrix_raw), order.by = "freq", keep.order = T, 
      sets = colnames(upset_matrix_raw[,!colnames(upset_matrix_raw) %in% c(outgroup)]))
dev.off()








