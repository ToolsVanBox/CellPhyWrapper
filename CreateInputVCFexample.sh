#!/bin/bash

# Load modules
module load bedtools 
module load bcftools

# Select all mutations after PTATO filtering 
bedtools intersect -a /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt2283D/ptato_vcfs/pt2283D/pt2283D.ptato.merged.vcf.gz -b /hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/pt2283D/snvs/pt2283D/*/*vcf.gz -wa -header -u > /hpc/pmc_vanboxtel/projects/TallFlow/1_Input/TreeBuilding/pt2283D/pt2283D_filtered.vcf
# Select only SNVs
bcftools view -m2 -M2 -v snps /hpc/pmc_vanboxtel/projects/TallFlow/1_Input/TreeBuilding/pt2283D/pt2283D_filtered.vcf > /hpc/pmc_vanboxtel/projects/TallFlow/1_Input/TreeBuilding/pt2283D/pt2283D_filtered_snv.vcf
# Grep the header of the vcf file
grep -P "^##" /hpc/pmc_vanboxtel/projects/TallFlow/1_Input/TreeBuilding/pt2283D/pt2283D_filtered_snv.vcf > /hpc/pmc_vanboxtel/projects/TallFlow/1_Input/TreeBuilding/pt2283D/pt2283D_filtered_snv_NoBulk.vcf
# Grep all selected columns, cutting out the bulk sample (location 10 in this example)
grep -P "^#CHROM" /hpc/pmc_vanboxtel/projects/TallFlow/1_Input/TreeBuilding/pt2283D/pt2283D_filtered_snv.vcf | cut -f -9,11- >> /hpc/pmc_vanboxtel/projects/TallFlow/1_Input/TreeBuilding/pt2283D/pt2283D_filtered_snv_NoBulk.vcf
# Grep all mutations of the file, again cutting out the bulk sample (location 10 in this example)
grep -v "^#" /hpc/pmc_vanboxtel/projects/TallFlow/1_Input/TreeBuilding/pt2283D/pt2283D_filtered_snv.vcf | cut -f -9,11- >> /hpc/pmc_vanboxtel/projects/TallFlow/1_Input/TreeBuilding/pt2283D/pt2283D_filtered_snv_NoBulk.vcf



