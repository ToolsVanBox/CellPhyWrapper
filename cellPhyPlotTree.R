#!/usr/bin/env Rscript

.getSourceDir <- function() {
  cmdArgs <- commandArgs(trailingOnly=FALSE)
  fileArg <- "--file="
  match <- grep(fileArg, cmdArgs)
  if (length(match) > 0) {
    sourcedir <- dirname(gsub(fileArg, "", cmdArgs[match]))
  }
  return( sourcedir )
}

sourcedir <- .getSourceDir()
print(sourcedir)

# Load libraries
source(paste(sourcedir, "/cellPhyUtils.R",sep=""), chdir = T)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(VariantAnnotation)
ref_genome <- "BSgenome.Hsapiens.NCBI.GRCh38"
library(ref_genome, character.only = TRUE)

# get command arguments
args = commandArgs(trailingOnly = TRUE)
cellphydir = args[1]
vcf_file = args[2]
ptato_dir = args[3]
outgroup = args[4]

# read in the PTATO VCFs
if (ptato_dir != "NONE") {
  ptato_grl = GRangesList(xxxx) # HERE THE FILTERED PTATO VCFS SHOULD BE EXTRACTED FROM THE DIR AND PUT INTO ONE GRANGESLIST
} else {
  ptato_grl = NA
}

# read in the SMuRF VCF
vcf = readVcf(vcf_file)
seqlevels(vcf, pruning.mode = 'tidy') = c(1:22, "X", "Y")

# If no VAF column is present create one 
if (is.null(geno(vcf)$VAF)){
  print("Empty VAF")
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

# load the tree, all annotations (muts, boots), and filter end-branches
cell_phy_tree <- load_tree_with_info(
  dir = cellphydir,
  outgr = outgroup,
  vcf = vcf,
  ptato_grl = ptato_grl, 
  cellphy_rm_non1 = TRUE,
  mutation_soure = 'cellphy', 
  norm_pres_max = 0.95,  high_frac_min = 0.05, min_frac_all = 0.5
)

# plot the tree
tree_plot = ggtree(cell_phy_tree, branch.length = 'branch_length') + #, aes(color = group)) +
  geom_nodelab(aes(label = n_boot), geom = 'label', color = 'grey50', fill = rgb(1,1,1,0.7), size = 2) +
  geom_text(aes(x = branch, label = branch_length), color = 'black', nudge_y = -0.5, nudge_x = 1) +
  geom_tiplab(size=3) + 
  coord_cartesian(clip = 'off') +
  theme(plot.margin = unit(c(10,100,10,10), 'points'))

ggsave(plot = tree_plot, filename = paste0(cellphydir, "/CellPhyWrapperTree.pdf"), 
       width = 10, height = 7)

saveRDS(cell_phy_tree, file = paste0(cellphydir, "/TreeObject.RDS"))

### HERE POSSIBLY CODE WITH REFITTING STANDARD COSMIC SIGNATURES TO BRANCHES


#### ToDo for monday: 
# Add your changed to the actual script
# Save RDS object of tree
# Add PTATO filtering
# Can we create branch specific VCF files from CellPhy? 


