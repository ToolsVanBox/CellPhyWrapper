########## THIS SHOULD NOT BE HARD CODED
source("/hpc/pmc_vanboxtel/tools/external/cellphy_scripts/cellPhyUtils")
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(VariantAnnotation)

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

# load the tree, all annotations (muts, boots), and filter end-branches
load_tree_with_info(
  dir = cellphydir,
  outgr = outgroup,
  vcf = vcf,
  ptato_grl = ptato_grl, 
  cellphy_rm_non1 = TRUE,
  mutation_soure = 'cellphy', 
  norm_pres_max = 0.95,  high_frac_min = 0.05, min_frac_all = 0.5
)

# plot the tree
tree_plot = ggtree(tree_list_cp2[[n_tree]], branch.length = 'branch_length', aes(color = group)) +
  geom_text(aes(x = branch, label = n_boot), color = 'grey50', nudge_y = 0.5, nudge_x = 1) +
  geom_text(aes(x = branch, label = branch_length), color = 'black', nudge_y = -0.5, nudge_x = 1) +
  geom_tiplab(size=3) + 
  coord_cartesian(clip = 'off') +
  theme(plot.margin = unit(c(10,100,10,10), 'points'))

ggsave(plot = tree_plot, filename = paste0(cellphydir, "/CellPhyTreeInitial.pdf"), 
       width = 10, height = 7)


### HERE POSSIBLY CODE WITH REFITTING STANDARD COSMIC SIGNATURES TO BRANCHES




