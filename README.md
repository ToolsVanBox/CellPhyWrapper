# CellPhyWrapper
Wrapper for running and plotting trees while using the algorithm of CellPhy. CellPhy is really good in determining the shape of a phylogetic tree, however it assigns mutations to tell where certain cells split of. In a cancer setting it is important that all mutations are assigned correctly. That is why this wrapper has been written. It consists of 2 parts: 1) running CellPhy 2) Reassigning mutations and plotting trees. 


# Processing input: 
Input for Cellphy wrapper is a multi sample vcf which only contains Single Nucleotide Variants (SNVs). All data should be of the same data type, with the only exception being your normal-control sample. Meaning that if you have PTA samples and Bulk sequencing, you should remove the bulk sample from you multisample VCF file. An example of how you could do this is in: CreateInputVCFexample.sh


# help function CellPhyWrapper 
If you want to know which function are available run: 
'''
cellPhyPipeline.sh --help
'''


# Running CellPhyWrapper 
An example of how to run CellPhyWrapper. The minimal required arguments are:
'''
sbatch cellPhyPipeline.sh --input input.vcf --outgr CTRL-SAMPLE --outputdir /Your/Location/
'''
If you want to redo the plotting you could change "--plottingonly" to true and this will skip the computations heavy part of creating the trees
'''
sbatch cellPhyPipeline.sh --input input.vcf --outgr CTRL-SAMPLE --outputdir /Your/Location/ --prefix SampelName --percentage 0.4 --plottingonly false
'''

# Add-on scripts 
For the additional function for plotty and making your figures prettier install the add-on scripts: 
'''
devtools::install_local("./AddOnScripts", force = T)
'''

# Ultrametric trees
In order to create ultrametric from the tree objects install the following functions: 
'''
devtools::install_local("./Ultrametric", force = T)
'''


# Publications:
Original paper where CellPhy is described: 
Alexey Kozlov, João M Alves, Alexandros Stamatakis, and David Posada (2022) **CellPhy: accurate and fast probabilistic inference of single-cell phylogenies from scDNA-seq data** *Genome Biol 23, 37* doi: [10.1186/s13059-021-02583-w](https://doi.org/10.1186/s13059-021-02583-w)

