# required libraries
require(ggtree)
require(ape)
require(UpSetR)
require(tidyverse)

# read the initial --------------------------------------------
read_ggtree <- function(file) {
  tree_text = readLines(file)
  tree_text = gsub('(:[0-9\\.eE\\+\\-]+)\\[(\\d+)\\]', '\\1\\[&&NHX:N=\\2\\]', tree_text)
  tree_text = gsub(';', '\\[&&NHX:N=-1\\];', tree_text)
  # print(tree_text)
  tree = treeio::read.nhx(textConnection(tree_text))
  return(tree)
}

# root the tree --------------------------------------------
root_tree <- function(tree, outgr, add_info = FALSE) {
  if (outgr != "NONE") {
    tree = treeio::root(tree, outgroup = outgr, resolve.root = TRUE, edgelabel = TRUE)
    tree@phylo$tip.label = tree@data$tip.label[order(tree@data$node)][1:length(tree@phylo$tip.label)]
  } else {
    tree = tree
  }
  if (add_info) { tree = add_sample_info(tree) }
  return(tree)
}

# merge the cosmic/dbsnp doulbe-annotated mutation names -------------------------------
merge_cosmic_muts <- function(mutations, vcf_names = NULL, vcf = NULL) {
  stopifnot(!is.null(vcf_names) || !is.null(vcf))
  if (is.null(vcf_names)) { vcf_names = gsub("_.*", "", names(vcf)) }
  not_in = which(!mutations %in% vcf_names)
  add = c()
  remove = c()
  for (i in not_in) {
    if ((i+1) %in% not_in & grepl("^rs", mutations[i]) & grepl("COSM", mutations[i+1])) {
      add = paste0(mutations[i], ';', mutations[i+1])
      remove = c(remove, i, i+1)
    }
  }
  if (length(add) > 0) {
    c(mutations[-remove], add)
  } else {
    mutations
  }
}

# add the info of the samples per branch to the tree  -------------------------------
add_sample_info <- function(tree) {
  ape_tree = tree@phylo
  # get initial edge data
  df = cbind.data.frame(node1 = ape_tree$edge[,1], 
                        node2 = ape_tree$edge[,2], 
                        tip.label = replace(rep(NA, length(ape_tree$edge.length)),
                                            ape_tree$edge[,2] <= length(ape_tree$tip.label),
                                            ape_tree$tip.label[ape_tree$edge[ape_tree$edge[,2] <= length(ape_tree$tip.label), 2]])
  )
  # add samples per edge
  df$samples = sapply(df$node2, function(node2) {
    Env <- new.env(parent = emptyenv())
    Env <- list2env(list(this_tip_labels = vector()), envir = Env)
    check_label_or_add(df = df, subnode = node2, Env = Env)
    out = paste(Env$this_tip_labels, collapse = "|")
  })
  # determine n_samples and node_id
  df = df %>%
    dplyr::select(-c(node1)) %>%
    dplyr::rename(node = node2) %>%
    mutate(n_samples = lengths(strsplit(samples, "\\|")),
           node_id = paste0(.$node, "_n", n_samples))
  # add data to the tree@data
  match = match(tree@data$node, df$node)
  columns_to_add = c('samples', 'n_samples', 'node_id', 'tip.label')
  tree@data = as_tibble(cbind(tree@data[ ,setdiff(colnames(tree@data), columns_to_add)] , df[match,columns_to_add]))

  # if root branches are not present, add them
  not_in_tree = setdiff(df$node, tree@data$node)
  same_columns = intersect(colnames(df), colnames(tree@data))
  new_columns = setdiff(colnames(tree@data), colnames(df))
  if (length(not_in_tree) > 0) {
    for (i in not_in_tree) {
      new_entry = df[df$node == i, same_columns]
      add = matrix(NA, nrow = 1, ncol = length(new_columns), dimnames = list('', new_columns))
      tree@data = rbind(tree@data, cbind.data.frame(new_entry, add))
    }
  }
  return(tree)
}

# a assistent function for determining the samples per branch --------------------------------
get_sub_nodes <- function(df, node, Env) {
  subnodes = df$node2[df$node1 == node]
  for(subnode in subnodes) {
    check_label_or_add(df, subnode, Env)
  }
}
# another assistent function for determining the samples per branch --------------------------------
check_label_or_add <- function(df, subnode, Env) {
  if (!is.na(df$tip.label[df$node2 == subnode])) {
    Env$this_tip_labels <- c(Env$this_tip_labels, df$tip.label[df$node2 == subnode])
  } else {
    get_sub_nodes(df, subnode, Env) 
  }
}

# in the most simple way, assign mutations to branches (do all expected cells have it, and none of the others?) ------------------
assign_muts_to_branches_simple <- function(vcf, tree, add1missing = TRUE, ptato_grl = NA) {
  # get mutation VAF
  vcf_names = gsub("_.*", "", names(vcf))
  vaf_all = vcf@assays@data$VAF %>% apply(2, unlist)
  vaf_all[is.na(vaf_all)] = 0
  colnames(vaf_all) = samples(header(vcf))
  
  # make tree table -> above 0 = present
  tt = lapply(colnames(vaf_all), function(cn) {
    ifelse(vaf_all[ ,cn] > 0, 1, 0)
  }) %>% do.call(cbind, .) %>% `colnames<-`(colnames(vaf_all))
  tt = tt[ ,tree@phylo$tip.label]
  
  # assign mutations that fit exactly to each branch
  dsct = sapply(tree@data$samples, function(label) {
    samples = strsplit(label, "\\|")[[1]]
    ifelse(colnames(tt) %in% samples, 1, 0)
  }) %>% 
    t %>%
    `colnames<-`(colnames(tt)) %>%
    `rownames<-`(tree@data$node_id)
  wbranch <- apply(tt, 1, function(trow) {
    same_bool = apply(dsct, 1, function(drow) {
      all.equal(trow, drow) == TRUE
    })
    rownames(dsct)[same_bool]
  })
  
  # include if can be fitted to a single branch if assuming that 1 is missing
  if (add1missing) {
    wbranchmiss <- apply(tt, 1, function(t_numb) {
      onemissedfit = apply(dsct, 1, function(d_numb) {
        sum(t_numb) == (sum(d_numb) - 1) &
          length(setdiff(names(t_numb)[t_numb == 0], names(d_numb)[d_numb == 0])) == 1 &
          length(intersect(names(t_numb)[t_numb == 1], names(d_numb)[d_numb == 1])) == (sum(d_numb) - 1)
      })
      rownames(dsct)[onemissedfit]
    })
    wbranch[lengths(wbranch) == 0 & lengths(wbranchmiss) == 1] = unlist(wbranchmiss[lengths(wbranch) == 0 & lengths(wbranchmiss) == 1])
  }
  
  # label non-hits and return
  for (n in which(lengths(wbranch) == 0)) { wbranch[[n]] = "NOBRANCH" }
  wbranch = unlist(wbranch)
  message('keeping ', sum(wbranch != "NOBRANCH"), ' from ', length(wbranch), ' mutations that fit the tree')
  
  if (length(ptato_grl) == 1 && is.na(ptato_grl)) {
    return(wbranch)
  } else {
    # filter PTATO if present
    n1_samps = setNames(tree@data$node_id, tree@data$tip.label) %>%
      grep("n1", ., value = T)
    n1_samps = n1_samps[grepl("PTA", names(n1_samps))]
    
    for (br_name in n1_samps) {
      pta_samp = names(n1_samps)[n1_samps == br_name]
      sc_muts = which(wbranch == br_name) # single cell
      sc_gr = granges(vcf[sc_muts])
      seqlevels(sc_gr) = paste0('chr', seqlevels(sc_gr))
      ptato_gr = ptato_grl[[pta_samp]]
      which_muts_keep = queryHits(findOverlaps(sc_gr, ptato_gr)) %>% unique %>% sort
      remove_single_muts = !seq(length(sc_gr)) %in% which_muts_keep
      message('filtering out ', sum(remove_single_muts), ' of ', length(remove_single_muts), ' non-PTATO muts from ', pta_samp)
      wbranch[sc_muts][remove_single_muts] = "NOBRANCH"
    }
    return(wbranch)
  }
}

# load_tree_with_support --------------------------------------------
load_tree_with_info <- function(dir, outgr = "NONE", prefix = NULL, 
                               mutation_soure = 'cellphy', 
                               vcf = NULL, ptato_grl = NA, cellphy_rm_non1 = TRUE,
                               norm_pres_max = 0.95,  high_frac_min = 0.05, min_frac_all = 0.5) {
  # define the input files
  if (is.null(prefix)) {
    files = list.files(dir)
    starttree = grep(".Tree.raxml.startTree", files, value = T)
    prefix = gsub("Tree.raxml.startTree", "", starttree)
  }
  treef = paste0(dir, '/', prefix, 'Support.raxml.mutationMapTree')
  trees = paste0(dir, '/', prefix, 'Support.raxml.support')
  mutf = paste0(dir, '/', prefix, 'Support.raxml.mutationMapList')
  
  # load tree
  message('reading and processing tree...')
  tree = read_ggtree(treef)
  data = read.table(mutf, head=F, fill=T, col.names=c("edgeID", "NumberOfMutations", "MutationList"))
  gene_names = NULL
  
  # add samples
  tree = add_sample_info(tree)
  
  # add bootstraps
  tree_supp = read_ggtree(trees)
  all.equal(tree@phylo$edge, tree_supp@phylo$edge)
  tree@data$n_boot = NA
  interm_nodes = (length(tree@phylo$tip.label) + 1):max(tree@phylo$edge)
  tree@data$n_boot[match(interm_nodes, tree@data$node)] = tree_supp@phylo$node.label
  
  # add the mutationList and branch_lengths
  tree@data$mutation_names = NA
  tree@data$branch_length = NA
  message('fitting mutations to the tree...')
  if (mutation_soure == 'cellphy') {
    message('cellphy method selected')
    for (i in 1:length(tree@data$N)) {
      if (tree@data$N[i] != -1) {
        names = strsplit(subset(data[,3], data$edgeID == tree@data$N[i]), ',')[[1]]
        names = merge_cosmic_muts(names, vcf = vcf)
        names = grep("^Un|^chrUn|^MT|^chrMT", names, value = TRUE, invert = TRUE)
        tree@data$branch_length[i] = length(names)
        tree@data$mutation_names[i] = paste(names, collapse = "|")
      } else {
        tree@data$mutation_names[i] = ""
        tree@data$branch_length[i] = 0
      }
    }
    
    
    # check if "OUTGROUP" branch mutations are correct, or belong to the other branch
    outgr_branch = which(sapply(tree@data$samples, function(s) { all(strsplit(s, "\\|")[[1]] == outgr) }))
    tree = root_tree(tree, outgr, add_info = TRUE)
    
    # First select the samples, then do the split 
    normal_muts = strsplit(tree@data$mutation_names, "\\|")[[outgr_branch]]
    print("Did it pass my error?")
    if (length(normal_muts) > 0) {
      vcf_names = gsub("_.*", "", names(vcf))
      normal_vcf = vcf[vcf_names %in% normal_muts]
      
      # extract VAF and calculate how frequent in outgroup
      vafs = normal_vcf@assays@data$VAF %>% apply(2, unlist)
      vafs[is.na(vafs)] = 0
      colnames(vafs) = samples(header(vcf))
      normal_pres_frac = sum(vafs[ ,outgr] > 0)/nrow(vafs)
      n_samp = rowSums(vafs > 0)
      high_samp_frac = sum(n_samp > 2 & vafs[ ,outgr] == 0)/nrow(vafs)
      other_vaf = vafs[ ,setdiff(colnames(vafs), outgr)]
      n_others = rowSums(other_vaf > 0)
      frac_all_others = sum(n_others == ncol(other_vaf))/nrow(vafs)
      other_branch = which(sapply(strsplit(tree@data$samples, '\\|'), function(s) all(colnames(other_vaf) %in% s)))
      
      if (high_samp_frac > high_frac_min && normal_pres_frac < norm_pres_max && frac_all_others > min_frac_all) {
        message("moving mutations from outgroup to main branch")
        message("this often happens when all cells in the sample share mutations")
        message(format(normal_pres_frac * 100, digits = 5), "% of muts were found in outgroup")
        message(format(high_samp_frac * 100, digits = 5), "% of muts were in more than two mutations other than the outgroup")
        message("on average found in ", format(mean(n_others), digits = 2), " out of ", ncol(other_vaf), " other samples")
        tree@data$mutation_names[other_branch] = tree@data$mutation_names[outgr_branch]
        tree@data$branch_length[other_branch] = tree@data$branch_length[outgr_branch]
        tree@data$mutation_names[outgr_branch] = NA
        tree@data$branch_length[outgr_branch] = 0
      }
    }
    
    # filter PTATO if present
    if (length(ptato_grl) > 1 && any(class(ptato_grl) == "CompressedGRangesList")) {
      n1_samps = setNames(tree@data$node_id, tree@data$tip.label) %>%
        grep("n1", ., value = T)
      n1_samps = n1_samps[grepl("PTA", names(n1_samps))]
      vcf_names = gsub("_.*", "", names(vcf))
      for (br_name in n1_samps) {
        pta_samp = names(n1_samps)[n1_samps == br_name]
        muts = strsplit(tree@data$mutation_names, "\\|")[[which(tree@data$node_id == br_name)]]
        sc_muts = which(vcf_names %in% muts) # single-cell
        sc_gr = granges(vcf[sc_muts])
        seqlevels(sc_gr) = paste0('chr', seqlevels(sc_gr))
        ptato_gr = ptato_grl[[pta_samp]]
        which_muts_keep = queryHits(findOverlaps(sc_gr, ptato_gr)) %>% unique %>% sort
        remove_single_muts = !seq(length(sc_gr)) %in% which_muts_keep
        message('filtering out ', sum(remove_single_muts), ' of ', length(remove_single_muts), ' non-PTATO muts from ', pta_samp)
        tree@data$mutation_names[which(tree@data$node_id == br_name)] = paste(gsub("_.*", "", names(vcf))[sc_muts][!remove_single_muts], collapse = "|")
      }
    }
    
    # remove endbranch mutations that are found in more than 1 sample
    if (cellphy_rm_non1) {
      n1_samps = setNames(tree@data$node_id, tree@data$tip.label) %>%
        grep("n1", ., value = T)
      vcf_names = gsub("_.*", "", names(vcf))
      vaf = vcf@assays@data$VAF %>% 
        apply(2, unlist)
      vaf[is.na(vaf)] = 0
      rownames(vaf) = vcf_names
      for (br_name in n1_samps) {
        samp = names(n1_samps)[n1_samps == br_name]
        muts = strsplit(tree@data$mutation_names, "\\|")[[which(tree@data$node_id == br_name)]]
        muts = intersect(muts, vcf_names)
        if (length(muts) > 0) {
          if (length(muts) > 1 || !is.na(muts)) {
            remove = rowSums(vaf[muts, , drop = F] > 0) > 1
            message('filtering out ', sum(remove), ' of ', length(remove), ' non-individual muts from ', samp)
            tree@data$mutation_names[which(tree@data$node_id == br_name)] = paste(muts[!remove], collapse = "|")
          }
        }
      }
    }
    return(tree)
  } else if (mutation_soure == 'simple' & !is.null(vcf)) {
    # do the same for the "fitting" mutations
    message('VAF method selected')
    tree = root_tree(tree, outgr, add_info = TRUE)
    wbranch = assign_muts_to_branches_simple(vcf, tree, ptato_grl = ptato_grl)
    vcf_names = gsub("_.*", "", names(vcf))
    tree@data$mutation_names = sapply(tree@data$node_id, function(n_node) {
      paste(vcf_names[wbranch == n_node], collapse =  "|")
    })
    tree@data$branch_length = lengths(strsplit(tree@data$mutation_names, "\\|"))
  } else {
    stop('no valid mutation_source, or no VCF provided')
  }
  return(tree)
}