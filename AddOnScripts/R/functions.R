# prepare_tree ----------------------------------------
#' Adds branch_ids to a tree
#'
#' @param tree a \code{\link[treeio]{treedata}} object created by "cellPhyWrapper"
#'
#' @return a tree, with a \code{branch_id} column for branch identification
#' @export
prepare_tree <- function(tree) {
  ltrs = c(LETTERS, letters) # must be longer
  symbols = lapply(c('', 2:10), function(x) paste0(ltrs, x)) %>% do.call(c, .)
  tree@data$branch_id = symbols[seq(nrow(tree@data))]
  return(tree)
}

# prepare_tree ----------------------------------------
#' change the tip.labels of your tree
#'
#' @param tree a \code{\link[treeio]{treedata}} object created by "cellPhyWrapper"
#' @param id_change_tb a data.frame with 2 columns, the first the current tip.labels (to change), and the second, the new corresponding IDs
#' @return a tree, with a \code{branch_id} column for branch identification
#' @export
change_tip_labels <- function(tree, id_change_tb) {
  stopifnot(any(class(id_change_tb) %in% c('matrix', 'data.frame', 'tbl_df')))
  message('using column ', colnames(id_change_tb)[1], 'as old IDs, and ', colnames(id_change_tb)[2], 'as new IDs')
  stopifnot(any(tree@phylo$tip.label %in% pull(id_change_tb, 1)))
  if (!all(tree@phylo$tip.label %in% pull(id_change_tb, 1))) {
    warning('Warning: only replacing ', sum(tree@phylo$tip.label %in% pull(id_change_tb, 1)), ' tip labels')
  }
  search_for = intersect(tree@phylo$tip.label, id_change_tb$old)
  tree@phylo$tip.label[match(search_for, tree@phylo$tip.label)] = id_change_tb$new[match(search_for, id_change_tb$old)]
  tree@data$tip.label[match(search_for, tree@data$tip.label)] = id_change_tb$new[match(search_for, id_change_tb$old)]
  return(tree)
}

# find_common_sample_string ----------------------------------------
#' find the longest common substring in the sample IDs
#'
#' @inheritParams prepare_tree
#'
#' @return the longest common substring in the sample IDs
find_common_sample_string <- function(tree) {
  tipl = tree@data$tip.label
  tipl = tipl[!is.na(tipl)]
  for (n_tip in 2:length(tipl)) {
    prev_tip = n_tip -1
    if (prev_tip == 1) { output_prev = tipl[1] }
    sb <- stri_sub(tipl[n_tip], 1, 1:nchar(tipl[n_tip]))
    ## extract them from 'a' if they exist
    output_prev <- stri_extract_all_coll(output_prev, sb) %>% lapply('[[', 1) %>% unlist %>% na.omit
    ## match the longest one
    output_prev = output_prev[which.max(nchar(output_prev))]
  }
  common_str = output_prev
}
# add_x_axis ----------------------------------------
#' add an x-axis to a ggtree
add_x_axis <- function() {
  theme(axis.line.x = element_line(),
        axis.text.x = element_text(),
        axis.ticks.x = element_line()
  )
}
# extract_vcf_per_branch ----------------------------------------
#' Title
#'
#' @inheritParams prepare_tree
#' @param vcf the VCF used to build the tree as a \code{CollapsedVCF}, i.e., imported by \code{\link[VariantAnnotation]{readVcf}}
#' @param make_granges a boolean: whether to output a list of \code{CollapsedVCFs} (\code{FALSE}) or a \code{\link[GenomicRanges]{GRangesList}} (\code{TRUE})
#' @param ref_genome a character string: the reference BSgenome
#' @param incl_chroms a numeric: which chromosomes to include from the ref_genome
#'
#' @return the mutations per branch, either as a list of \code{CollapsedVCFs} or a \code{\link[GenomicRanges]{GRangesList}}
#' @export
extract_vcf_per_branch <- function(tree, vcf, make_granges = FALSE, ref_genome = NA, incl_chroms = 1:24) {
  if (!'branch_id' %in% colnames(tree@data)) { stop("please run prepare_tree before this function") }
  ref = get(ref_genome)
  chroms = seqlevels(ref)[incl_chroms]
  if (grepl("chr", chroms[1]) & !grepl("chr", seqlevels(vcf)[1])) { seqlevels(vcf) = paste0('chr', seqlevels(vcf)) }
  if (!grepl("chr", chroms[1]) & grepl("chr", seqlevels(vcf)[1])) { seqlevels(vcf) = str_remove(seqlevels(vcf), 'chr') }
  seqlevels(vcf, pruning.mode = 'coarse') = chroms

  mut_name_list = strsplit(tree@data$mutation_names, "\\|")
  vcf_names = str_remove(names(vcf), "_.*")
  # check if all mutation are found in the tree
  no_overlap = length(setdiff(unique(unlist(mut_name_list)), vcf_names))
  no_ov_perc =  (no_overlap/sum(lengths(mut_name_list)))
  if (no_ov_perc > 0.02) { warning(paste(no_ov_perc*100, "percent of your tree mutations is not found in your VCF. Are you sure you loaded the correct VCF?")) }
  # make vcf
  branch_vcf_list = lapply(mut_name_list, function(mut_names) {
    vcf[vcf_names %in% mut_names]
  })
  names(branch_vcf_list) = tree@data$branch_id

  # make GRL
  if (make_granges) {
    branch_vcf_list = convert_vcf_to_granges(branch_vcf_list, ref_genome, incl_chroms)
  }
  return(branch_vcf_list)
}
# convert_vcf_to_granges ----------------------------------------
#' convert a list of VCFs to a GRangesList
#'
#' @inheritParams extract_vcf_per_branch
#' @param branch_vcf_list a list of \code{CollapsedVCFs} to convert to a \code{\link[GenomicRanges]{GRangesList}}
#'
#' @return a GRangesList
#' @export
convert_vcf_to_granges <- function(branch_vcf_list, ref_genome, incl_chroms = 1:24) {
  if (!is.na(ref_genome)) {
    # load ref genome
    ref = get(ref_genome)
    chroms = seqlevels(ref)[incl_chroms]
    branch_grl = GRangesList(lapply(branch_vcf_list, function(v) {
      gr = granges(v)
      genome(gr) = genome(ref)[1]
      # make seq levels correct
      if (grepl("chr", chroms[1]) & !grepl("chr", seqlevels(gr)[1])) { seqlevels(gr) = paste0('chr', seqlevels(gr)) }
      if (!grepl("chr", chroms[1]) & grepl("chr", seqlevels(gr)[1])) { seqlevels(gr) = str_remove(seqlevels(gr), 'chr') }
      seqlevels(gr, pruning.mode = 'coarse') = chroms
      return(gr)
    }))
  } else {
    stop("no ref_genome given, please specify")
  }
  return(branch_grl)
}
# add_info_to_tree ----------------------------------------
#' add additional branch/node informationi to the tree
#'
#' @inheritParams prepare_tree
#' @param by the column in \code{tree@@data} to join by
#' @param ... info to add to the tree.
#' Either data.frames/tibbles with a column specified in \code{by}, containing the branch ids (as present in \code{tree@@data[ ,by]})
#' or a numeric/factor/character named with the branch ids found in \code{tree@@data[ ,by]}
#'
#' @return a tree with the additional stored information
#' @export
add_info_to_tree <- function(tree, by = 'branch_id', ...) {
  new_info_list = list(...)
  if(is.null(names(new_info_list)) || any(names(new_info_list) == '')) {
    stop("please provide a name for every variable: 'tree <- add_info_to_tree(tree, id = patient_id, number = 1:23)', etc.")
  }
  by_column_tree = pull(tree@data, by)
  # for each variable, check info and add
  for (n_info in seq(length(new_info_list))) {
    info = new_info_list[[n_info]]
    info_name = names(new_info_list)[[n_info]]
    # data.frame-like
    if (any(class(info) %in% c('tibble', 'data.frame', 'matrix'))) {
      if (by %in% colnames(info)) {
        by_column_new = pull(info, by)
      } else {
        by_column_new = rownames(info)
        info[ ,by] = rownames(info)
        warning(paste0("using rownames of '", info_names, "' for merging. It is better to add column '", by, "'."))
      }
      if (length(intersect(by_column_new, by_column_tree)) == 0) {
        stop(paste0("none of the values of the '", by, "' column of '", info_name, "', overlap with tree@data$", by))
      } else if (length(setdiff(by_column_new, by_column_tree)) > 0) {
        stop(paste0("some of the values of the '", by, "' column of '", info_name, "', are not found in tree@data$", by, ". Please filter your input data"))
      } else {
        keep_columns = setdiff(colnames(tree), setdiff(colnames(info), by))
        tree@data = left_join(tree@data[ ,keep_columns], info)
      }
    # vector-like
    } else if (any(class(info) %in% c('character', 'numeric', 'factor'))) {
      if (is.null(names(info)) || length(intersect(names(info), by_column_tree)) == 0) {
        stop(paste0("none of the names of '", names(new_info_list)[n_info], "', overlap with tree@data$", by, ". Please add names or change the 'by' parameter"))
      } else if (length(setdiff(names(info), by_column_tree)) > 0) {
        stop(paste0("some of the names of '", names(new_info_list)[n_info], "', are not found in tree@data$", by, ". Please filter your input data"))
      } else {
        if (!setequal(names(info), by_column_tree)) {
          info2 = setNames(rep(NA, length(by_column_tree)), by_column_tree)
          info2[names(info)] = unname(info)
          info = info2
        }
        tree@data[ ,info_name] = info[match(by_column_tree, names(info))]
      }
    }
  }
  return(tree)
}
# check_0_mm ----------------------------------------
#' remove columns with 0 (or another defined number of) mutations from a mutational matrix
#'
#' @param mut_matrix a mutational matrix with samples in the columns and mutations in the rows
#' @param remove_min columns with equal or fewer mutations that this number will be removed
#'
check_0_mm <- function(mut_matrix, remove_min = 0) {
  if (any(colSums(mut_matrix) <= remove_min)) {
    miss_colnames = colnames(mut_matrix)[colSums(mut_matrix) <= remove_min]
    warning(paste('removing the following profiles with equal or fewer than', remove_min, 'mutations:', paste0(miss_colnames, collapse = ', ')))
    mut_matrix = mut_matrix[ ,!colnames(mut_matrix) %in% miss_colnames]
  }
  return(mut_matrix)
}
# fit_to_signatures_strict_tree ----------------------------------------
#' A wrapper that runs and filters the output of \code{\link[MutationalPatterns]{fit_to_signatures_strict}}
#'
#' @param mut_matrix a mutational matrix with samples in the columns and mutations in the rows
#' @param signatures a matrix with signature names in the columns and mutations in the rows
#' @param max_delta see \code{\link[MutationalPatterns]{fit_to_signatures_strict}}
#' @inheritParams check_0_mm
#'
#' @return a contribution matrix with signatures in the rows and samples in the columns
#' @export
fit_to_signatures_strict_tree <- function(mut_matrix, signatures, max_delta = 0.01, remove_min = 20) {
  mut_matrix = check_0_mm(mut_matrix, remove_min)
  contri = fit_to_signatures_strict(mut_matrix, signatures, max_delta = max_delta)$fit_res$contribution
  contri = contri[rowSums(contri) > 0, ]
}
# check_reconstructed_cosine ----------------------------------------
#' plotting whether a reconstructed profile (after signature refitting) is similar enough to the original profiel
#'
#' @param contribution a contribution matrix with signatures in the rows and samples in the columns
#' @inheritParams check_0_mm
#' @inheritParams fit_to_signatures_strict_tree
#' @inheritParams prepare_tree
#'
#' @return a plot comparing the number of mutation per profile to the cosine similarity between the reconstructed vs the original profile
#' @export
check_reconstructed_cosine <- function(contribution, mut_matrix, signatures, tree) {
  recon = signatures %*% contribution
  recon_cosine = diag(cos_sim_matrix(recon, mut_matrix[ ,colnames(contribution)]))
  tibble(
    recon = recon_cosine,
    n_mut = tree@data$branch_length[match(colnames(contribution), tree@data$branch_id)],
    branch = colnames(contribution)
  ) %>% ggplot(aes(x = log10(n_mut + 1), y = recon)) +
    geom_point() +
    geom_hline(yintercept = 0.85, lty = 3) +
    geom_vline(xintercept = log10(201), lty = 3)+
    ggrepel::geom_text_repel(data = ~ subset(.x, n_mut > 200 & recon < 0.85), aes(label = branch)) +
    labs(x = 'log10 # of mutations', y = 'cosine sim. reconstructed vs original profile')
}
# add_tree_drivers ---------------------------------------------
#' add information on (potential) driver mutations to a tree
#'
#' @inheritParams prepare_tree
#' @param branch_vcf_list a list of \code{CollapsedVCFs}, one for each branch
#' @param add_mod_high a boolean: whether or not to include MODERATE/HIGH impact mutations
#' @param add_cosmic a boolean: whether or not to include mutations with a COSMIC ID (i.e., at least 1 sample registered in the COSMIC database)
#' @param blacklist a (regex) string which will be excluded from the driver list, e.g. "^MUC.*" to remove mutation in MUC genes
#' @param include_list a vector of gene names to include in the driver list (ignoring all other genes)
#' @return a tree with a new column "drivers" containting potential cancer drives
#' @export
add_tree_drivers <- function(tree, branch_vcf_list, add_mod_high = TRUE, add_cosmic = FALSE, blacklist = NULL, include_list = NULL) {
  drivers = sapply(branch_vcf_list, function(v) {
    ann = v@info$ANN
    driv = c()
    exonic = lengths(ann) > 0
    if (sum(exonic) > 0) {
      ann_ex = sapply(ann[exonic], '[[', 1)
      ann_split = strsplit(ann_ex, '\\|')

      slct =  c()
      mod_high = sapply(ann_split, '[[', 3) %in% c("MODERATE", "HIGH")
      if (add_mod_high) { slct = c(slct, which(mod_high)) }
      cosmic = grep("COSM", names(v[exonic]))
      if (add_cosmic) { slct = c(slct, cosmic) }
      if (length(slct) > 0) {
        ann_split_slct = ann_split[slct]
        type = sapply(ann_split_slct, '[[', 2) %>%
          str_replace('missense_variant', 'miss') %>%
          str_replace('stop_gained', 'stop') %>%
          str_replace('structural_interaction_variant', 'miss') %>%
          str_replace('frameshift', 'fs') %>%
          str_remove('&intron_variant') %>%
          str_replace_all('splice.*', 'splice')
        gene = sapply(ann_split_slct, '[[', 4)
        driv = c(driv, paste0(gene, "_", type))
      } else {
        driv = ''
      }
    } else {
      driv = ''
    }
    if (!is.null(blacklist)) { driv = driv[!driv %in% blacklist] }
    if (!is.null(include_list)) { driv = driv[driv %in% include_list] }
    driv = paste(driv, collapse = "\n")
  })
  tree@data$drivers = drivers
  return(tree)
}

# find_gene  ----------------------------------------
#' Find mutations in certain genes in a list of
#'
#' @inheritParams add_tree_drivers
#' @param find_gene
#'
#' @return a string of matching mutations in the genes, named by the branch_id
#' @export
find_gene <- function(branch_vcf_list, find_gene) {
  found = sapply(branch_vcf_list, function(v) {
    ann = v@info$ANN
    driv = c()
    exonic = lengths(ann) > 0
    if (sum(exonic) > 0) {
      ann_ex = sapply(ann[exonic], '[[', 1)
      if (any(grepl(find_gene, ann_ex))) {
        ann_split = strsplit(ann_ex[grepl(find_gene, ann_ex)], '\\|')
        names = names(v)[exonic][grepl(find_gene, ann_ex)]
        paste0(names, "|", sapply(ann_split, '[[', 2), "_", sapply(ann_split, '[[', 4), collapse = "|")
      } else return(NA)
    } else return(NA)
  })
  found = found[!is.na(found)]
  if (length(found) == 0) {
    warning(paste("no branches with mutations in", find_gene, "were found"))
  } else {
    return(found)
  }
}
# add_contribution ----------------------------------------
#' add a contribution matrix to a tree
#'
#' @inheritParams prepare_tree
#' @param contribution a contribution matrix with signatures in the rows and samples in the columns
#' @inheritParams fit_to_signatures_strict_tree
#' @inheritParams check_0_mm
#' @param ... parameters passed on to \code{\link[MutationalPatterns]{fit_to_signatures_strict}}
#'
#' @return a tree with a new data column for each signature in "contribution"
#' @export
add_contribution <- function(tree, contribution = NULL, signatures = NULL, mut_matrix = NULL, ...) {
  if((is.null(contribution) && is.null(signatures)) ||
     (!is.null(contribution) && !is.null(signatures))) stop("please provide either a contribution matrix, or the signatures you want to fit on")
  if(!is.null(signatures) && is.null(mut_matrix)) stop("please also provide a mutation matrix of the branches")
  # fit signatures if not yet done
  if(!is.null(signatures)) { contribution = fit_to_signatures_strict_tree(mut_matrix, signatures, ...) }

  # add contribution to the tree
  if(!is.null(contribution)) {
    contribution = t(t(contribution)/colSums(contribution))
    add = as.data.frame(t(contribution))
    branches_missing = setdiff(tree@data$branch_id, colnames(contribution))
    add = rbind(add, as.data.frame(matrix(NA, nrow = length(branches_missing), ncol = ncol(add), dimnames = list(branches_missing, colnames(add))))) %>%
      rownames_to_column('branch_id') %>%
      as_tibble
    overlap = setdiff(intersect(colnames(tree@data), colnames(add)), 'branch_id')
    tree@data = full_join(tree@data[ ,setdiff(colnames(tree@data), overlap)], add, by = 'branch_id')
  }
  return(tree)
}
# grouped_tree_contributions ----------------------------------------
#' add the contribution matrix columns from a group together
#'
#' @inheritParams prepare_tree
#' @param group a column in \code{tree@@data} that specifies the group by which to summarize the contribution matrix
#' @param contribution a contribution matrix with signatures in the rows and samples in the columns
#' @inheritParams fit_to_signatures_strict_tree
#' @inheritParams check_0_mm
#' @param ... parameters passed on to \code{\link[MutationalPatterns]{fit_to_signatures_strict}}
#'
#' @return a contribution matrix with one column per group, containing the sum of the mutation numbers in each group
#' @export
group_tree_contributions <- function(tree, group, contribution = NULL, signatures = NULL, mut_matrix = NULL, ...) {
  if((is.null(contribution) && is.null(signatures)) ||
     (!is.null(contribution) && !is.null(signatures))) stop("please provide either a contribution matrix, or the signatures you want to fit on")
  if(!is.null(signatures) && is.null(mut_matrix)) stop("please also provide a mutation matrix of the branches")
  # fit signatures if not yet done
  if(!is.null(signatures)) { contribution = fit_to_signatures_strict_tree(mut_matrix, signatures, ...) }
  if(!is.null(contribution)) {
    contri = as.data.frame(contribution)
    branches_missing = setdiff(tree@data$branch_id, colnames(contribution))
    contri = cbind(contri, as.data.frame(matrix(0, ncol = length(branches_missing), nrow = nrow(contri),
                                                dimnames = list(rownames(contri), branches_missing))))
    group = pull(tree@data, group)
    if (!any(class(group) == 'factor')) { group = factor(group) }
    grouped_contri = lapply(levels(group), function(gr) { rowSums(contri[ ,tree@data$branch_id[group == gr],drop=FALSE]) }) %>%
      do.call(cbind, .) %>% `colnames<-`(levels(group))
  }
  return(grouped_contri)
}

# group_tree_mut_matrix ----------------------------------------
#' add the mutational matrix columns from a group together
#'
#' @inheritParams prepare_tree
#' @param group a column in \code{tree@@data} that specifies the group by which to summarize the mutational matrix
#' @inheritParams check_0_mm
#'
#' @return a mutational matrix with one column per group, containing the sum of the mutation numbers in each group
#' @export
group_tree_mut_matrix <- function(tree, group, mut_matrix) {
  mm = as.data.frame(mut_matrix)
  branches_missing = setdiff(tree@data$branch_id, colnames(mm))
  mm_full = cbind(mm, as.data.frame(matrix(0, ncol = length(branches_missing), nrow = nrow(mm),
                                              dimnames = list(rownames(mm), branches_missing))))
  group = pull(tree@data, group)
  if (!any(class(group) == 'factor')) { group = factor(group) }
  grouped_mm = lapply(levels(group), function(gr) { rowSums(mm_full[ ,tree@data$branch_id[group == gr], drop=FALSE]) }) %>%
    do.call(cbind, .) %>% `colnames<-`(levels(group))
  return(grouped_mm)
}

# plot_gg_tree ----------------------------------------
#' plot a tree
#'
#' @inheritParams prepare_tree
#' @param common_name a string to substract from the sample names and use to plot above the plot. (default: the longest common sample substring)
#' @param add_x a boolean: whether to plot an x-axis
#' @param add_tip_label a boolean: whether to plot the sample names at the tips of the leaves (i.e., end branches)
#' @param add_bootstrap a boolean: whether to plot the bootstrap information on the nodes
#' @param add_branch_length a boolean: whether to plot the branch length numbers
#' @param add_branch_id a boolean: whether to plot the branch IDs
#' @param branch_text_param a character: which additional text to plot at each branch.
#' Should be one of the column names of \code{tree@@data}.
#' @param branch_text_size the text size of the 'branch_text_param' parameter
#' @param only_shared_branches a boolean: whether to only plot the shared branches
#' @param plot_margin the plot margin (default: \code{unit(c(10,100,10,10), 'points')})
#' @param branch_color_param a character: how to color each branch.
#' Should be one of the column names of \code{tree@@data}.
#' @param branch_colors the colors used for coloring the branches
#' @param add_title a character with the title to plot. Or a boolean: whether to plot the 'common_name' as a title
#' @param legend_pos a specification of the legend position, e.g. "bottom". See \code{\link[ggplot2]{theme}}, "legend.position"
#'
#' @return a ggplot
#' @export
plot_gg_tree <- function(tree, common_name = NA, add_x = TRUE,
                         add_tip_label = 'short_tip_lab', add_bootstrap = FALSE,
                         add_branch_length = FALSE, add_branch_id = FALSE,
                         branch_text_param = NULL, branch_text_size = 3,
                         only_shared_branches = FALSE, plot_margin = unit(c(10,100,10,10), 'points'),
                         branch_color_param = NULL, branch_colors = NA, add_title = TRUE, legend_pos = NULL,
                         branchid_nudge=0, branchlength_nudge=0) {
  # adjust labels
  if (is.na(common_name)) { common_name = find_common_sample_string(tree) }
  if (!is.null(common_name)) {
    tree@data$short_tip_lab = str_remove(tree@data$tip.label, common_name)
  } else {
    tree@data$short_tip_lab = tree@data$tip.label
  }

  # filter branches for plotting
  if (only_shared_branches) {
    tree@data[tree@data$n_samples == 1, setdiff(colnames(tree@data), 'node')] = NA
    add_tip_label = NA
    plot_margin = unit(c(10,10,10,10), 'points')
  }

  # plot tree
  if (!is.null(branch_color_param)) {
    gg = ggtree(tree, branch.length = 'branch_length', mapping = aes(color = !!sym(branch_color_param))) +
      coord_cartesian(clip = 'off') +
      theme(plot.margin = plot_margin)
    if (length(branch_colors) == 1 && is.na(branch_colors)) {
      if (any(c('factor', 'character') %in% class(pull(tree@data, branch_color_param)))) {
        classes = length(unique(tree@data[ ,branch_color_param]))
        if (classes <= 22) { cols = dist_cols } else
          if (classes <= 31) { cols = dist_cols31 } else
            if (classes <= 50) { cols = dist_cols31 } else
              { cols = scales::hue_pal()(classes) }
        gg = gg + scale_color_manual(values = cols)
      } else {
        gg = gg + scale_color_gradientn(colors = grad_cols)
      }
    } else {
      if (any(c('factor', 'character') %in% class(tree@data[ ,branch_color_param]))) {
        gg = gg + scale_color_gradientn(colors = branch_colors)
      } else {
        gg = gg + scale_color_manual(values = branch_colors)
      }
    }
  } else {
    gg = ggtree(tree, branch.length = 'branch_length') +
      coord_cartesian(clip = 'off') +
      theme(plot.margin = plot_margin)
  }

  # add any info
  if (add_branch_id) gg = gg + geom_text(aes(x=branch+branchid_nudge, label=branch_id), nudge_y = 0.5, size = 2, na.rm = T)
  if (add_branch_length) gg = gg + geom_text(aes(x=branch+branchlength_nudge, label=branch_length), nudge_y = 0.5, size = 2, na.rm = T)
  if (add_bootstrap) gg = gg + geom_nodelab(aes(label = n_boot), geom = 'label', na.rm = T, fill = rgb(1,1,1,0.8), size = 2, color = 'grey70')
  if (!is.na(add_tip_label)) {
    if (any(class(add_tip_label) == 'character')) {
      gg = gg + geom_tiplab(aes(label = !!sym(add_tip_label)), size = 3, color = 'grey70')
    } else if (any(class(add_tip_label) == 'logical')) {
      if (add_tip_label) {
        gg = gg + geom_tiplab(size = 3, color = 'grey70')
      }
    }
  }
    labs(title = common_name)
  if (!is.null(branch_text_param)) gg = gg + geom_text(aes(x=branch, label = !!sym(branch_text_param)), nudge_y = 0.1, size = branch_text_size, color = 'black')
  if (add_x) gg = gg + add_x_axis()
  if (is.logical(add_title) && add_title) gg = gg + labs(title = common_name)
  if (is.character(add_title)) gg = gg + labs(title = add_title)
  if (!is.null(legend_pos)) gg = gg + theme(legend.position = legend_pos)
  return(gg)
}

# plot_gg_tree_base -----------------------------------------------------------
#' the plot_gg_tree with some default settings, like plotting the branch IDs
#'
#' @inheritParams prepare_tree
#' @inheritParams plot_gg_tree
#' @param ... parameters passed on to \code{\link{plot_gg_tree}}
#'
#' @return a ggplot
#' @export
plot_gg_tree_base <- function(tree, add_tip_label = 'short_tip_lab', add_bootstrap = TRUE,
                               add_branch_length = TRUE, add_branch_id = TRUE, ...) {
  plot_gg_tree(tree, add_tip_label = add_tip_label, add_bootstrap = add_bootstrap,
               add_branch_length = add_branch_length, add_branch_id = add_branch_id,
               ...)
}

# plot_tree_contribution -----------------------------------------------
#' plot the tree with the contribution of signature per branch
#'
#' @inheritParams prepare_tree
#' @inheritParams plot_gg_tree
#' @param signature the signature to plot. Should have been add before with \link{add_contribution}
#' @param type a character. "pie" to plot a pie chart, 'color' to color the branches, 'color_dot' to add a colored dot
#' @param pie_size the size of the piechart per branch if type = 'pie'. (default = 1)
#' @param ... parameters to pass on to \link{plot_gg_tree}
#'
#' @return a ggplot
#' @export
plot_tree_contribution <- function(tree, signature, type = 'pie', pie_size = 1, ...) {
  if (type == 'pie') {
    p = plot_gg_tree(tree, ...) +
      new_scale_color()
    add_contri_pie_to_tree(plot = p, tree = tree, signature = signature, cex = pie_size)
  } else if (type == 'color_dot') {
    plot_gg_tree(tree, ...) +
      geom_point(aes(x = branch, color = !!sym(signature)), size = pie_size * 2) +
      scale_color_gradientn(colors = grad_cols)
  } else if (type == 'color') {
    plot_gg_tree(tree, branch_color_param = signature, ...)
  } else { stop("type must be 'pie', 'color', or 'color_dot'") }
}

# plot_tree_contribution -----------------------------------------------
#' plot the tree with the contribution of signature per branch, as bars in the tree
#'
#' @inheritParams prepare_tree
#' @inheritParams plot_gg_tree
#' @inheritParams group_tree_contributions
#' @param scaling how to scale the bars. If the bars are too big/small, this can be changed to adjust them
#' @param remove_min columns with equal or fewer mutations that this number will be removed
#' @param ... parameters to pass on to \link{plot_gg_tree}
#'
#' @return a ggplot
#' @export
plot_tree_contribution_bars <- function(tree, signatures = NULL, contribution = NULL,
                                        mut_matrix = NULL, signature_colors = NULL, scaling = 1, remove_min = 20, ...) {
  if((is.null(contribution) && is.null(signatures)) ||
     (!is.null(contribution) && !is.null(signatures))) stop("please provide either a contribution matrix, or the signatures you want to fit on")
  if(!is.null(signatures) && is.null(mut_matrix)) stop("please also provide a mutation matrix of the branches")
  if(!is.null(signatures)) { contribution = fit_to_signatures_strict_tree(mut_matrix, signatures, remove_min = remove_min, ...) }

  # colors
  if (is.null(signature_colors)) {
    n_sig = nrow(contribution)
    if (n_sig <= 22) { cols = dist_cols } else
      if (n_sig <= 31) { cols = dist_cols31 } else
        if (n_sig <= 50) { cols = dist_cols31 } else
        { cols = scales::hue_pal()(n_sig) }
    signature_colors = setNames(cols, rownames(contribution))
    signature_colors = signature_colors[!is.na(names(signature_colors))]
  }
  if (is.null(names(signature_colors))) {
    names(signature_colors) = rownames(contribution)
    signature_colors = signature_colors[!is.na(names(signature_colors))]
  }

  # calculate the fractions
  fraction = contribution %>%
    t() %>% as.data.frame %>%
    rownames_to_column('node') %>%
    pivot_longer(cols = -starts_with('node'), names_to = 'sig', values_to = 'contribution')
  fraction$node = tree@data$node[match(fraction$node, tree@data$branch_id)] # use node_id, this is what the "inset" function need later

  # this part of the code has been adapted from ggtree::nodepie
  # all credits go to Guangchuang Yu: https://bioconductor.org/packages/devel/bioc/vignettes/ggtree/inst/doc/ggtree.html
  fractions_split <- split(fraction, fraction$node)
  bars = lapply(fractions_split, function(df) {
    ggplot(df, aes(x = '', y = contribution, fill = sig)) +
      geom_col() +
      coord_flip() +
      scale_fill_manual(values = signature_colors) +
      theme_void() +
      theme(legend.position = 'none')
  })
  nodes_slct = fractions_split %>% sapply(function(x) pull(x, node)[1]) %>% unname

  ### OLD CODE ###
  # determine the bar lengths
  # bar_lengths = tree@data$branch_length/max(tree@data$branch_length)
  # bar_lengths = bar_lengths[match(nodes_slct, tree@data$node)] * scaling
  # plot = plot_gg_tree(tree, ...)
  # tree_plot = ggtree::inset(tree_view = plot, insets = bars, width = bar_lengths, height = 0.05, x = 'branch')
  ### END ###

  # determine the bar lengths
  plot = plot_gg_tree(tree, ...)
  max_x <- max(layer_scales(plot)$x$range$range)
  bar_lengths = tree@data$branch_length/max_x
  bar_lengths = bar_lengths[match(nodes_slct, tree@data$node)] * scaling
  for ( i in c(1:length(bars))) {
    tree_plot = ggtree::inset(tree_view = plot, insets = bars[i], width = bar_lengths[i], height = 0.05, x = 'branch')
    plot <- tree_plot
  }

  # add legend
  fig_for_leg = ggplot(fraction, aes(x = '', y = contribution, fill = sig)) +
    geom_col() +
    scale_fill_manual(values = signature_colors)
  legend = cowplot::get_legend(fig_for_leg)

  cowplot::plot_grid(tree_plot, legend, nrow = 1, rel_widths = c(1, 0.2))
}

# add_contri_pie_to_tree ----------------------------------------------------
#' add piecharts to a pre-plotted tree
#'
#' @param plot a plot produced by ggtree to add the piecharts to
#' @inheritParams prepare_tree
#' @param signature the signature for which to plot the piecharts
#' @param sig_matrix the signature
#' @param cex the size of the piechart
#'
#' @return the ggplot with the piecharts added
add_contri_pie_to_tree <- function(plot, tree, signature, cex) {
  fraction = pull(tree@data, signature)

  # make the plotting df
  perc_df = tibble(
    node = tree@data$node,
    fraction = fraction
  )
  perc_df = perc_df[!is.na(perc_df$fraction), ]
  perc_df$other = 1 - perc_df$fraction

  # Plot labels
  # this part of the code has been adapted from ggtree::nodepie
  # all credits go to Guangchuang Yu: https://bioconductor.org/packages/devel/bioc/vignettes/ggtree/inst/doc/ggtree.html
  type <- value <- NULL
  cols = 2:3
  ldf <- gather(perc_df, type, value, !!cols) %>% split(., .$node)
  pies = lapply(ldf, function(df) {
    ggplot(df, aes(x = 1, y = value)) +
      geom_segment(y = 0, yend = 1, x = 1.45, xend = 1.45, size = 1) +
      geom_bar(stat = "identity", aes(fill = type)) +
      scale_fill_manual(values = c('fraction' = 'black', 'other' = 'white')) +
      scale_y_continuous(limits = c(0,1)) +
      coord_polar(theta = "y", clip = 'off') +
      theme_inset()
  })
  ggtree::inset(plot, pies, width = 0.08*cex, height = 0.08*cex, x = 'branch')
}

# color ------------------------------------------------------------
dist_cols <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4',
               '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff',
               '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
               '#000075', '#808080', '#b8ccd1', '#000000')
dist_cols31 <- c('#ff0000', '#7f2200', '#ffe1bf', '#474d00', '#ccff00', '#00593c', '#00ccff', '#002e73', '#220033',
                 '#ff00aa', '#ff4059', '#330000', '#ffa280', '#8c5e00', '#fbffbf', '#71b359', '#00ffcc', '#005580',
                 '#0044ff', '#cc00ff', '#7f0044', '#806460', '#ff8800', '#ffcc00', '#64664d', '#00ff22', '#7ca69d',
                 '#0088ff', '#d0bfff', '#9553a6', '#ffbfd9')
dist_cols50 <- c('#808080', '#c0c0c0', '#2f4f4f', '#556b2f', '#8b4513', '#2e8b57', '#7f0000', '#191970',
                 '#006400', '#808000', '#5f9ea0', '#b8860b', '#4682b4', '#d2691e', '#9acd32', '#cd5c5c',
                 '#00008b', '#32cd32', '#8fbc8f', '#8b008b', '#d2b48c', '#ff4500', '#ffa500', '#ffd700',
                 '#0000cd', '#00ff00', '#9400d3', '#00ff7f', '#4169e1', '#dc143c', '#00ffff', '#00bfff',
                 '#adff2f', '#ff6347', '#ff00ff', '#db7093', '#f0e68c', '#ffff54', '#6495ed', '#dda0dd',
                 '#87ceeb', '#ff1493', '#ffa07a', '#ee82ee', '#98fb98', '#7fffd4', '#ff69b4', '#fffacd', '#e0ffff', '#ffb6c1')
grad_cols <- c(gplots::colorpanel(n = 50, low = '#2f2bad', mid = '#93cfe2', high = '#fffaaa'),
               gplots::colorpanel(n = 50, low = '#fffaaa', mid = "#ffad66", high = '#d60000'))

# plot_96_tree ----------------------------------------
#' an adjustment on MutationalPatterns' plot_96_profile
#' @return a ggplot
#' @export
plot_96_profile2 <- function(mut_matrix,
                         colors = NA,
                         ymax = NULL,
                         condensed = FALSE,
                         addLine = TRUE,
                         small = FALSE,
                         axisnormal = FALSE,
                         normalize = TRUE,
                         ann_size = 3,
                         label_x = "A[T>C]A",
                         ncol = 1) {
  mut_matrix = check_0_mm(mut_matrix)

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  freq <- full_context <- substitution <- context <- NULL
  #browser()
  # Check color vector length
  # Colors for plotting
  if (MutationalPatterns:::.is_na(colors)) {
    colors <- MutationalPatterns:::COLORS6
  }
  if (length(colors) != 6) {
    stop("Provide colors vector with length 6", call. = FALSE)
  }

  # Make contribution relative
  if (normalize) {
    norm_mut_matrix <- apply(mut_matrix, 2, function(x) x / sum(x))
    if (is.null(ymax)) ymax = 0.2
    yaxis_step = 0.1
  } else {
    norm_mut_matrix = mut_matrix
    if (is.null(ymax)) ymax = 0.2
    yaxis_step = 5
  }

  # Get substitution and context from rownames and make long.
  #browser()
  tb <- norm_mut_matrix %>%
    as.data.frame() %>%
    tibble::rownames_to_column("full_context") %>%
    dplyr::mutate(
      substitution = stringr::str_replace(full_context, "\\w\\[(.*)\\]\\w", "\\1"),
      context = stringr::str_replace(full_context, "\\[.*\\]", "\\.")
    ) %>%
    dplyr::select(-full_context) %>%
    tidyr::pivot_longer(c(-substitution, -context), names_to = "sample", values_to = "freq") %>%
    dplyr::mutate(sample = factor(sample, levels = unique(sample)),
                  plot_x = factor(rep(rownames(norm_mut_matrix), each = ncol(mut_matrix)), levels = rownames(norm_mut_matrix)))

  ucon = unique(tb$context)
  tb_text = tb %>% filter(plot_x == label_x)

  # Change plotting parameters based on whether plot should be condensed.
  if (condensed == TRUE) {
    width <- 1
    spacing <- 0
  } else {
    width <- 0.6
    spacing <- 0.1
  }
  # Create figure
  plot = ggplot(mapping = aes(
    fill = substitution,
    color = substitution,
    width = width
  )) + scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    geom_bar(data = tb, aes(x = plot_x, y = freq), stat = "identity", size = .2) +
    geom_text(data = tb_text, mapping = aes(label = sample, x = plot_x, y = ymax*0.7), color = 'black', size = ann_size)

  plot = plot +
    facet_wrap(~sample, ncol = ncol, dir = 'v') +
    labs(x = '', y = "relative contribution") +
    coord_cartesian(ylim = c(0, (ymax+0.2*ymax))) +
    scale_y_continuous(breaks = seq(0, ymax, yaxis_step), expand = c(0,0)) +
    scale_x_discrete(labels = rep(ucon, 6)) +
    guides(fill = FALSE, color = FALSE) +
    theme(
      panel.spacing.x = unit(spacing, "lines"),
      strip.background = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.line	= element_line(colour = "black", linewidth = 1),
      axis.title.y = element_text(vjust = 1),
      axis.ticks = element_line(colour = "black"),
      strip.background.y = element_blank(),
      strip.text = element_blank(),
      panel.grid.major.x = element_blank()
    )
  if (!axisnormal) plot = plot + theme(axis.line = element_line(linewidth = 0.5))
  plot = plot + theme(text = element_text(color = 'black'),
                      axis.text.y = element_text(color = 'black'),
                      axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5, color = 'black'))
  return(plot)
}
