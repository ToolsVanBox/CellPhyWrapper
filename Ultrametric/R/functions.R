###' functions to make ultrametric trees


#' callable loci correction

#' Function to calculate the product of (1 - sensitivity) for each sample in a node
calculate_product_sensitivity <- function(samples, df_sensitivity) {
  #' Split the 'samples' string by '|'
  sample_list <- strsplit(samples, "\\|")[[1]]

  #' Retrieve sensitivities for these samples
  sensitivities <- df_sensitivity$sensitivity[df_sensitivity$SAMPLE %in% sample_list]

  #' Calculate the product of (1 - sensitivity)
  if (length(sensitivities) == 0) return(1)  # fallback to avoid division by zero
  return(1 - prod(1 - sensitivities, na.rm = TRUE))
}

#' function to take a tree and correct the branch lengths based on the product sensitivity
#' requires samples (samples column from output of cellphywrapper object: tree@data$samples .
#' requires a CALLABLE df (containing columns SAMPLES and CALLABLE (in nr of bases)
#' Function to calculate the product of (1 - sensitivity) for each sample in a node
correct_branches <- function(tree, callable_df){
  #' get sensitivity df by taking fraction of total possible callable loci
  max_callable <- 2745186691
  callable_df$sensitivity <- callable_df$CALLABLE / max_callable
  sensitivity_df <- callable_df[c('SAMPLE','sensitivity')]
  print(sensitivity_df)
  #' calculate the sensitivity per node/branch using the calculate_product_sensitivity function
  tree@data$product_sensitivity <- sapply(tree@data$samples, calculate_product_sensitivity, sensitivity_df)

  #' correct the branch lengths
  tree@data$corr_branch_lengths <- tree@data$branch_length / tree@data$product_sensitivity

  tree@data$corr_branch_lengths[is.na(tree@data$corr_branch_lengths)] <- 0 #' for merged germline sample

  #' store also in the phylo object, required for downstream analyses
  tree@phylo$edge.length <- tree@data[order(tree@data$node),]$corr_branch_lengths

  return(tree)
}


####' ultrametric tree making step 1: remove prenatal region of tree

#' Function to recursively subtract mutations as we traverse the tree from the root >> remove prenatal section of the tree
#' use find root node as initial input for 'node' and the mut load at birth as initial 'remaining_mutations'
subtract_mutations <- function(tree, node, remaining_mutations) {
  #' Get the edges that descend from the current node
  child_edges <- which(tree$edge[, 1] == node)

  for (edge_index in child_edges) {

    #' Get the length of edge
    edge_length <- tree$edge.length[edge_index]

    if (remaining_mutations > 0) { #' in this part store a 'new' remaining_mutations2 so that the loop retains the original value!
      if (edge_length <= remaining_mutations) {
        #' If the edge length is shorter than or equal to the 'remaining mutations', set the whole edge to length 0
        remaining_mutations2 <- remaining_mutations - edge_length
        tree$edge.length[edge_index] <- 0
      } else {
        #' If the edge length is longer than the remaining mutations, subtract the 'remaining mutations'
        tree$edge.length[edge_index] <- edge_length - remaining_mutations
        remaining_mutations2 <- 0
      }

      #' Recursively process the child node (if there are remaining mutations)
      child_node <- tree$edge[edge_index, 2]

      for (single_child in child_node){
        tree <- subtract_mutations(tree, single_child, remaining_mutations2)
      }

    }
  }

  return(tree)
}


#' ultrametric making ##' FROM FABRE ET AL 2022 NATURE https://doi.org/10.1038/s41586-022-04785-z
make.ultrametric.tree <- function(tree) {
  root.number <- length(tree$tip.label) + 1
  ultra.tree <- length.normalise(tree, tree, root.number, 1)

  ultra.tree$edge.length[ultra.tree$edge.length == Inf] <- 0 ##' the outgroup gives an Inf value after ultrametric making, fix like this

  return(ultra.tree)
}

length.normalise <- function(orig.tree, new.tree, curr.node, remaining.stick) {
  curr.node.children <- unlist(Descendants(orig.tree, curr.node, 'children')) #' 'children' mode?

  for (j in curr.node.children) {
    index <- which(orig.tree$edge[,2] == j & orig.tree$edge[,1] == curr.node, arr.ind = TRUE)

    if (j %in% orig.tree$tip.label) {
      new.tree$edge.length[index] <- remaining.stick
    } else {
      curr.node.tips <- unlist(Descendants(orig.tree, j, 'tips')) #' 'tips' mode?
      curr.dist <- find.distance(orig.tree, curr.node, j)

      if (curr.dist == 0) {curr.dist <- 0.01} #' So that no edge lengths are zero
      desc.lengths <- sapply(curr.node.tips, FUN = find.distance, tree = orig.tree, from = curr.node)
      new.tree$edge.length[index] <- remaining.stick * curr.dist / mean(desc.lengths)
      shorter.stick <- remaining.stick - new.tree$edge.length[index]

      #' Call as recursive function
      new.tree <- length.normalise(orig.tree, new.tree, j, shorter.stick)
    }
  }
  return(new.tree)
}

find.distance <- function(tree, from, to) {
  path <- nodepath(tree, from, to)
  res <- 0
  for (j in 2:length(path)) {
    index <- which(tree$edge[,2] == path[j] & tree$edge[,1] == path[j-1], arr.ind = TRUE)
    res <- res + tree$edge.length[index]
  }
  return(res)
}

####' from here on Lucca functions again
####' age scaling of ultrametric trees (ultrametric is initially scaled from 0-1)
##' perform age scaling and store branch lengths in the data layer again
perform_age_scaling <- function(tree_obj, age_in_years){
  tree_obj@data$branch_length <- round(tree_obj@phylo$edge.length[names(tree_obj@data$branch_length)] * age_in_years, 1)
  return(tree_obj)
}


###' ultrametric making workflow for postnatal part of the tree (takes output of correct_branches as input)

make_postnatal_ultra <- function(tree_input, mutload_birth, age_at_sampling){
  #' Get the root node (in a rooted tree, root is node Ntip + 1)
  root_node <- Ntip(tree_input@phylo) + 1

  ##' use subtract_mutations function to remove prenatal section of tree. Start from the root and remove 79 mutations
  tree_output <- tree_input #' store copy of tree to work with
  tree_output@phylo <- subtract_mutations(tree_output@phylo, root_node, mutload_birth)

  #' perform scaling of post-natal tree to make ultrametric
  tree_output@phylo <- make.ultrametric.tree(tree_output@phylo)

  #' make age corrected
  tree_output <- perform_age_scaling(tree_output, age_at_sampling)

  return(tree_output)
}


####' ultrametric tree making step 2: isolate only prenatal region of tree
#' Function to recursively traverse the tree and cut the max distance traversed to 'remaining mutations' (start at mutation load at birth)
traverse_and_cut <- function(tree, node, remaining_mutations = 79) {
  #' Get the edges that descend from the current node
  child_edges <- which(tree$edge[, 1] == node)

  for (edge_index in child_edges) {

    #' Get the length of edge
    edge_length <- tree$edge.length[edge_index]

    if (remaining_mutations > 0) { #' in this part store a 'new' remaining_mutations2 so that the loop retains the original value!
      if (edge_length <= remaining_mutations) {
        #' If the edge length is shorter than or equal to the 'remaining mutations', keep whole edge length
        #' update remaining mutations to subtract the edge length of this edge
        remaining_mutations2 <- remaining_mutations - edge_length

        #' Recursively process the child node (if there are remaining mutations)
        child_node <- tree$edge[edge_index, 2]

        for (single_child in child_node){
          tree <- traverse_and_cut(tree, single_child, remaining_mutations2)
        }

      } else {

        #' If the edge length is longer than the remaining mutations, set to 'remaining mutations'
        tree$edge.length[edge_index] <- remaining_mutations
        remaining_mutations2 <- 0

        #' set length of descending branches to 0 as these are all cut off
        ##' first get the index of all the nodes and tips downstream of current node
        if (tree$edge[edge_index, 2] > Ntip(tree)) { #this excludes tips as these have no downstream data

          #' get descending nodes
          descendants <- Descendants(tree, tree$edge[edge_index, 2], type = 'all')

          #' set length of the edges belonging to those nodes to 0
          tree$edge.length[which(tree$edge[, 2] %in% descendants)] <- 0

        }

      }

    }

  }

  return(tree)
}

###' function to create ultrametric prenatal tree
make_prenatal_ultra <- function(raw_tree, mutations_at_birth){

  #' first retrieve prenatal tree
  raw_tree@phylo <- traverse_and_cut(raw_tree@phylo, Ntip(raw_tree) + 1, mutations_at_birth)

  #' make prenatal tree ultrametric
  raw_tree@phylo <- make.ultrametric.tree(raw_tree@phylo)
  raw_tree@data$branch_length <- round(raw_tree@phylo$edge.length[names(raw_tree@data$branch_length)] * 0.75, 2)

  return(raw_tree)
}

#######' function to merge the postnatal and prenatal tree
merge_prenatal_postnatal <- function(prenatal_tree, postnatal_tree){

  #' retrieve lengths to be added to postnatal tree
  lengths_to_add <- prenatal_tree@data$branch_length[prenatal_tree@data$branch_length > 0]
  lengths_to_add <- lengths_to_add[!is.na(lengths_to_add)]

  #' add lengths to postnatal tree
  for (nodename in names(lengths_to_add)){
    postnatal_tree@data[postnatal_tree@data$node_id == nodename, 'branch_length'] <- postnatal_tree@data[postnatal_tree@data$node_id == nodename, 'branch_length'] + lengths_to_add[nodename]

  }

  #' also add this data to the phylo object
  postnatal_tree@phylo$edge.length <- postnatal_tree@data[order(postnatal_tree@data$node),]$branch_length

  return(postnatal_tree)
}


####' full wrapper for ultrametric tree building from 'correct_branches' onwards.
####' uses output of correct_branches

create_wholelife_tree <- function(corr_tree, age_HSC_at_FU, mutload_birth = 79){

  postnatal_tree <- make_postnatal_ultra(corr_tree, mutload_birth, age_HSC_at_FU)
  prenatal_tree <- prenatal_final_tree <- make_prenatal_ultra(corr_tree, mutload_birth)

  wholelife_tree <- merge_prenatal_postnatal(prenatal_tree, postnatal_tree)

  return(wholelife_tree)

}


##################' retrieving the distances from root to branch end
##' gives output for upper, lower, and estimated nr of mutations at birth
##' input: output of correct_branches

###' plotting interval of branch timing (based on different possible mutation loads at birth)
get_interval_branch <- function(CLcorrTree, branch_of_interest, age_HSC_at_FU,
                                ave_mutload_birth = 79,
                                lower_mutload_birth = 50,
                                higher_mutload_birth=108.9){

  #' get trees with upper and lower boundary of HSC load at birth
  final_tree <- create_wholelife_tree(CLcorrTree, age_HSC_at_FU, ave_mutload_birth)
  final_treeHigh <- create_wholelife_tree(CLcorrTree, age_HSC_at_FU, lower_mutload_birth)
  final_treeLow <- create_wholelife_tree(CLcorrTree, age_HSC_at_FU, higher_mutload_birth)

  #' get node
  final_tree = prepare_tree(final_tree)
  node_of_interest = final_tree@data[final_tree@data$branch_id == branch_of_interest, ][['node']]

  #' get values
  middleTime <- find.distance(final_tree@phylo, Ntip(final_tree@phylo) + 1, node_of_interest)
  lowerTime <- find.distance(final_treeLow@phylo, Ntip(final_treeLow@phylo) + 1, node_of_interest)
  higherTime <- find.distance(final_treeHigh@phylo, Ntip(final_treeHigh@phylo) + 1, node_of_interest)
  return(c(middleTime, lowerTime, higherTime))
}
