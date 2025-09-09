# check_tree_size()
# Check that a tree is an expected size before calculating simulation results
# tree_list: A list of trees to check
# expected_tree_size: The expected size of each tree
# Returns TRUE if all trees have tip size == expected_tree_size, and false otherwise.
check_tree_size = function(tree_list, expected_tree_size){
  expected = TRUE
  
  #Iterator through tree list.
  for(i in 1:length(tree_list)){
    if(length(tree_list[[i]]$tip.label) != expected_tree_size){
      expected = FALSE
    }
  }
  return(expected)
}

# calculate_mast_random_pairs()
# Calculate MAST between random pairings of tree lists. 
# trees1 and trees2 are lists of trees with identical tips and same number of trees in each list
# cores: numbers of cores to use.
# Returns a list of MAST distances between random pairs of trees in the lists.
calculate_mast_random_pairs = function(trees1, trees2, cores){
  
  print("Registering cores...")
  num_cores <- cores
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  print("Registering cores...DONE")
  
  # Generate random pairings between trees of each list.
  pair_indices <- data.frame(
    trees1_idx = sample(seq_along(trees1), replace = FALSE), 
    trees2_idx = sample(seq_along(trees2), replace = FALSE)
  )
  
  # Parallel computation of mast distance
  mast_dist <- foreach(i = 1:nrow(pair_indices), .combine = c, .packages = "phangorn", .export = c("find_clades_size_n", "get_clade_tips", "matching_clades")) %dopar% {
    length(mast(x = phangorn::midpoint(trees1[[pair_indices[i, 1]]]), y = phangorn::midpoint(trees2[[pair_indices[i, 2]]]), tree = FALSE, rooted = TRUE))
  }
  
  # Stop the cluster
  stopCluster(cl)
  return(mast_dist)
}

# calculate_shared_subtrees_dist()
# Returns number of matching clades between trees that are size n - all possible pairings
# trees1 and trees2: Lists of two different trees. Trees should have same tips, and lists should be equal lengths
# cores: number of cores to use for computation
# n: Clade size (# of tips) that should be checked for equality across trees.
# Returns a list of percent of clades of size n that were present in both trees for each pairing of trees. 
calculate_shared_subtrees_dist = function(trees1, trees2, cores, n){
  
  trees1 = phangorn::midpoint(trees1)
  trees2 = phangorn::midpoint(trees2)
  
  print("Registering cores...")
  num_cores <- cores
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  print("Registering cores...DONE")
  
  # Generate all unique tree pairs
  pair_indices <- expand.grid(seq_along(trees1), seq_along(trees2))
  
  # Parallel computation of clade similarity scores.
  clade_dist <- foreach(i = 1:nrow(pair_indices), .combine = c, .packages = "phangorn", .export = c("find_clades_size_n", "get_clade_tips", "matching_clades")) %dopar% {
    matching_clades(trees1[[pair_indices[i, 1]]], trees2[[pair_indices[i, 2]]], n = n)
  }
  
  # Stop the cluster
  stopCluster(cl)
  return(clade_dist)
  
}

# calculate_shared_subtrees_dist_random_pairs()
# Returns number of matching clades between trees that are size n - random pairings between lists.
# trees1 and trees2: Lists of two different trees. Trees should have same tips, and lists should be equal lengths
# cores: number of cores to use for computation
# n: Clade size (# of tips) that should be checked for equality across trees.
# Returns a list of percent of clades of size n that were present in both trees for eachh pairing of trees. 
calculate_shared_subtrees_dist_random_pairs = function(trees1, trees2, cores, n){
  
  print("Registering cores...")
  num_cores <- cores
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  print("Registering cores...DONE")
  
  # Generate all unique tree pairs
  pair_indices <- data.frame(
    trees1_idx = sample(seq_along(trees1), replace = FALSE), 
    trees2_idx = sample(seq_along(trees2), replace = FALSE)
  )
  
  # Parallel computation of matching_clades
  clade_dist <- foreach(i = 1:nrow(pair_indices), .combine = c, .packages = "phangorn", .export = c("find_clades_size_n", "get_clade_tips", "matching_clades")) %dopar% {
    matching_clades(phangorn::midpoint(trees1[[pair_indices[i, 1]]]), phangorn::midpoint(trees2[[pair_indices[i, 2]]]), n = n)
  }
  
  # Stop the cluster
  stopCluster(cl)
  return(clade_dist)
  
}

# # calculate_shared_subtrees_dist_random_pairs_no_midpoint()
# Returns number of matching clades between trees that are size n - random pairings between lists.
# Trees are not midpoint rooted first, making this computationally intractable for large trees.
# trees1 and trees2: Lists of two different trees. Trees should have same tips, and lists should be equal lengths
# cores: number of cores to use for computation
# n: Clade size (# of tips) that should be checked for equality across trees.
# Returns a list of percent of clades of size n that were present in both trees for pairing of trees. 
calculate_shared_subtrees_dist_random_pairs_no_midpoint = function(trees1, trees2, cores, n){
  
  print("Registering cores...")
  num_cores <- cores
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  print("Registering cores...DONE")
  
  # Generate random pairings of trees between lists. 
  pair_indices <- data.frame(
    trees1_idx = sample(seq_along(trees1), replace = FALSE), 
    trees2_idx = sample(seq_along(trees2), replace = FALSE)
  )
  
  # Parallel computation matching clades. 
  clade_dist <- foreach(i = 1:nrow(pair_indices), .combine = c, .packages = "phangorn", .export = c("find_clades_size_n", "get_clade_tips", "matching_clades")) %dopar% {
    matching_clades(trees1[[pair_indices[i, 1]]], trees2[[pair_indices[i, 2]]], n = n)
  }
  
  # Stop the cluster
  stopCluster(cl)
  return(clade_dist)
  
}

# find_clades_size_n()
# Function to find clades of size n
# Returns the number of clades in the tree which are size n
# tree: A phylogenetic tree
# n: the clade size to look for
find_clades_size_n <- function(tree, n) {
  node_sizes <- sapply(Ntip(tree) + 1:tree$Nnode, function(node) {
    length(Descendants(tree, node, type = "tips")[[1]])
  })
  
  clade_nodes <- which(node_sizes == n) + Ntip(tree)  # Convert to node indices
  return(clade_nodes)
}

# get_clade_tips()
# Function to extract tip labels from a given node
# tree: A phylogenetic tree
# node: the node to extract tips from. 
# Returns a list of tips which are descendants from that node.
get_clade_tips <- function(tree, node) {
  tips <- Descendants(tree, node, type = "tips")[[1]]
  return(sort(tree$tip.label[tips]))  # Sort for consistent comparison
}

# matching_clades()
# Returns the fraction of clades of size n which are shared between tree1 and tree2.
# tree1 and tree2: Two phylogenetic trees. Tips should be identical between the two trees. 
# n: Clade size to look for matching clades.
matching_clades = function(tree1, tree2, n){
  
  # Get clade nodes of size n for both trees
  print("Getting clade sizes...")
  clade_nodes_tree1 <- find_clades_size_n(tree1, n)
  clade_nodes_tree2 <- find_clades_size_n(tree2, n)
  print("Getting clade sizes...DONE")
  
  # Extract clade tips from Tree 1
  print("Extracting tips from relevant nodes...")
  clades_tree1 <- lapply(clade_nodes_tree1, function(node) get_clade_tips(tree1, node))
  
  # Extract clade tips from Tree 2
  clades_tree2 <- lapply(clade_nodes_tree2, function(node) get_clade_tips(tree2, node))
  
  print("Extracting tips from relevant nodes...DONE")
  
  # Compare clades between trees
  print("Determining Matching clades...")
  matching_clades <- sapply(clades_tree1, function(clade1) {
    any(sapply(clades_tree2, function(clade2) identical(clade1, clade2)))
  })
  print("Determining Matching clades...DONE")
  
  # Return matching score.
  n_matching = sum(matching_clades)
  return((2*n_matching / (length(clade_nodes_tree1) + length(clade_nodes_tree2))))
  
}

# n_matching_clades()
# Returns the number of matching clades between two trees of size n.
# tree1 and tree2: Two phylogenetic trees. Tips should be identical between the two trees. 
# n: Clade size to look for matching clades.
n_matching_clades = function(tree1, tree2, n){
  
  # Get clade nodes of size n
  print("Getting clade sizes...")
  clade_nodes_tree1 <- find_clades_size_n(tree1, n)
  clade_nodes_tree2 <- find_clades_size_n(tree2, n)
  print("Getting clade sizes...DONE")
  
  # Extract clade tips from Tree 1
  print("Extracting tips from relevant nodes...")
  clades_tree1 <- lapply(clade_nodes_tree1, function(node) get_clade_tips(tree1, node))
  
  # Extract clade tips from Tree 2
  clades_tree2 <- lapply(clade_nodes_tree2, function(node) get_clade_tips(tree2, node))
  
  print("Extracting tips from relevant nodes...DONE")
  
  # Compare clades between trees
  print("Determining Matching clades...")
  matching_clades <- sapply(clades_tree1, function(clade1) {
    any(sapply(clades_tree2, function(clade2) identical(clade1, clade2)))
  })
  print("Determining Matching clades...DONE")
  
  # Summ matching clades.
  n_matching = sum(matching_clades)
  return(n_matching)
  
}

# nanostructure_dist_percent_tips()
# Calculates the percent of tips which in clades of size n which have matching membership between two trees 
# and then converts to a distance metric
# tree1 and tree2: Two phylogenetic trees with identical tips. 
# cores: number of cores to use for the computation. 
# tree size: number of tips that are in both trees. 
nanostructure_dist_percent_tips = function(trees1, trees2, cores, n, tree_size){
  
  print("Registering cores...")
  num_cores <- cores
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  print("Registering cores...DONE")
  
  # Generate all unique tree pairs
  pair_indices <- expand.grid(seq_along(trees1), seq_along(trees2))
  
  # Parallel computation 
  clade_dist <- foreach(i = 1:nrow(pair_indices), .combine = c, .packages = "phangorn", .export = c("find_clades_size_n", "get_clade_tips", "n_matching_clades")) %dopar% {
    1 - n*(n_matching_clades(trees1[[pair_indices[i, 1]]], trees2[[pair_indices[i, 2]]], n = n))/tree_size
    
  }
  
  # Stop the cluster
  stopCluster(cl)
  return(clade_dist)
}

# calaculate_msd_distance_pairs()
# Check msd distance of all possible pairwise comparisons.
# trees1: a list of trees. 
# trees2: a second list of tress.
# Cores: the number of cores to use.
# RETURN: a list of MSD distances.
calculate_msd_distance = function(trees1, trees2, cores){
  print("Registering cores...")
  
  # Set up parallel backend
  num_cores <- cores
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  print("Registering cores...DONE")
  
  # Generate all unique tree pairs
  pair_indices <- expand.grid(seq_along(trees1), seq_along(trees2))
  
  # Parallel computation of MSD distances
  msd_distance <- foreach(i = 1:nrow(pair_indices), .combine = c, .packages = "phangorn") %dopar% {
    BigTreeDist::MatchingSplitDistance(trees1[[pair_indices[i, 1]]], trees2[[pair_indices[i, 2]]])
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  return(msd_distance)
  
}

# calaculate_msd_distance_random_pairs()
# Check msd distance between random pairings of trees. 
# trees1: a list of trees. 
# trees2: a second list of tress.
# Cores: the number of cores to use.
# RETURN: a list of MSD distances.
calculate_msd_distance_random_pairs = function(trees1, trees2, cores){
  print("Registering cores...")
  
  # Set up parallel backend
  num_cores <- cores
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  print("Registering cores...DONE")
  
  # Generate all unique tree pairs
  pair_indices <- data.frame(
    trees1_idx = sample(seq_along(trees1), replace = FALSE), 
    trees2_idx = sample(seq_along(trees2), replace = FALSE)
  )
  
  # Parallel computation of MSD distances
  msd_distance <- foreach(i = 1:nrow(pair_indices), .combine = c, .packages = "phangorn") %dopar% {
    BigTreeDist::MatchingSplitDistance(trees1[[pair_indices[i, 1]]], trees2[[pair_indices[i, 2]]])
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  return(msd_distance)
  
}

# calaculate_mcd_distance()
# Mutual clustering distance (MCD)
# Check mcd distance of all possible pairwise comparisons
# trees1: a list of trees. 
# trees2: a second list of tress.
# Cores: the number of cores to use.
# RETURN: a list of MCD distances.
calculate_mcd_distance = function(trees1, trees2, cores){
  print("Registering cores...")
  
  # Set up parallel backend
  num_cores <- cores
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  print("Registering cores...DONE")
  
  # Generate all unique tree pairs
  pair_indices <- expand.grid(seq_along(trees1), seq_along(trees2))
  
  # Parallel computation of MCD distances
  mcd_distance <- foreach(i = 1:nrow(pair_indices), .combine = c, .packages = "phangorn") %dopar% {
    BigTreeDist::ClusteringInfoDistance(trees1[[pair_indices[i, 1]]], trees2[[pair_indices[i, 2]]], normalize = TRUE)
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  return(mcd_distance)
  
}

# calaculate_mcd_distance_pairs()
# Mutual clustering distance (MCD)
# Check mcd distance of random pairs of trees in trees1 and trees2
# trees1: a list of trees. 
# trees2: a second list of tress.
# Cores: the number of cores to use
# RETURN: a list of MCD distances.
calculate_mcd_distance_random_pairs = function(trees1, trees2, cores){
  print("Registering cores...")
  
  # Set up parallel backend
  num_cores <- cores
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  print("Registering cores...DONE")
  
  # Generate all unique tree pairs
  pair_indices <- data.frame(
    trees1_idx = sample(seq_along(trees1), replace = FALSE), 
    trees2_idx = sample(seq_along(trees2), replace = FALSE)
  )
  # Parallel computation of MCD distances
  mcd_distance <- foreach(i = 1:nrow(pair_indices), .combine = c, .packages = "phangorn") %dopar% {
    BigTreeDist::ClusteringInfoDistance(trees1[[pair_indices[i, 1]]], trees2[[pair_indices[i, 2]]], normalize = TRUE)
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  return(mcd_distance)
  
}


# calaculate_nye_distance_pairs()
# Return Nye distance of all possible pairwise comparisons
# trees1: a list of trees. 
# trees2: a second list of tress.
# Cores: the number of cores to use
calculate_nye_distance = function(trees1, trees2, cores){
  
  
  print("Registering cores...")
  
  # Set up parallel backend
  num_cores <- cores
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  print("Registering cores...DONE")
  
  # Generate all unique tree pairs
  pair_indices <- expand.grid(seq_along(trees1), seq_along(trees2))
  
  # Parallel computation of RF distances
  nye_distance <- foreach(i = 1:nrow(pair_indices), .combine = c, .packages = "phangorn") %dopar% {
    BigTreeDist::NyeSimilarity(trees1[[pair_indices[i, 1]]], trees2[[pair_indices[i, 2]]], similarity = FALSE, normalize = TRUE)
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  return(nye_distance)
  
}

# calaculate_nye_distance_pairs()
# Return Nye distance random pairings
# trees1: a list of trees. 
# trees2: a second list of tress.
# Cores: the number of cores to use
calculate_nye_distance_random_pairs = function(trees1, trees2, cores){
  
  print("Registering cores...")
  
  # Set up parallel backend
  num_cores <- cores
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  print("Registering cores...DONE")
  
  # Generate all unique tree pairs
  pair_indices <- data.frame(
    trees1_idx = sample(seq_along(trees1), replace = FALSE), 
    trees2_idx = sample(seq_along(trees2), replace = FALSE)
  )
  
  # Parallel computation of Nye distances
  nye_distance <- foreach(i = 1:nrow(pair_indices), .combine = c, .packages = c("phangorn", "BigTreeDist")) %dopar% {
    BigTreeDist::NyeSimilarity(trees1[[pair_indices[i, 1]]], trees2[[pair_indices[i, 2]]], similarity = FALSE, normalize = TRUE)
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  return(nye_distance)
  
}

# calculate_jrf_dist()
# Return JRF  distance between all possible pairings of trees.
# trees1 and trees2: a list of trees. 
# Cores: the number of cores to utilize
calculate_jrf_distance = function(trees1, trees2, cores){
  
  print("Registering cores...")
  
  # Set up parallel backend
  num_cores <- cores
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  print("Registering cores...DONE")
  
  # Generate all unique tree pairs
  pair_indices <- expand.grid(seq_along(trees1), seq_along(trees2))
  
  # Parallel computation of RF distances
  jrf_distance <- foreach(i = 1:nrow(pair_indices), .combine = c, .packages = "BigTreeDist") %dopar% {
    BigTreeDist::JaccardRobinsonFoulds(trees1[[pair_indices[i, 1]]], trees2[[pair_indices[i, 2]]])
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  return(jrf_distance)
  
}

# calculate_jrf_dist()
# Return list of JRF distances between random pairs of trees in trees1 and trees2 lists.
# trees1 and trees2: a list of trees. Trees should have identical tips
# Cores: the number of cores to utilize
calculate_jrf_distance_random_pairs = function(trees1, trees2, cores){
  
  print("Registering cores...")
  
  # Set up parallel backend
  num_cores <- cores
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  # 
  print("Registering cores...DONE")
  
  # Generate all unique tree pairs
  pair_indices <- data.frame(
    trees1_idx = sample(seq_along(trees1), replace = FALSE), 
    trees2_idx = sample(seq_along(trees2), replace = FALSE)
  )
  # Parallel computation of JRF distances
  jrf_distance <- foreach(i = 1:nrow(pair_indices), .combine = c, .packages = "BigTreeDist") %dopar% {
    BigTreeDist::JaccardRobinsonFoulds(trees1[[pair_indices[i, 1]]], trees2[[pair_indices[i, 2]]])
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  return(jrf_distance)
  
}


# calculate_path_dist()
# Check path distance between every possible pairing of trees in list 1 and 2.
# trees1 and trees2: a list of trees. 
# Cores: the number of cores to utilize
# RETURN: a list of path distances distances.
calculate_path_dist = function(trees1, trees2, cores = 20){
  
  if(length(trees1$tip.labels) != length(trees2$tip.labels)){
    stop("Trees are not of the same size.")
  }
  print("Registering cores...")
  
  # Set up parallel backend
  num_cores <- cores
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  print("Registering cores...DONE")
  
  # Generate all unique tree pairs
  pair_indices <- expand.grid(seq_along(trees1), seq_along(trees2))
  
  # Parallel computation of Path distances
  path_distance <- foreach(i = 1:nrow(pair_indices), .combine = c, .packages = "TreeDist") %dopar% {
    PathDist(trees1[[pair_indices[i, 1]]], trees2[[pair_indices[i, 2]]])
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  return(path_distance)
  
}

# calculate_path_dist_random_pairs()
# Returns path distance for random pairings of trees in list1 and list2
# trees1 and trees2: a list of trees. 
# Cores: the number of cores to utilize
calculate_path_dist_random_pairs = function(trees1, trees2, cores = 20){
  
  if(length(trees1$tip.labels) != length(trees2$tip.labels)){
    stop("Trees are not of the same size.")
  }
  print("Registering cores...")
  
  # Set up parallel backend
  num_cores <- cores
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  print("Registering cores...DONE")
  
  # Generate random pairings between tree lists.
  pair_indices <- data.frame(
    trees1_idx = sample(seq_along(trees1), replace = FALSE), 
    trees2_idx = sample(seq_along(trees2), replace = FALSE)
  )
  
  # Parallel computation of path distances.
  path_distance <- foreach(i = 1:nrow(pair_indices), .combine = c, .packages = "BigTreeDist") %dopar% {
    BigTreeDist::PathDist(trees1[[pair_indices[i, 1]]], trees2[[pair_indices[i, 2]]])
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  return(path_distance)
  
}

# calculate_path_dist()
# Returns a list of path distance between all possible pairings of tree in list one to list two
# trees1 and trees2: a list of trees. 
# Cores: the number of cores to utilize
calculate_path_dist = function(trees1, trees2, cores = 20){
  
  if(length(trees1$tip.labels) != length(trees2$tip.labels)){
    stop("Trees are not of the same size.")
  }
  print("Registering cores...")
  
  # Set up parallel backend
  num_cores <- cores
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  print("Registering cores...DONE")
  
  # Generate all unique tree pairs
  pair_indices <- expand.grid(seq_along(trees1), seq_along(trees2))
  
  # Parallel computation of path distances
  path_distance <- foreach(i = 1:nrow(pair_indices), .combine = c, .packages = "TreeDist") %dopar% {
    PathDist(trees1[[pair_indices[i, 1]]], trees2[[pair_indices[i, 2]]])
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  return(path_distance)
  
}

# write_ref_dropped_trees()
# Returns a set of filepaths to newly written trees where guides have been removed.
# file_list: a list of file paths to trees to read
# pattern: the pattern to use to remove the guide sequences.
# Writes out trees with guide reference tips pruned.
write_ref_dropped_trees = function(file_list, pattern = "^REF."){
  
  trees.path = vector()
  
  for(i in 1:length(file_list)){
    cur_file = file_list[i]
    cur_tree = read.tree(cur_file)
    
    ref_tips = cur_tree$tip.label[stringr::str_detect(cur_tree$tip.label, pattern = pattern)]
    cur_tree = drop.tip(cur_tree, ref_tips)
    cur_dir = dirname(cur_file)
    cur_filename = stringr::str_replace(basename(cur_file), pattern = "align.tre", "_refdropped.tre")
    cur_path = file.path(cur_dir, cur_filename)
    
    if(!file.exists(cur_path)){
      write.tree(phy = cur_tree, file = cur_path)
    }
    trees.path = c(cur_path, trees.path)
    
  }
  
  return(trees.path)
}

# calculate_quartet_distances()
# Returns a list of quartet distances between all possible pairings of trees in list 1 and lsit 2.
# paths.trees1: List of tree paths
# paths.trees2: List of tree paths to compare against
# cores; number of cores to use
calculate_quartet_distances = function(paths.trees1, paths.trees2, cores = 20){
  
  print("Registering cores...")
  
  # Set up parallel backend
  num_cores <- cores
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  print("Registering cores...DONE")
  
  # Generate all unique tree pairs
  pair_indices <- expand.grid(seq_along(paths.trees1), seq_along(paths.trees2))
  
  # Parallel computation of quartet distances
  quartet_dists <- foreach(i = 1:nrow(pair_indices), .combine = c, .export = "quartet_distance") %dopar% {
    quartet_distance(
      path.tree1 = paths.trees1[pair_indices[i, 1]], 
      path.tree2 = paths.trees2[pair_indices[i, 2]]
    )
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  return(quartet_dists)
  
}

# calculate_quartet_distances_random_pairs()
# Returns a list of quartet distances between random pairings of trees in list 1 and lsit 2.
# paths.trees1: List of tree paths
# paths.trees2: List of tree paths to compare against
# cores; number of cores to use
calculate_quartet_distances_random_pairs = function(paths.trees1, paths.trees2, cores = 20){
  
  print("Registering cores...")
  
  # Set up parallel backend
  num_cores <- cores
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  print("Registering cores...DONE")
  
  # Generate all unique tree pairs
  pair_indices <- data.frame(
    trees1_idx = sample(seq_along(paths.trees1), replace = FALSE), 
    trees2_idx = sample(seq_along(paths.trees2), replace = FALSE)
  )
  
  # Parallel computation of quartet distances
  quartet_dists <- foreach(i = 1:nrow(pair_indices), .combine = c, .export = "quartet_distance") %dopar% {
    quartet_distance(
      path.tree1 = paths.trees1[pair_indices[i, 1]], 
      path.tree2 = paths.trees2[pair_indices[i, 2]]
    )
  }
  
  
  # Stop the cluster
  stopCluster(cl)
  
  return(quartet_dists)
  
}

# quartet_distance()
# Calculates quartet distance between two trees provided.
# path.tree1: Provide path to first tree
# path.tree2: Provide path to second tree
# Returns a quartet distance as calculated by quartet_dist from
# https://www.birc.au.dk/~cstorm/software/tqdist/
quartet_distance = function(path.tree1, path.tree2){
  quartet_distance = NA
  if(!file.exists(path.tree1) | !file.exists(path.tree2)){
    stop("Tree files do not exist.")
  }else{
    
    cmd = paste0(c("/nfs3/Sharpton_Lab/tmp/src/arnoldhk/bin/tqDist/bin/quartet_dist ",
                   path.tree1, 
                   " ",
                   path.tree2), sep = "", collapse = "")
    
    quartet_dist = system(cmd, intern = TRUE)
  }
  return(as.numeric(quartet_dist))
}


# read_trees
# Reads in a list of trees. If drop = TRUE (e.g. the case of guides), 
# tips matching the pattern supplied are dropped out of the tree
# Returns a list of trees.
# file_list: a list of file paths to trees to read
read_trees = function(file_list, drop = FALSE, pattern = "^REF."){
  trees = list()
  
  for(i in 1:length(file_list)){
    cur_tree = read.tree(file_list[i])
    
    if(drop){
      ref_tips = cur_tree$tip.label[stringr::str_detect(cur_tree$tip.label, pattern = pattern)]
      cur_tree = drop.tip(cur_tree, ref_tips)
    }
    trees[[basename(file_list[i])]]= cur_tree
    
  }
  return(trees)
}

# calculate_rf_distance()
# Check RF distance of all possible comparisons of tree i in list 1 to trees in list2. 
# trees1 and trees2: lists of trees. All trees should have the same tips present.
# Cores: the number of cores to use
calculate_rf_distance = function(trees1, trees2, cores = 20){
  
  if(length(trees1$tip.labels) != length(trees2$tip.labels)){
    stop("Trees are not of the same size.")
  }
  print("Registering cores...")
  
  # Set up parallel backend
  num_cores <- cores
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  print("Registering cores...DONE")
  
  # Generate all unique tree pairs
  pair_indices <- expand.grid(seq_along(trees1), seq_along(trees2))
  
  # Parallel computation of RF distances
  rf_distance <- foreach(i = 1:nrow(pair_indices), .combine = c, .packages = "phangorn") %dopar% {
    RF.dist(trees1[[pair_indices[i, 1]]], trees2[[pair_indices[i, 2]]], normalize = TRUE)
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  return(rf_distance)
  
}

# calaculate_rf_distance_random_pairs()
# Check RF distance between random comparisons of trees in list 1 to trees in list2. 
# trees1 and trees2: lists of trees. All trees should have the same tips present.
# Cores: the number of cores to use
calculate_rf_distance_random_pairs = function(trees1, trees2, cores = 20){
  
  if(length(trees1$tip.labels) != length(trees2$tip.labels)){
    stop("Trees are not of the same size.")
  }
  
  print("Registering cores...")
  
  # Set up parallel backend
  num_cores <- cores
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  print("Registering cores...DONE")
  
  # Generate all random pairings
  pair_indices <- data.frame(
    trees1_idx = sample(seq_along(trees1), replace = FALSE), 
    trees2_idx = sample(seq_along(trees2), replace = FALSE)
  )
  
  # Parallel computation of RF distances
  rf_distance <- foreach(i = 1:nrow(pair_indices), .combine = c, .packages = "phangorn") %dopar% {
    RF.dist(trees1[[pair_indices[i, 1]]], trees2[[pair_indices[i, 2]]], normalize = TRUE)
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  return(rf_distance)
  
}

# calaculate_rf_distance_self()
# Returns a list of RF distances of RF distances of all trees calculated against all others in the list.
# tree_list: a list of trees. 
# Cores: the number of cores to do.
calculate_rf_distance_self = function(tree_list, cores = 20){
  
  # Set up parallel backend
  num_cores <- cores
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # Generate all unique tree pairs
  pair_indices <- combn(seq_along(tree_list), 2, simplify = FALSE)
  pair_indices
  
  # Parallel computation of RF distances
  rf_distance <- foreach(pair = pair_indices, .combine = c, .packages = "phangorn") %dopar% {
    RF.dist(tree_list[[pair[1]]], tree_list[[pair[2]]])
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  return(rf_distance)
  
}

