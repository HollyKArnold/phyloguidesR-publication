################################################################################
#AUTHOR: ARNOLD
#DAY: Jan 24th, 2025
#SCRIPT: sim_analysis_nanostructure_random_pairs.R
#DESCRIPTION: Calculate similarity of nanostructures 2, ... , 10

################################################################################
################################################################################
# SET CONSTANTS
################################################################################
CORES = 10

################################################################################
# LOAD LIBRARIES
################################################################################
library(phangorn)
library(parallel)
library(foreach)
library(doParallel)
library(reshape2)
library(ggplot2)
library(progressr)
library(tidyr)
library(phylotools)
library(ggpubr)
library(ape)
library(dplyr)


################################################################################
# SET DIRECTORIES
################################################################################
# Set Home
HOME = "/nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration"
HOME_TMP = "/nfs3/Sharpton_Lab/tmp/projects/arnoldh/2025_hvr_guide_phylogenetic_integration/"

# Make directories
setwd(HOME)
dir.scripts = file.path(HOME, "2025_hvr_guide_phylogenetic_integration/scripts/")
dir.ref = file.path(HOME, "reference_data/")
dir.out = file.path(HOME, "2025_hvr_guide_phylogenetic_integration/data_visulization/")
dir.sim.out = file.path(HOME, "2025_hvr_guide_phylogenetic_integration/out/experimental_demo//")

setwd(HOME_TMP)
dir.demo = file.path(HOME_TMP, "experimental_demo/")

setwd(dir.scripts)
source("sim_analysis_functions_V2.R")

################################################################################
# SET CONSTANTS
################################################################################

# Check that all the tree files are present. 
directories <- c(list.dirs(dir.demo, full.names = TRUE, recursive = FALSE))
names(directories) = basename(directories)
basename(directories)

################################################################################
# 0. SIM PREPREP: READ INPUT FILES
################################################################################

setwd(dir.demo)
dir.demos = directories[c("short_guides_optimized", "short", "long", "long_control")] 
dir.demos
long = read_trees(file.path(dir.demos["long"], file_list = list.files(path = dir.demos["long"], pattern = ".tre")))
long_control = read_trees(file.path(dir.demos["long_control"], file_list = list.files(path = dir.demos["long_control"], pattern = ".tre")))
short = read_trees(file.path(dir.demos["short"], file_list = list.files(path = dir.demos["short"], pattern = ".tre")))
short_guides = read_trees(file.path(dir.demos["short_guides_optimized"], file_list = list.files(path = dir.demos["short_guides_optimized"], pattern = "align.tre")), drop = TRUE)

# S2
short_subtree_size_s2 = calculate_shared_subtrees_dist_random_pairs(trees1 = long, trees2 = short, cores = CORES, n = 2)
short_guides_subtree_size_s2 = calculate_shared_subtrees_dist_random_pairs(trees1 = long, trees2 = short_guides, cores = CORES, n = 2)
long_control_subtree_size_s2 = calculate_shared_subtrees_dist_random_pairs(trees1 = long, trees2 = long_control, cores = CORES, n = 2)


experimental_demo_subtree_sizes = 
  data.frame("s2_short" = short_subtree_size_s2,
             "s2_short_guides" = short_guides_subtree_size_s2,
             "s2_long_control" = long_control_subtree_size_s2)

setwd(dir.out)
saveRDS(file = "experimental_demo_subtrees_random_pairs.rds", experimental_demo_subtree_sizes, compress = TRUE)

