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
#
################################################################################

setwd(dir.demo)
dir.demos = directories[c("short_guides_optimized", "short", "long", "long_control")] 
dir.demos
long = read_trees(file.path(dir.demos["long"], file_list = list.files(path = dir.demos["long"], pattern = ".tre")))
long_control = read_trees(file.path(dir.demos["long_control"], file_list = list.files(path = dir.demos["long_control"], pattern = ".tre")))
short = read_trees(file.path(dir.demos["short"], file_list = list.files(path = dir.demos["short"], pattern = ".tre")))
short_guides = read_trees(file.path(dir.demos["short_guides_optimized"], file_list = list.files(path = dir.demos["short_guides_optimized"], pattern = "align.tre")), drop = TRUE)


short_path = calculate_path_dist_random_pairs(trees1 = long, trees2 = short)
short_guides_path = calculate_path_dist_random_pairs(trees1 = long, trees2 = short_guides)
long_control = calculate_path_dist_random_pairs(trees1 = long, trees2 = long_control)

# Make Data Frame Of normalized RF
exp_demo_path =
  data.frame("short" = short_path,
             "short_guides" = short_guides_path,
             "long_control" = long_control) 

setwd(dir.out)
saveRDS(object = exp_demo_path, "experimental_demo_path_random_pairings.rds", compress = TRUE)
