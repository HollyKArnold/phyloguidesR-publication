################################################################################
#AUTHOR: ARNOLD
#DAY: Jan 24th, 2025
#SCRIPT: experimental_demo_mast_random_pairs.R
#DESCRIPTION: Calculate MAST 

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
# 
################################################################################

setwd(dir.demo)
dir.demos = directories[c("short_guides_optimized", "short", "long", "long_control")] 
dir.demos
long = read_trees(file.path(dir.demos["long"], file_list = list.files(path = dir.demos["long"], pattern = ".tre")))
long_control = read_trees(file.path(dir.demos["long_control"], file_list = list.files(path = dir.demos["long_control"], pattern = ".tre")))
short = read_trees(file.path(dir.demos["short"], file_list = list.files(path = dir.demos["short"], pattern = ".tre")))
short_guides = read_trees(file.path(dir.demos["short_guides_optimized"], file_list = list.files(path = dir.demos["short_guides_optimized"], pattern = "align.tre")), drop = TRUE)

short_micro = calculate_mast_random_pairs(trees1 = long, trees2 = short, cores = CORES)
short_guides_micro = calculate_mast_random_pairs(trees1 = long, trees2 = short_guides, cores = CORES)
long_control_micro = calculate_mast_random_pairs(trees1 = long, trees2 = long_control, cores = CORES)

expeirmental_demo_mast_random_pairs = data.frame("short" = short_micro,
                                                 "short_guides" = short_guides_micro,
                                                 "long_control_micro" = long_control_micro)

setwd(dir.out)
saveRDS(file = "experimental_demo_mast_random_pairs.rds", expeirmental_demo_mast_random_pairs, compress = TRUE)
