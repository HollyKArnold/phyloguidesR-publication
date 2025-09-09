################################################################################
#AUTHOR: ARNOLD
#DAY: Jan 24th, 2025
#SCRIPT: sim_analysis_nanostructure_random_pairs.R
#DESCRIPTION: Sim 1: Calculate similarity of nanostructures 2, ... , 10

################################################################################
################################################################################
# SET CONSTANTS
################################################################################
CORES = 100

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
dir.sim.out = file.path(HOME, "2025_hvr_guide_phylogenetic_integration/out/sim_bp/")

setwd(HOME_TMP)
dir.simulation = file.path(HOME_TMP, "length_simulations/")
dir.mix = file.path(HOME_TMP, "mix_experiment/")

setwd(dir.scripts)
source("sim_analysis_functions_V2.R")

################################################################################
# SET CONSTANTS
################################################################################
# Check that all the tree files are present. 
directories <- c(list.dirs(dir.mix, full.names = TRUE, recursive = FALSE),
                 list.dirs(dir.simulation, full.names = TRUE, recursive = FALSE))
names(directories) = basename(directories)
basename(directories)

################################################################################
# 0. SIM PREPREP: READ INPUT FILES
################################################################################

setwd(dir.simulation)
dir.hvrs = directories[c("V1_HVR", "V2_HVR", "V3_HVR", "V4_HVR", 
                         "V5_HVR", "V6_HVR", "V7_HVR", "V8_HVR", "V9_HVR", 
                         "V1_HVR_guides", "V2_HVR_guides", "V3_HVR_guides", 
                         "V4_HVR_guides", "V5_HVR_guides", "V6_HVR_guides", 
                         "V7_HVR_guides", "V8_HVR_guides", "V9_HVR_guides", 
                         "VFull", "VFull_guides", "VFull_control")]


VFull = read_trees(file.path(dir.hvrs["VFull"], file_list = list.files(path = dir.hvrs["VFull"], pattern = ".tre")))
VFull_guides = read_trees(file.path(dir.hvrs["VFull_guides"], file_list = list.files(path = dir.hvrs["VFull_guides"], pattern = "align.tre")), drop = TRUE)
VFull_control = read_trees(file.path(dir.hvrs["VFull_control"], file_list = list.files(path = dir.hvrs["VFull_control"], pattern = "align.tre")))

V1 = read_trees(file.path(dir.hvrs["V1_HVR"], file_list = list.files(path = dir.hvrs["V1_HVR"], pattern = ".tre")))
V2 = read_trees(file.path(dir.hvrs["V2_HVR"], file_list = list.files(path = dir.hvrs["V2_HVR"], pattern = ".tre")))
V3 = read_trees(file.path(dir.hvrs["V3_HVR"], file_list = list.files(path = dir.hvrs["V3_HVR"], pattern = ".tre")))
V4 = read_trees(file.path(dir.hvrs["V4_HVR"], file_list = list.files(path = dir.hvrs["V4_HVR"], pattern = ".tre")))
V5 = read_trees(file.path(dir.hvrs["V5_HVR"], file_list = list.files(path = dir.hvrs["V5_HVR"], pattern = ".tre")))
V6 = read_trees(file.path(dir.hvrs["V6_HVR"], file_list = list.files(path = dir.hvrs["V6_HVR"], pattern = ".tre")))
V7 = read_trees(file.path(dir.hvrs["V7_HVR"], file_list = list.files(path = dir.hvrs["V7_HVR"], pattern = ".tre")))
V8 = read_trees(file.path(dir.hvrs["V8_HVR"], file_list = list.files(path = dir.hvrs["V8_HVR"], pattern = ".tre")))
V9 = read_trees(file.path(dir.hvrs["V9_HVR"], file_list = list.files(path = dir.hvrs["V9_HVR"], pattern = ".tre")))

V1_guides = read_trees(file.path(dir.hvrs["V1_HVR_guides"], file_list = list.files(path = dir.hvrs["V1_HVR_guides"], pattern = "align.tre")), drop = TRUE)
V2_guides = read_trees(file.path(dir.hvrs["V2_HVR_guides"], file_list = list.files(path = dir.hvrs["V2_HVR_guides"], pattern = "align.tre")), drop = TRUE)
V3_guides = read_trees(file.path(dir.hvrs["V3_HVR_guides"], file_list = list.files(path = dir.hvrs["V3_HVR_guides"], pattern = "align.tre")), drop = TRUE)
V4_guides = read_trees(file.path(dir.hvrs["V4_HVR_guides"], file_list = list.files(path = dir.hvrs["V4_HVR_guides"], pattern = "align.tre")), drop = TRUE)
V5_guides = read_trees(file.path(dir.hvrs["V5_HVR_guides"], file_list = list.files(path = dir.hvrs["V5_HVR_guides"], pattern = "align.tre")), drop = TRUE)
V6_guides = read_trees(file.path(dir.hvrs["V6_HVR_guides"], file_list = list.files(path = dir.hvrs["V6_HVR_guides"], pattern = "align.tre")), drop = TRUE)
V7_guides = read_trees(file.path(dir.hvrs["V7_HVR_guides"], file_list = list.files(path = dir.hvrs["V7_HVR_guides"], pattern = "align.tre")), drop = TRUE)
V8_guides = read_trees(file.path(dir.hvrs["V8_HVR_guides"], file_list = list.files(path = dir.hvrs["V8_HVR_guides"], pattern = "align.tre")), drop = TRUE)
V9_guides = read_trees(file.path(dir.hvrs["V9_HVR_guides"], file_list = list.files(path = dir.hvrs["V9_HVR_guides"], pattern = "align.tre")), drop = TRUE)

# 
# SUBTREE_SIZE = 2
# v1_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V1, cores = CORES, n = SUBTREE_SIZE)
# v2_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V2, cores = CORES, n = SUBTREE_SIZE)
# v3_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V3, cores = CORES, n = SUBTREE_SIZE)
# v4_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V4, cores = CORES, n = SUBTREE_SIZE)
# v5_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V5, cores = CORES, n = SUBTREE_SIZE)
# v6_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V6, cores = CORES, n = SUBTREE_SIZE)
# v7_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V7, cores = CORES, n = SUBTREE_SIZE)
# v8_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V8, cores = CORES, n = SUBTREE_SIZE)
# v9_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V9, cores = CORES, n = SUBTREE_SIZE)
# 
# v1_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V1_guides, cores = CORES, n = SUBTREE_SIZE)
# v2_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V2_guides, cores = CORES, n = SUBTREE_SIZE)
# v3_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V3_guides, cores = CORES, n = SUBTREE_SIZE)
# v4_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V4_guides, cores = CORES, n = SUBTREE_SIZE)
# v5_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V5_guides, cores = CORES, n = SUBTREE_SIZE)
# v6_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V6_guides, cores = CORES, n = SUBTREE_SIZE)
# v7_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V7_guides, cores = CORES, n = SUBTREE_SIZE)
# v8_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V8_guides, cores = CORES, n = SUBTREE_SIZE)
# v9_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V9_guides, cores = CORES, n = SUBTREE_SIZE)
# 
# vfull_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = VFull_control, cores = CORES, n = SUBTREE_SIZE)
# vfull_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = VFull_guides, cores = CORES, n = SUBTREE_SIZE)
# 
# s2 = as_tibble(data.frame("V1" = 1 - v1_subtree_size,
#                           "V2" = 1 - v2_subtree_size,
#                           "V3" = 1 - v3_subtree_size,
#                           "V4" = 1 - v4_subtree_size,
#                           "V5" = 1 - v5_subtree_size,
#                           "V6" = 1 - v6_subtree_size,
#                           "V7" = 1 - v7_subtree_size,
#                           "V8" = 1 - v8_subtree_size,
#                           "V9" = 1 - v9_subtree_size,
#                           "V1_guides" = 1 - v1_guides_subtree_size,
#                           "V2_guides" = 1 - v2_guides_subtree_size,
#                           "V3_guides" = 1 - v3_guides_subtree_size,
#                           "V4_guides" = 1 - v4_guides_subtree_size,
#                           "V5_guides" = 1 - v5_guides_subtree_size,
#                           "V6_guides" = 1 - v6_guides_subtree_size,
#                           "V7_guides" = 1 - v7_guides_subtree_size,
#                           "V8_guides" = 1 - v8_guides_subtree_size,
#                           'V9_guides' = 1 - v9_guides_subtree_size,
#                           "C" = 1 - vfull_subtree_size,
#                           "C_guides" = 1 - vfull_guides_subtree_size))
# s2
# setwd(dir.out)
# saveRDS(file = "hvr_s2_random_pairs.rds", s2, compress = TRUE)
# 
# 
# SUBTREE_SIZE = 3
# v1_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V1, cores = CORES, n = SUBTREE_SIZE)
# v2_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V2, cores = CORES, n = SUBTREE_SIZE)
# v3_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V3, cores = CORES, n = SUBTREE_SIZE)
# v4_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V4, cores = CORES, n = SUBTREE_SIZE)
# v5_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V5, cores = CORES, n = SUBTREE_SIZE)
# v6_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V6, cores = CORES, n = SUBTREE_SIZE)
# v7_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V7, cores = CORES, n = SUBTREE_SIZE)
# v8_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V8, cores = CORES, n = SUBTREE_SIZE)
# v9_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V9, cores = CORES, n = SUBTREE_SIZE)
# 
# v1_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V1_guides, cores = CORES, n = SUBTREE_SIZE)
# v2_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V2_guides, cores = CORES, n = SUBTREE_SIZE)
# v3_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V3_guides, cores = CORES, n = SUBTREE_SIZE)
# v4_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V4_guides, cores = CORES, n = SUBTREE_SIZE)
# v5_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V5_guides, cores = CORES, n = SUBTREE_SIZE)
# v6_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V6_guides, cores = CORES, n = SUBTREE_SIZE)
# v7_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V7_guides, cores = CORES, n = SUBTREE_SIZE)
# v8_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V8_guides, cores = CORES, n = SUBTREE_SIZE)
# v9_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V9_guides, cores = CORES, n = SUBTREE_SIZE)
# 
# vfull_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = VFull_control, cores = CORES, n = SUBTREE_SIZE)
# vfull_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = VFull_guides, cores = CORES, n = SUBTREE_SIZE)
# 
# s3 = as_tibble(data.frame("V1" = 1 - v1_subtree_size,
#                           "V2" = 1 - v2_subtree_size,
#                           "V3" = 1 - v3_subtree_size,
#                           "V4" = 1 - v4_subtree_size,
#                           "V5" = 1 - v5_subtree_size,
#                           "V6" = 1 - v6_subtree_size,
#                           "V7" = 1 - v7_subtree_size,
#                           "V8" = 1 - v8_subtree_size,
#                           "V9" = 1 - v9_subtree_size,
#                           "V1_guides" = 1 - v1_guides_subtree_size,
#                           "V2_guides" = 1 - v2_guides_subtree_size,
#                           "V3_guides" = 1 - v3_guides_subtree_size,
#                           "V4_guides" = 1 - v4_guides_subtree_size,
#                           "V5_guides" = 1 - v5_guides_subtree_size,
#                           "V6_guides" = 1 - v6_guides_subtree_size,
#                           "V7_guides" = 1 - v7_guides_subtree_size,
#                           "V8_guides" = 1 - v8_guides_subtree_size,
#                           'V9_guides' = 1 - v9_guides_subtree_size, 
#                           "C" = 1 - vfull_subtree_size,
#                           "C_guides" = 1 - vfull_guides_subtree_size))
# s3
# setwd(dir.out)
# saveRDS(file = "hvr_s3_random_pairs.rds", s3, compress = TRUE)


SUBTREE_SIZE = 4
v1_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V1, cores = CORES, n = SUBTREE_SIZE)
v2_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V2, cores = CORES, n = SUBTREE_SIZE)
v3_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V3, cores = CORES, n = SUBTREE_SIZE)
v4_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V4, cores = CORES, n = SUBTREE_SIZE)
v5_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V5, cores = CORES, n = SUBTREE_SIZE)
v6_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V6, cores = CORES, n = SUBTREE_SIZE)
v7_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V7, cores = CORES, n = SUBTREE_SIZE)
v8_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V8, cores = CORES, n = SUBTREE_SIZE)
v9_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V9, cores = CORES, n = SUBTREE_SIZE)

v1_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V1_guides, cores = CORES, n = SUBTREE_SIZE)
v2_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V2_guides, cores = CORES, n = SUBTREE_SIZE)
v3_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V3_guides, cores = CORES, n = SUBTREE_SIZE)
v4_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V4_guides, cores = CORES, n = SUBTREE_SIZE)
v5_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V5_guides, cores = CORES, n = SUBTREE_SIZE)
v6_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V6_guides, cores = CORES, n = SUBTREE_SIZE)
v7_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V7_guides, cores = CORES, n = SUBTREE_SIZE)
v8_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V8_guides, cores = CORES, n = SUBTREE_SIZE)
v9_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V9_guides, cores = CORES, n = SUBTREE_SIZE)

vfull_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = VFull_control, cores = CORES, n = SUBTREE_SIZE)
vfull_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = VFull_guides, cores = CORES, n = SUBTREE_SIZE)

s4 = as_tibble(data.frame("V1" = 1 - v1_subtree_size,
                          "V2" = 1 - v2_subtree_size,
                          "V3" = 1 - v3_subtree_size,
                          "V4" = 1 - v4_subtree_size,
                          "V5" = 1 - v5_subtree_size,
                          "V6" = 1 - v6_subtree_size,
                          "V7" = 1 - v7_subtree_size,
                          "V8" = 1 - v8_subtree_size,
                          "V9" = 1 - v9_subtree_size,
                          "V1_guides" = 1 - v1_guides_subtree_size,
                          "V2_guides" = 1 - v2_guides_subtree_size,
                          "V3_guides" = 1 - v3_guides_subtree_size,
                          "V4_guides" = 1 - v4_guides_subtree_size,
                          "V5_guides" = 1 - v5_guides_subtree_size,
                          "V6_guides" = 1 - v6_guides_subtree_size,
                          "V7_guides" = 1 - v7_guides_subtree_size,
                          "V8_guides" = 1 - v8_guides_subtree_size,
                          'V9_guides' = 1 - v9_guides_subtree_size, 
                          "C" = 1 - vfull_subtree_size,
                          "C_guides" = 1 - vfull_guides_subtree_size))
s4
setwd(dir.out)
saveRDS(file = "hvr_s4_random_pairs.rds", s4, compress = TRUE)


SUBTREE_SIZE = 5
v1_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V1, cores = CORES, n = SUBTREE_SIZE)
v2_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V2, cores = CORES, n = SUBTREE_SIZE)
v3_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V3, cores = CORES, n = SUBTREE_SIZE)
v4_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V4, cores = CORES, n = SUBTREE_SIZE)
v5_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V5, cores = CORES, n = SUBTREE_SIZE)
v6_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V6, cores = CORES, n = SUBTREE_SIZE)
v7_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V7, cores = CORES, n = SUBTREE_SIZE)
v8_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V8, cores = CORES, n = SUBTREE_SIZE)
v9_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V9, cores = CORES, n = SUBTREE_SIZE)

v1_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V1_guides, cores = CORES, n = SUBTREE_SIZE)
v2_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V2_guides, cores = CORES, n = SUBTREE_SIZE)
v3_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V3_guides, cores = CORES, n = SUBTREE_SIZE)
v4_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V4_guides, cores = CORES, n = SUBTREE_SIZE)
v5_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V5_guides, cores = CORES, n = SUBTREE_SIZE)
v6_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V6_guides, cores = CORES, n = SUBTREE_SIZE)
v7_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V7_guides, cores = CORES, n = SUBTREE_SIZE)
v8_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V8_guides, cores = CORES, n = SUBTREE_SIZE)
v9_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V9_guides, cores = CORES, n = SUBTREE_SIZE)

vfull_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = VFull_control, cores = CORES, n = SUBTREE_SIZE)
vfull_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = VFull_guides, cores = CORES, n = SUBTREE_SIZE)

s5 = as_tibble(data.frame("V1" = 1 - v1_subtree_size,
                          "V2" = 1 - v2_subtree_size,
                          "V3" = 1 - v3_subtree_size,
                          "V4" = 1 - v4_subtree_size,
                          "V5" = 1 - v5_subtree_size,
                          "V6" = 1 - v6_subtree_size,
                          "V7" = 1 - v7_subtree_size,
                          "V8" = 1 - v8_subtree_size,
                          "V9" = 1 - v9_subtree_size,
                          "V1_guides" = 1 - v1_guides_subtree_size,
                          "V2_guides" = 1 - v2_guides_subtree_size,
                          "V3_guides" = 1 - v3_guides_subtree_size,
                          "V4_guides" = 1 - v4_guides_subtree_size,
                          "V5_guides" = 1 - v5_guides_subtree_size,
                          "V6_guides" = 1 - v6_guides_subtree_size,
                          "V7_guides" = 1 - v7_guides_subtree_size,
                          "V8_guides" = 1 - v8_guides_subtree_size,
                          'V9_guides' = 1 - v9_guides_subtree_size, 
                          "C" = 1 - vfull_subtree_size,
                          "C_guides" = 1 - vfull_guides_subtree_size))
s5
setwd(dir.out)
saveRDS(file = "hvr_s5_random_pairs.rds", s5, compress = TRUE)

# SUBTREE_SIZE = 6
# v1_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V1, cores = CORES, n = SUBTREE_SIZE)
# v2_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V2, cores = CORES, n = SUBTREE_SIZE)
# v3_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V3, cores = CORES, n = SUBTREE_SIZE)
# v4_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V4, cores = CORES, n = SUBTREE_SIZE)
# v5_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V5, cores = CORES, n = SUBTREE_SIZE)
# v6_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V6, cores = CORES, n = SUBTREE_SIZE)
# v7_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V7, cores = CORES, n = SUBTREE_SIZE)
# v8_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V8, cores = CORES, n = SUBTREE_SIZE)
# v9_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V9, cores = CORES, n = SUBTREE_SIZE)
# 
# v1_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V1_guides, cores = CORES, n = SUBTREE_SIZE)
# v2_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V2_guides, cores = CORES, n = SUBTREE_SIZE)
# v3_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V3_guides, cores = CORES, n = SUBTREE_SIZE)
# v4_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V4_guides, cores = CORES, n = SUBTREE_SIZE)
# v5_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V5_guides, cores = CORES, n = SUBTREE_SIZE)
# v6_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V6_guides, cores = CORES, n = SUBTREE_SIZE)
# v7_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V7_guides, cores = CORES, n = SUBTREE_SIZE)
# v8_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V8_guides, cores = CORES, n = SUBTREE_SIZE)
# v9_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V9_guides, cores = CORES, n = SUBTREE_SIZE)
# 
# vfull_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = VFull_control, cores = CORES, n = SUBTREE_SIZE)
# vfull_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = VFull_guides, cores = CORES, n = SUBTREE_SIZE)
# 
# s6 = as_tibble(data.frame("V1" = 1 - v1_subtree_size,
#                           "V2" = 1 - v2_subtree_size,
#                           "V3" = 1 - v3_subtree_size,
#                           "V4" = 1 - v4_subtree_size,
#                           "V5" = 1 - v5_subtree_size,
#                           "V6" = 1 - v6_subtree_size,
#                           "V7" = 1 - v7_subtree_size,
#                           "V8" = 1 - v8_subtree_size,
#                           "V9" = 1 - v9_subtree_size,
#                           "V1_guides" = 1 - v1_guides_subtree_size,
#                           "V2_guides" = 1 - v2_guides_subtree_size,
#                           "V3_guides" = 1 - v3_guides_subtree_size,
#                           "V4_guides" = 1 - v4_guides_subtree_size,
#                           "V5_guides" = 1 - v5_guides_subtree_size,
#                           "V6_guides" = 1 - v6_guides_subtree_size,
#                           "V7_guides" = 1 - v7_guides_subtree_size,
#                           "V8_guides" = 1 - v8_guides_subtree_size,
#                           'V9_guides' = 1 - v9_guides_subtree_size, 
#                           "C" = 1 - vfull_subtree_size,
#                           "C_guides" = 1 - vfull_guides_subtree_size))
# s6
# setwd(dir.out)
# saveRDS(file = "hvr_s6_random_pairs.rds", s6, compress = TRUE)
# 
# SUBTREE_SIZE = 7
# v1_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V1, cores = CORES, n = SUBTREE_SIZE)
# v2_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V2, cores = CORES, n = SUBTREE_SIZE)
# v3_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V3, cores = CORES, n = SUBTREE_SIZE)
# v4_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V4, cores = CORES, n = SUBTREE_SIZE)
# v5_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V5, cores = CORES, n = SUBTREE_SIZE)
# v6_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V6, cores = CORES, n = SUBTREE_SIZE)
# v7_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V7, cores = CORES, n = SUBTREE_SIZE)
# v8_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V8, cores = CORES, n = SUBTREE_SIZE)
# v9_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V9, cores = CORES, n = SUBTREE_SIZE)
# 
# v1_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V1_guides, cores = CORES, n = SUBTREE_SIZE)
# v2_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V2_guides, cores = CORES, n = SUBTREE_SIZE)
# v3_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V3_guides, cores = CORES, n = SUBTREE_SIZE)
# v4_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V4_guides, cores = CORES, n = SUBTREE_SIZE)
# v5_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V5_guides, cores = CORES, n = SUBTREE_SIZE)
# v6_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V6_guides, cores = CORES, n = SUBTREE_SIZE)
# v7_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V7_guides, cores = CORES, n = SUBTREE_SIZE)
# v8_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V8_guides, cores = CORES, n = SUBTREE_SIZE)
# v9_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V9_guides, cores = CORES, n = SUBTREE_SIZE)
# 
# vfull_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = VFull_control, cores = CORES, n = SUBTREE_SIZE)
# vfull_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = VFull_guides, cores = CORES, n = SUBTREE_SIZE)
# 
# s7 = as_tibble(data.frame("V1" = 1 - v1_subtree_size,
#                           "V2" = 1 - v2_subtree_size,
#                           "V3" = 1 - v3_subtree_size,
#                           "V4" = 1 - v4_subtree_size,
#                           "V5" = 1 - v5_subtree_size,
#                           "V6" = 1 - v6_subtree_size,
#                           "V7" = 1 - v7_subtree_size,
#                           "V8" = 1 - v8_subtree_size,
#                           "V9" = 1 - v9_subtree_size,
#                           "V1_guides" = 1 - v1_guides_subtree_size,
#                           "V2_guides" = 1 - v2_guides_subtree_size,
#                           "V3_guides" = 1 - v3_guides_subtree_size,
#                           "V4_guides" = 1 - v4_guides_subtree_size,
#                           "V5_guides" = 1 - v5_guides_subtree_size,
#                           "V6_guides" = 1 - v6_guides_subtree_size,
#                           "V7_guides" = 1 - v7_guides_subtree_size,
#                           "V8_guides" = 1 - v8_guides_subtree_size,
#                           'V9_guides' = 1 - v9_guides_subtree_size, 
#                           "C" = 1 - vfull_subtree_size,
#                           "C_guides" = 1 - vfull_guides_subtree_size))
# s7
# setwd(dir.out)
# saveRDS(file = "hvr_s7_random_pairs.rds", s7, compress = TRUE)
# 
# SUBTREE_SIZE = 8
# v1_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V1, cores = CORES, n = SUBTREE_SIZE)
# v2_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V2, cores = CORES, n = SUBTREE_SIZE)
# v3_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V3, cores = CORES, n = SUBTREE_SIZE)
# v4_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V4, cores = CORES, n = SUBTREE_SIZE)
# v5_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V5, cores = CORES, n = SUBTREE_SIZE)
# v6_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V6, cores = CORES, n = SUBTREE_SIZE)
# v7_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V7, cores = CORES, n = SUBTREE_SIZE)
# v8_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V8, cores = CORES, n = SUBTREE_SIZE)
# v9_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V9, cores = CORES, n = SUBTREE_SIZE)
# 
# v1_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V1_guides, cores = CORES, n = SUBTREE_SIZE)
# v2_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V2_guides, cores = CORES, n = SUBTREE_SIZE)
# v3_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V3_guides, cores = CORES, n = SUBTREE_SIZE)
# v4_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V4_guides, cores = CORES, n = SUBTREE_SIZE)
# v5_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V5_guides, cores = CORES, n = SUBTREE_SIZE)
# v6_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V6_guides, cores = CORES, n = SUBTREE_SIZE)
# v7_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V7_guides, cores = CORES, n = SUBTREE_SIZE)
# v8_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V8_guides, cores = CORES, n = SUBTREE_SIZE)
# v9_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V9_guides, cores = CORES, n = SUBTREE_SIZE)
# 
# vfull_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = VFull_control, cores = CORES, n = SUBTREE_SIZE)
# vfull_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = VFull_guides, cores = CORES, n = SUBTREE_SIZE)
# 
# s8 = as_tibble(data.frame("V1" = 1 - v1_subtree_size,
#                           "V2" = 1 - v2_subtree_size,
#                           "V3" = 1 - v3_subtree_size,
#                           "V4" = 1 - v4_subtree_size,
#                           "V5" = 1 - v5_subtree_size,
#                           "V6" = 1 - v6_subtree_size,
#                           "V7" = 1 - v7_subtree_size,
#                           "V8" = 1 - v8_subtree_size,
#                           "V9" = 1 - v9_subtree_size,
#                           "V1_guides" = 1 - v1_guides_subtree_size,
#                           "V2_guides" = 1 - v2_guides_subtree_size,
#                           "V3_guides" = 1 - v3_guides_subtree_size,
#                           "V4_guides" = 1 - v4_guides_subtree_size,
#                           "V5_guides" = 1 - v5_guides_subtree_size,
#                           "V6_guides" = 1 - v6_guides_subtree_size,
#                           "V7_guides" = 1 - v7_guides_subtree_size,
#                           "V8_guides" = 1 - v8_guides_subtree_size,
#                           'V9_guides' = 1 - v9_guides_subtree_size, 
#                           "C" = 1 - vfull_subtree_size,
#                           "C_guides" = 1 - vfull_guides_subtree_size))
# s8
# setwd(dir.out)
# saveRDS(file = "hvr_s8_random_pairs.rds", s8, compress = TRUE)
# 
# SUBTREE_SIZE = 9
# v1_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V1, cores = CORES, n = SUBTREE_SIZE)
# v2_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V2, cores = CORES, n = SUBTREE_SIZE)
# v3_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V3, cores = CORES, n = SUBTREE_SIZE)
# v4_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V4, cores = CORES, n = SUBTREE_SIZE)
# v5_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V5, cores = CORES, n = SUBTREE_SIZE)
# v6_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V6, cores = CORES, n = SUBTREE_SIZE)
# v7_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V7, cores = CORES, n = SUBTREE_SIZE)
# v8_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V8, cores = CORES, n = SUBTREE_SIZE)
# v9_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V9, cores = CORES, n = SUBTREE_SIZE)
# 
# v1_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V1_guides, cores = CORES, n = SUBTREE_SIZE)
# v2_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V2_guides, cores = CORES, n = SUBTREE_SIZE)
# v3_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V3_guides, cores = CORES, n = SUBTREE_SIZE)
# v4_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V4_guides, cores = CORES, n = SUBTREE_SIZE)
# v5_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V5_guides, cores = CORES, n = SUBTREE_SIZE)
# v6_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V6_guides, cores = CORES, n = SUBTREE_SIZE)
# v7_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V7_guides, cores = CORES, n = SUBTREE_SIZE)
# v8_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V8_guides, cores = CORES, n = SUBTREE_SIZE)
# v9_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V9_guides, cores = CORES, n = SUBTREE_SIZE)
# 
# vfull_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = VFull_control, cores = CORES, n = SUBTREE_SIZE)
# vfull_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = VFull_guides, cores = CORES, n = SUBTREE_SIZE)
# 
# s9 = as_tibble(data.frame("V1" = 1 - v1_subtree_size,
#                           "V2" = 1 - v2_subtree_size,
#                           "V3" = 1 - v3_subtree_size,
#                           "V4" = 1 - v4_subtree_size,
#                           "V5" = 1 - v5_subtree_size,
#                           "V6" = 1 - v6_subtree_size,
#                           "V7" = 1 - v7_subtree_size,
#                           "V8" = 1 - v8_subtree_size,
#                           "V9" = 1 - v9_subtree_size,
#                           "V1_guides" = 1 - v1_guides_subtree_size,
#                           "V2_guides" = 1 - v2_guides_subtree_size,
#                           "V3_guides" = 1 - v3_guides_subtree_size,
#                           "V4_guides" = 1 - v4_guides_subtree_size,
#                           "V5_guides" = 1 - v5_guides_subtree_size,
#                           "V6_guides" = 1 - v6_guides_subtree_size,
#                           "V7_guides" = 1 - v7_guides_subtree_size,
#                           "V8_guides" = 1 - v8_guides_subtree_size,
#                           'V9_guides' = 1 - v9_guides_subtree_size, 
#                           "C" = 1 - vfull_subtree_size,
#                           "C_guides" = 1 - vfull_guides_subtree_size))
# s9
# setwd(dir.out)
# saveRDS(file = "hvr_s9_random_pairs.rds", s9, compress = TRUE)
# 
# SUBTREE_SIZE = 10
# v1_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V1, cores = CORES, n = SUBTREE_SIZE)
# v2_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V2, cores = CORES, n = SUBTREE_SIZE)
# v3_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V3, cores = CORES, n = SUBTREE_SIZE)
# v4_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V4, cores = CORES, n = SUBTREE_SIZE)
# v5_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V5, cores = CORES, n = SUBTREE_SIZE)
# v6_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V6, cores = CORES, n = SUBTREE_SIZE)
# v7_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V7, cores = CORES, n = SUBTREE_SIZE)
# v8_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V8, cores = CORES, n = SUBTREE_SIZE)
# v9_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V9, cores = CORES, n = SUBTREE_SIZE)
# 
# v1_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V1_guides, cores = CORES, n = SUBTREE_SIZE)
# v2_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V2_guides, cores = CORES, n = SUBTREE_SIZE)
# v3_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V3_guides, cores = CORES, n = SUBTREE_SIZE)
# v4_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V4_guides, cores = CORES, n = SUBTREE_SIZE)
# v5_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V5_guides, cores = CORES, n = SUBTREE_SIZE)
# v6_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V6_guides, cores = CORES, n = SUBTREE_SIZE)
# v7_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V7_guides, cores = CORES, n = SUBTREE_SIZE)
# v8_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V8_guides, cores = CORES, n = SUBTREE_SIZE)
# v9_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = V9_guides, cores = CORES, n = SUBTREE_SIZE)
# 
# vfull_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = VFull_control, cores = CORES, n = SUBTREE_SIZE)
# vfull_guides_subtree_size = calculate_shared_subtrees_dist_random_pairs(trees1 = VFull, trees2 = VFull_guides, cores = CORES, n = SUBTREE_SIZE)
# 
# s10 = as_tibble(data.frame("V1" = 1 - v1_subtree_size,
#                            "V2" = 1 - v2_subtree_size,
#                            "V3" = 1 - v3_subtree_size,
#                            "V4" = 1 - v4_subtree_size,
#                            "V5" = 1 - v5_subtree_size,
#                            "V6" = 1 - v6_subtree_size,
#                            "V7" = 1 - v7_subtree_size,
#                            "V8" = 1 - v8_subtree_size,
#                            "V9" = 1 - v9_subtree_size,
#                            "V1_guides" = 1 - v1_guides_subtree_size,
#                            "V2_guides" = 1 - v2_guides_subtree_size,
#                            "V3_guides" = 1 - v3_guides_subtree_size,
#                            "V4_guides" = 1 - v4_guides_subtree_size,
#                            "V5_guides" = 1 - v5_guides_subtree_size,
#                            "V6_guides" = 1 - v6_guides_subtree_size,
#                            "V7_guides" = 1 - v7_guides_subtree_size,
#                            "V8_guides" = 1 - v8_guides_subtree_size,
#                            'V9_guides' = 1 - v9_guides_subtree_size, 
#                            "C" = 1 - vfull_subtree_size,
#                            "C_guides" = 1 - vfull_guides_subtree_size))
# s10
# setwd(dir.out)
# saveRDS(file = "hvr_s10_random_pairs.rds", s10, compress = TRUE)

