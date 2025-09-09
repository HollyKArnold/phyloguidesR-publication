# Author Arnold
# Sim 1: Path distance short-read and long-read trees.
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
TREE_SIZE = 6846
RF_MAX = 2*(TREE_SIZE -3)
QUARTET_MAX = choose(TREE_SIZE, 4)

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
length(VFull)
VFull_guides = read_trees(file.path(dir.hvrs["VFull_guides"], file_list = list.files(path = dir.hvrs["VFull_guides"], pattern = "align.tre")), drop = TRUE)
length(VFull_guides)
VFull_control = read_trees(file.path(dir.hvrs["VFull_control"], file_list = list.files(path = dir.hvrs["VFull_control"], pattern = "align.tre")))
length(VFull_control)

V1 = read_trees(file.path(dir.hvrs["V1_HVR"], file_list = list.files(path = dir.hvrs["V1_HVR"], pattern = ".tre")))
length(V1)
V2 = read_trees(file.path(dir.hvrs["V2_HVR"], file_list = list.files(path = dir.hvrs["V2_HVR"], pattern = ".tre")))
length(V2)
V3 = read_trees(file.path(dir.hvrs["V3_HVR"], file_list = list.files(path = dir.hvrs["V3_HVR"], pattern = ".tre")))
length(V3)
V4 = read_trees(file.path(dir.hvrs["V4_HVR"], file_list = list.files(path = dir.hvrs["V4_HVR"], pattern = ".tre")))
length(V4)
V5 = read_trees(file.path(dir.hvrs["V5_HVR"], file_list = list.files(path = dir.hvrs["V5_HVR"], pattern = ".tre")))
length(V5)
V6 = read_trees(file.path(dir.hvrs["V6_HVR"], file_list = list.files(path = dir.hvrs["V6_HVR"], pattern = ".tre")))
length(V6)
V7 = read_trees(file.path(dir.hvrs["V7_HVR"], file_list = list.files(path = dir.hvrs["V7_HVR"], pattern = ".tre")))
length(V7)
V8 = read_trees(file.path(dir.hvrs["V8_HVR"], file_list = list.files(path = dir.hvrs["V8_HVR"], pattern = ".tre")))
length(V8)
V9 = read_trees(file.path(dir.hvrs["V9_HVR"], file_list = list.files(path = dir.hvrs["V9_HVR"], pattern = ".tre")))
length(V9)

V1_guides = read_trees(file.path(dir.hvrs["V1_HVR_guides"], file_list = list.files(path = dir.hvrs["V1_HVR_guides"], pattern = "align.tre")), drop = TRUE)
length(V1_guides)
V2_guides = read_trees(file.path(dir.hvrs["V2_HVR_guides"], file_list = list.files(path = dir.hvrs["V2_HVR_guides"], pattern = "align.tre")), drop = TRUE)
length(V2_guides)
V3_guides = read_trees(file.path(dir.hvrs["V3_HVR_guides"], file_list = list.files(path = dir.hvrs["V3_HVR_guides"], pattern = "align.tre")), drop = TRUE)
length(V3_guides)
V4_guides = read_trees(file.path(dir.hvrs["V4_HVR_guides"], file_list = list.files(path = dir.hvrs["V4_HVR_guides"], pattern = "align.tre")), drop = TRUE)
length(V4_guides)
V5_guides = read_trees(file.path(dir.hvrs["V5_HVR_guides"], file_list = list.files(path = dir.hvrs["V5_HVR_guides"], pattern = "align.tre")), drop = TRUE)
length(V5_guides)
V6_guides = read_trees(file.path(dir.hvrs["V6_HVR_guides"], file_list = list.files(path = dir.hvrs["V6_HVR_guides"], pattern = "align.tre")), drop = TRUE)
length(V6_guides)
V7_guides = read_trees(file.path(dir.hvrs["V7_HVR_guides"], file_list = list.files(path = dir.hvrs["V7_HVR_guides"], pattern = "align.tre")), drop = TRUE)
length(V7_guides)
V8_guides = read_trees(file.path(dir.hvrs["V8_HVR_guides"], file_list = list.files(path = dir.hvrs["V8_HVR_guides"], pattern = "align.tre")), drop = TRUE)
length(V8_guides)
V9_guides = read_trees(file.path(dir.hvrs["V9_HVR_guides"], file_list = list.files(path = dir.hvrs["V9_HVR_guides"], pattern = "align.tre")), drop = TRUE)
length(V9_guides)

v1_path = calculate_path_dist_random_pairs(trees1 = VFull, trees2 = V1)
v2_path = calculate_path_dist_random_pairs(trees1 = VFull, trees2 = V2)
v3_path = calculate_path_dist_random_pairs(trees1 = VFull, trees2 = V3)
v4_path = calculate_path_dist_random_pairs(trees1 = VFull, trees2 = V4)
v5_path = calculate_path_dist_random_pairs(trees1 = VFull, trees2 = V5)
v6_path = calculate_path_dist_random_pairs(trees1 = VFull, trees2 = V6)
v7_path = calculate_path_dist_random_pairs(trees1 = VFull, trees2 = V7)
v8_path = calculate_path_dist_random_pairs(trees1 = VFull, trees2 = V8)
v9_path = calculate_path_dist_random_pairs(trees1 = VFull, trees2 = V9)

vfull_guides_path = calculate_path_dist_random_pairs(trees1 = VFull, trees2 = VFull_guides, cores = CORES)
vfull_path = calculate_path_dist_random_pairs(trees1 = VFull, trees2 = VFull_control, cores = CORES)

v1_guides_path = calculate_path_dist_random_pairs(trees1 = VFull, trees2 = V1_guides, cores = CORES)
v2_guides_path = calculate_path_dist_random_pairs(trees1 = VFull, trees2 = V2_guides, cores = CORES)
v3_guides_path = calculate_path_dist_random_pairs(trees1 = VFull, trees2 = V3_guides, cores = CORES)
v4_guides_path = calculate_path_dist_random_pairs(trees1 = VFull, trees2 = V4_guides, cores = CORES)
v5_guides_path = calculate_path_dist_random_pairs(trees1 = VFull, trees2 = V5_guides, cores = CORES)
v6_guides_path = calculate_path_dist_random_pairs(trees1 = VFull, trees2 = V6_guides, cores = CORES)
v7_guides_path = calculate_path_dist_random_pairs(trees1 = VFull, trees2 = V7_guides, cores = CORES)
v8_guides_path = calculate_path_dist_random_pairs(trees1 = VFull, trees2 = V8_guides, cores = CORES)
v9_guides_path = calculate_path_dist_random_pairs(trees1 = VFull, trees2 = V9_guides, cores = CORES)


# Make Data Frame Of normalized RF
hvr_path =
  data.frame("V1" = v1_path,
             "V2" = v2_path,
             "V3" = v3_path,
             "V4" = v4_path,
             "V5" = v5_path,
             "V6" = v6_path,
             "V7" = v7_path,
             "V8" = v8_path,
             "V9" = v9_path,
             "VFull" = vfull_path,
             "V1_guides" = v1_guides_path,
             "V2_guides" = v2_guides_path,
             "V3_guides" = v3_guides_path,
             "V4_guides" = v4_guides_path,
             "V5_guides" = v5_guides_path,
             "V6_guides" = v6_guides_path,
             "V7_guides" = v7_guides_path,
             "V8_guides" = v8_guides_path,
             "V9_guides" = v9_guides_path,
             "VFull_guides" = vfull_guides_path)


setwd(dir.out)
saveRDS(object = hvr_path, "hvr_path_random_pairings.rds", compress = TRUE)
