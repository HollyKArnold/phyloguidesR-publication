################################################################################
#AUTHOR: ARNOLD
#DAY: Jan 24th, 2025
#SCRIPT: sim_analysis_hvr_mast.R
#DESCRIPTION: Calculate MAST 

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

v1_micro = calculate_mast_random_pairs(trees1 = VFull, trees2 = V1, cores = CORES)
v2_micro = calculate_mast_random_pairs(trees1 = VFull, trees2 = V2, cores = CORES)
v3_micro = calculate_mast_random_pairs(trees1 = VFull, trees2 = V3, cores = CORES)
v4_micro = calculate_mast_random_pairs(trees1 = VFull, trees2 = V4, cores = CORES)
v5_micro = calculate_mast_random_pairs(trees1 = VFull, trees2 = V5, cores = CORES)
v6_micro = calculate_mast_random_pairs(trees1 = VFull, trees2 = V6, cores = CORES)
v7_micro = calculate_mast_random_pairs(trees1 = VFull, trees2 = V7, cores = CORES)
v8_micro = calculate_mast_random_pairs(trees1 = VFull, trees2 = V8, cores = CORES)
v9_micro = calculate_mast_random_pairs(trees1 = VFull, trees2 = V9, cores = CORES)

v1_guides_micro = calculate_mast_random_pairs(trees1 = VFull, trees2 = V1_guides, cores = CORES)
v2_guides_micro = calculate_mast_random_pairs(trees1 = VFull, trees2 = V2_guides, cores = CORES)
v3_guides_micro = calculate_mast_random_pairs(trees1 = VFull, trees2 = V3_guides, cores = CORES)
v4_guides_micro = calculate_mast_random_pairs(trees1 = VFull, trees2 = V4_guides, cores = CORES)
v5_guides_micro = calculate_mast_random_pairs(trees1 = VFull, trees2 = V5_guides, cores = CORES)
v6_guides_micro = calculate_mast_random_pairs(trees1 = VFull, trees2 = V6_guides, cores = CORES)
v7_guides_micro = calculate_mast_random_pairs(trees1 = VFull, trees2 = V7_guides, cores = CORES)
v8_guides_micro = calculate_mast_random_pairs(trees1 = VFull, trees2 = V8_guides, cores = CORES)
v9_guides_micro = calculate_mast_random_pairs(trees1 = VFull, trees2 = V9_guides, cores = CORES)

vfull_micro = calculate_mast_random_pairs(trees1 = VFull, trees2 = VFull_control, cores = CORES)
vfull_guides_micro = calculate_mast_random_pairs(trees1 = VFull, trees2 = VFull_guides, cores = CORES)


hvr_mast_random_pairs = as_tibble(data.frame("V1" = v1_micro,
                                             "V2" = v2_micro,
                                             "V3" = v3_micro,
                                             "V4" = v4_micro,
                                             "V5" = v5_micro,
                                             "V6" = v6_micro,
                                             "V7" = v7_micro,
                                             "V8" = v8_micro,
                                             "V9" = v9_micro,
                                             "V1_guides" = v1_guides_micro,
                                             "V2_guides" = v2_guides_micro,
                                             "V3_guides" = v3_guides_micro,
                                             "V4_guides" = v4_guides_micro,
                                             "V5_guides" = v5_guides_micro,
                                             "V6_guides" = v6_guides_micro,
                                             "V7_guides" = v7_guides_micro,
                                             "V8_guides" = v8_guides_micro,
                                             'V9_guides' = v9_guides_micro, 
                                             "C" = vfull_micro,
                                             "C_guides" = vfull_guides_micro))

setwd(dir.out)
saveRDS(file = "hvr_mast_random_pairs.rds", hvr_mast_random_pairs, compress = TRUE)
