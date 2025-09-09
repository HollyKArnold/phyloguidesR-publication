# Author: Arnold
# Sim 1: MST distance.
################################################################################
# SET VARIABLES
################################################################################
CORES = 100
SAVE_FILE_NAME = "hvr_msd_random_pairs.rds"
################################################################################
# LOAD LIBRARIES
################################################################################
library(phangorn)
library(parallel)
library(foreach)
library(doParallel)
library(reshape2)
library(progressr)
library(tidyr)
library(phylotools)
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

print("Working on V1...")
v1 = calculate_msd_distance_random_pairs(trees1 = VFull, trees2 = V1, cores = CORES)
print("Working on V1...DONE")

print("Working on V2...")
v2 = calculate_msd_distance_random_pairs(trees1 = VFull, trees2 = V2, cores = CORES)
print("Working on V2...DONE")

print("Working on V3...")
v3 = calculate_msd_distance_random_pairs(trees1 = VFull, trees2 = V3, cores = CORES)
print("Working on V3...DONE")

print("Working on V4...")
v4 = calculate_msd_distance_random_pairs(trees1 = VFull, trees2 = V4, cores = CORES)
print("Working on V4...DONE")

print("Working on V5...")
v5 = calculate_msd_distance_random_pairs(trees1 = VFull, trees2 = V5, cores = CORES)
print("Working on V5...DONE")

print("Working on V6...")
v6 = calculate_msd_distance_random_pairs(trees1 = VFull, trees2 = V6, cores = CORES)
print("Working on V6...DONE")

print("Working on V7...")
v7 = calculate_msd_distance_random_pairs(trees1 = VFull, trees2 = V7, cores = CORES)
print("Working on V7...DONE")

print("Working on V8...")
v8 = calculate_msd_distance_random_pairs(trees1 = VFull, trees2 = V8, cores = CORES)
print("Working on V8...DONE")

print("Working on V9...")
v9 = calculate_msd_distance_random_pairs(trees1 = VFull, trees2 = V9, cores = CORES)
print("Working on V9...DONE")

print("Working on VFull_Guides...")
vfull_guides = calculate_msd_distance_random_pairs(trees1 = VFull, trees2 = VFull_guides, cores = CORES)
print("Working on VFull_Guides...DONE")

print("Working on VFull_control...")
vfull = calculate_msd_distance_random_pairs(trees1 = VFull, trees2 = VFull_control, cores = CORES)
print("Working on VFull_control...DONE")

print("Working on V1 Guides...")
v1_guides = calculate_msd_distance_random_pairs(trees1 = VFull, trees2 = V1_guides, cores = CORES)
print("Working on V1 Guides...DONE")

print("Working on V2 Guides...")
v2_guides = calculate_msd_distance_random_pairs(trees1 = VFull, trees2 = V2_guides, cores = CORES)
print("Working on V2 Guides...DONE")

print("Working on V3 Guides...")
v3_guides = calculate_msd_distance_random_pairs(trees1 = VFull, trees2 = V3_guides, cores = CORES)
print("Working on V3 Guides...DONE")

print("Working on V4 Guides...")
v4_guides = calculate_msd_distance_random_pairs(trees1 = VFull, trees2 = V4_guides, cores = CORES)
print("Working on V4 Guides...DONE")

print("Working on V5 Guides...")
v5_guides = calculate_msd_distance_random_pairs(trees1 = VFull, trees2 = V5_guides, cores = CORES)
print("Working on V5 Guides...DONE")

print("Working on V6 Guides...")
v6_guides = calculate_msd_distance_random_pairs(trees1 = VFull, trees2 = V6_guides, cores = CORES)
print("Working on V6 Guides...DONE")

print("Working on V7 Guides...")
v7_guides = calculate_msd_distance_random_pairs(trees1 = VFull, trees2 = V7_guides, cores = CORES)
print("Working on V7 Guides...DONE")

print("Working on V8 Guides...")
v8_guides = calculate_msd_distance_random_pairs(trees1 = VFull, trees2 = V8_guides, cores = CORES)
print("Working on V8 Guides...DONE")

print("Working on V9 Guides...")
v9_guides = calculate_msd_distance_random_pairs(trees1 = VFull, trees2 = V9_guides, cores = CORES)
print("Working on V9 Guides...DONE")

# # Make Data Frame Of normalized RF
hvr =
  data.frame("V1" = v1,
             "V2" = v2,
             "V3" = v3,
             "V4" = v4,
             "V5" = v5,
             "V6" = v6,
             "V7" = v7,
             "V8" = v8,
             "V9" = v9,
             "V1_guides" = v1_guides,
             "V2_guides" = v2_guides,
             "V3_guides" = v3_guides,
             "V4_guides" = v4_guides,
             "V5_guides" = v5_guides,
             "V6_guides" = v6_guides,
             "V7_guides" = v7_guides,
             "V8_guides" = v8_guides,
             "V9_guides" = v9_guides,
             "VFull" = vfull,
             "VFull_guides" = vfull_guides)


setwd(dir.out)
saveRDS(hvr, file = SAVE_FILE_NAME, compress = TRUE)
