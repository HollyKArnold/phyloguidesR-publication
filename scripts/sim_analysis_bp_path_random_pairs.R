# Author: Arnold
# Sim 2: Sim_analysis_bp_path_random_pairs.R
# Calculate Path distance between short-read and full-length control trees across different short-read lengths.
N = 100 # The number of trees to look at per simulation
CORES = 50

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
dir.save = "/nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/save_image_scripts"
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

setwd(dir.simulation)
dir.hvrs = directories[c("V1_L50", "V1_L100", "V1_L200", "V1_L300", "V1_L400", "V1_L500",
                         "V2_L50", "V2_L100", "V2_L200", "V2_L300", "V2_L400", "V2_L500",
                         "V3_L50", "V3_L100", "V3_L200", "V3_L300", "V3_L400", "V3_L500",
                         "V4_L50", "V4_L100", "V4_L200", "V4_L300", "V4_L400", "V4_L500",
                         "V5_L50", "V5_L100", "V5_L200", "V5_L300", "V5_L400", "V5_L500",
                         "V6_L50", "V6_L100", "V6_L200", "V6_L300", "V6_L400", 
                         "V7_L50", "V7_L100", "V7_L200", "V7_L300",  
                         "V8_L50", "V8_L100", "V8_L200", 
                         
                         "V1_L50_guides", "V1_L100_guides", "V1_L200_guides", "V1_L300_guides", "V1_L400_guides", "V1_L500_guides",
                         "V2_L50_guides", "V2_L100_guides", "V2_L200_guides", "V2_L300_guides", "V2_L400_guides", "V2_L500_guides",
                         "V3_L50_guides", "V3_L100_guides", "V3_L200_guides", "V3_L300_guides", "V3_L400_guides", "V3_L500_guides",
                         "V4_L50_guides", "V4_L100_guides", "V4_L200_guides", "V4_L300_guides", "V4_L400_guides", "V4_L500_guides",
                         "V5_L50_guides", "V5_L100_guides", "V5_L200_guides", "V5_L300_guides", "V5_L400_guides", "V5_L500_guides",
                         "V6_L50_guides", "V6_L100_guides", "V6_L200_guides", "V6_L300_guides", "V6_L400_guides", 
                         "V7_L50_guides", "V7_L100_guides", "V7_L200_guides", "V7_L300_guides",  
                         "V8_L50_guides", "V8_L100_guides", "V8_L200_guides", 
                         "VFull", 
                         "VFull_guides", 
                         "VFull_control")]

print("Reading in files...")
VFull = read_trees(file.path(dir.hvrs["VFull"], file_list = list.files(path = dir.hvrs["VFull"], pattern = ".tre")))[1:N]
VFull_guides = read_trees(file.path(dir.hvrs["VFull_guides"], file_list = list.files(path = dir.hvrs["VFull_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
VFull_control = read_trees(file.path(dir.hvrs["VFull_control"], file_list = list.files(path = dir.hvrs["VFull_control"], pattern = "align.tre")))[1:N]

V1_L50 = read_trees(file.path(dir.hvrs["V1_L50"], file_list = list.files(path = dir.hvrs["V1_L50"], pattern = ".tre")))[1:N]
V2_L50 = read_trees(file.path(dir.hvrs["V2_L50"], file_list = list.files(path = dir.hvrs["V2_L50"], pattern = ".tre")))[1:N]
V3_L50 = read_trees(file.path(dir.hvrs["V3_L50"], file_list = list.files(path = dir.hvrs["V3_L50"], pattern = ".tre")))[1:N]
V4_L50 = read_trees(file.path(dir.hvrs["V4_L50"], file_list = list.files(path = dir.hvrs["V4_L50"], pattern = ".tre")))[1:N]
V5_L50 = read_trees(file.path(dir.hvrs["V5_L50"], file_list = list.files(path = dir.hvrs["V5_L50"], pattern = ".tre")))[1:N]
V6_L50 = read_trees(file.path(dir.hvrs["V6_L50"], file_list = list.files(path = dir.hvrs["V6_L50"], pattern = ".tre")))[1:N]
V7_L50 = read_trees(file.path(dir.hvrs["V7_L50"], file_list = list.files(path = dir.hvrs["V7_L50"], pattern = ".tre")))[1:N]
V8_L50 = read_trees(file.path(dir.hvrs["V8_L50"], file_list = list.files(path = dir.hvrs["V8_L50"], pattern = ".tre")))[1:N]

V1_L100 = read_trees(file.path(dir.hvrs["V1_L100"], file_list = list.files(path = dir.hvrs["V1_L100"], pattern = ".tre")))[1:N]
V2_L100 = read_trees(file.path(dir.hvrs["V2_L100"], file_list = list.files(path = dir.hvrs["V2_L100"], pattern = ".tre")))[1:N]
V3_L100 = read_trees(file.path(dir.hvrs["V3_L100"], file_list = list.files(path = dir.hvrs["V3_L100"], pattern = ".tre")))[1:N]
V4_L100 = read_trees(file.path(dir.hvrs["V4_L100"], file_list = list.files(path = dir.hvrs["V4_L100"], pattern = ".tre")))[1:N]
V5_L100 = read_trees(file.path(dir.hvrs["V5_L100"], file_list = list.files(path = dir.hvrs["V5_L100"], pattern = ".tre")))[1:N]
V6_L100 = read_trees(file.path(dir.hvrs["V6_L100"], file_list = list.files(path = dir.hvrs["V6_L100"], pattern = ".tre")))[1:N]
V7_L100 = read_trees(file.path(dir.hvrs["V7_L100"], file_list = list.files(path = dir.hvrs["V7_L100"], pattern = ".tre")))[1:N]
V8_L100 = read_trees(file.path(dir.hvrs["V8_L100"], file_list = list.files(path = dir.hvrs["V8_L100"], pattern = ".tre")))[1:N]

V1_L200 = read_trees(file.path(dir.hvrs["V1_L200"], file_list = list.files(path = dir.hvrs["V1_L200"], pattern = ".tre")))[1:N]
V2_L200 = read_trees(file.path(dir.hvrs["V2_L200"], file_list = list.files(path = dir.hvrs["V2_L200"], pattern = ".tre")))[1:N]
V3_L200 = read_trees(file.path(dir.hvrs["V3_L200"], file_list = list.files(path = dir.hvrs["V3_L200"], pattern = ".tre")))[1:N]
V4_L200 = read_trees(file.path(dir.hvrs["V4_L200"], file_list = list.files(path = dir.hvrs["V4_L200"], pattern = ".tre")))[1:N]
V5_L200 = read_trees(file.path(dir.hvrs["V5_L200"], file_list = list.files(path = dir.hvrs["V5_L200"], pattern = ".tre")))[1:N]
V6_L200 = read_trees(file.path(dir.hvrs["V6_L200"], file_list = list.files(path = dir.hvrs["V6_L200"], pattern = ".tre")))[1:N]
V7_L200 = read_trees(file.path(dir.hvrs["V7_L200"], file_list = list.files(path = dir.hvrs["V7_L200"], pattern = ".tre")))[1:N]
V8_L200 = read_trees(file.path(dir.hvrs["V8_L200"], file_list = list.files(path = dir.hvrs["V8_L200"], pattern = ".tre")))[1:N]

V1_L300 = read_trees(file.path(dir.hvrs["V1_L100"], file_list = list.files(path = dir.hvrs["V1_L100"], pattern = ".tre")))[1:N]
V2_L300 = read_trees(file.path(dir.hvrs["V2_L300"], file_list = list.files(path = dir.hvrs["V2_L300"], pattern = ".tre")))[1:N]
V3_L300 = read_trees(file.path(dir.hvrs["V3_L300"], file_list = list.files(path = dir.hvrs["V3_L300"], pattern = ".tre")))[1:N]
V4_L300 = read_trees(file.path(dir.hvrs["V4_L300"], file_list = list.files(path = dir.hvrs["V4_L300"], pattern = ".tre")))[1:N]
V5_L300 = read_trees(file.path(dir.hvrs["V5_L300"], file_list = list.files(path = dir.hvrs["V5_L300"], pattern = ".tre")))[1:N]
V6_L300 = read_trees(file.path(dir.hvrs["V6_L300"], file_list = list.files(path = dir.hvrs["V6_L300"], pattern = ".tre")))[1:N]
V7_L300 = read_trees(file.path(dir.hvrs["V7_L300"], file_list = list.files(path = dir.hvrs["V7_L300"], pattern = ".tre")))[1:N]

V1_L400 = read_trees(file.path(dir.hvrs["V1_L100"], file_list = list.files(path = dir.hvrs["V1_L100"], pattern = ".tre")))[1:N]
V2_L400 = read_trees(file.path(dir.hvrs["V2_L400"], file_list = list.files(path = dir.hvrs["V2_L400"], pattern = ".tre")))[1:N]
V3_L400 = read_trees(file.path(dir.hvrs["V3_L400"], file_list = list.files(path = dir.hvrs["V3_L400"], pattern = ".tre")))[1:N]
V4_L400 = read_trees(file.path(dir.hvrs["V4_L400"], file_list = list.files(path = dir.hvrs["V4_L400"], pattern = ".tre")))[1:N]
V5_L400 = read_trees(file.path(dir.hvrs["V5_L400"], file_list = list.files(path = dir.hvrs["V5_L400"], pattern = ".tre")))[1:N]
V6_L400 = read_trees(file.path(dir.hvrs["V6_L400"], file_list = list.files(path = dir.hvrs["V6_L400"], pattern = ".tre")))[1:N]

V1_L500 = read_trees(file.path(dir.hvrs["V1_L100"], file_list = list.files(path = dir.hvrs["V1_L100"], pattern = ".tre")))[1:N]
V2_L500 = read_trees(file.path(dir.hvrs["V2_L500"], file_list = list.files(path = dir.hvrs["V2_L500"], pattern = ".tre")))[1:N]
V3_L500 = read_trees(file.path(dir.hvrs["V3_L500"], file_list = list.files(path = dir.hvrs["V3_L500"], pattern = ".tre")))[1:N]
V4_L500 = read_trees(file.path(dir.hvrs["V4_L500"], file_list = list.files(path = dir.hvrs["V4_L500"], pattern = ".tre")))[1:N]
V5_L500 = read_trees(file.path(dir.hvrs["V5_L500"], file_list = list.files(path = dir.hvrs["V5_L500"], pattern = ".tre")))[1:N]

V1_L50_guides = read_trees(file.path(dir.hvrs["V1_L50_guides"], file_list = list.files(path = dir.hvrs["V1_L50_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V2_L50_guides = read_trees(file.path(dir.hvrs["V2_L50_guides"], file_list = list.files(path = dir.hvrs["V2_L50_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V3_L50_guides = read_trees(file.path(dir.hvrs["V3_L50_guides"], file_list = list.files(path = dir.hvrs["V3_L50_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V4_L50_guides = read_trees(file.path(dir.hvrs["V4_L50_guides"], file_list = list.files(path = dir.hvrs["V4_L50_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V5_L50_guides = read_trees(file.path(dir.hvrs["V5_L50_guides"], file_list = list.files(path = dir.hvrs["V5_L50_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V6_L50_guides = read_trees(file.path(dir.hvrs["V6_L50_guides"], file_list = list.files(path = dir.hvrs["V6_L50_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V7_L50_guides = read_trees(file.path(dir.hvrs["V7_L50_guides"], file_list = list.files(path = dir.hvrs["V7_L50_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V8_L50_guides = read_trees(file.path(dir.hvrs["V8_L50_guides"], file_list = list.files(path = dir.hvrs["V8_L50_guides"], pattern = "align.tre")), drop = TRUE)[1:N]

V1_L100_guides = read_trees(file.path(dir.hvrs["V1_L100_guides"], file_list = list.files(path = dir.hvrs["V1_L100_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V2_L100_guides = read_trees(file.path(dir.hvrs["V2_L100_guides"], file_list = list.files(path = dir.hvrs["V2_L100_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V3_L100_guides = read_trees(file.path(dir.hvrs["V3_L100_guides"], file_list = list.files(path = dir.hvrs["V3_L100_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V4_L100_guides = read_trees(file.path(dir.hvrs["V4_L100_guides"], file_list = list.files(path = dir.hvrs["V4_L100_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V5_L100_guides = read_trees(file.path(dir.hvrs["V5_L100_guides"], file_list = list.files(path = dir.hvrs["V5_L100_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V6_L100_guides = read_trees(file.path(dir.hvrs["V6_L100_guides"], file_list = list.files(path = dir.hvrs["V6_L100_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V7_L100_guides = read_trees(file.path(dir.hvrs["V7_L100_guides"], file_list = list.files(path = dir.hvrs["V7_L100_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V8_L100_guides = read_trees(file.path(dir.hvrs["V8_L100_guides"], file_list = list.files(path = dir.hvrs["V8_L100_guides"], pattern = "align.tre")), drop = TRUE)[1:N]

V1_L200_guides = read_trees(file.path(dir.hvrs["V1_L200_guides"], file_list = list.files(path = dir.hvrs["V1_L200_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V2_L200_guides = read_trees(file.path(dir.hvrs["V2_L200_guides"], file_list = list.files(path = dir.hvrs["V2_L200_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V3_L200_guides = read_trees(file.path(dir.hvrs["V3_L200_guides"], file_list = list.files(path = dir.hvrs["V3_L200_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V4_L200_guides = read_trees(file.path(dir.hvrs["V4_L200_guides"], file_list = list.files(path = dir.hvrs["V4_L200_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V5_L200_guides = read_trees(file.path(dir.hvrs["V5_L200_guides"], file_list = list.files(path = dir.hvrs["V5_L200_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V6_L200_guides = read_trees(file.path(dir.hvrs["V6_L200_guides"], file_list = list.files(path = dir.hvrs["V6_L200_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V7_L200_guides = read_trees(file.path(dir.hvrs["V7_L200_guides"], file_list = list.files(path = dir.hvrs["V7_L200_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V8_L200_guides = read_trees(file.path(dir.hvrs["V8_L200_guides"], file_list = list.files(path = dir.hvrs["V8_L200_guides"], pattern = "align.tre")), drop = TRUE)[1:N]

V1_L300_guides = read_trees(file.path(dir.hvrs["V1_L300_guides"], file_list = list.files(path = dir.hvrs["V1_L300_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V2_L300_guides = read_trees(file.path(dir.hvrs["V2_L300_guides"], file_list = list.files(path = dir.hvrs["V2_L300_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V3_L300_guides = read_trees(file.path(dir.hvrs["V3_L300_guides"], file_list = list.files(path = dir.hvrs["V3_L300_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V4_L300_guides = read_trees(file.path(dir.hvrs["V4_L300_guides"], file_list = list.files(path = dir.hvrs["V4_L300_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V5_L300_guides = read_trees(file.path(dir.hvrs["V5_L300_guides"], file_list = list.files(path = dir.hvrs["V5_L300_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V6_L300_guides = read_trees(file.path(dir.hvrs["V6_L300_guides"], file_list = list.files(path = dir.hvrs["V6_L300_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V7_L300_guides = read_trees(file.path(dir.hvrs["V7_L300_guides"], file_list = list.files(path = dir.hvrs["V7_L300_guides"], pattern = "align.tre")), drop = TRUE)[1:N]

V1_L400_guides = read_trees(file.path(dir.hvrs["V1_L400_guides"], file_list = list.files(path = dir.hvrs["V1_L400_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V2_L400_guides = read_trees(file.path(dir.hvrs["V2_L400_guides"], file_list = list.files(path = dir.hvrs["V2_L400_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V3_L400_guides = read_trees(file.path(dir.hvrs["V3_L400_guides"], file_list = list.files(path = dir.hvrs["V3_L400_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V4_L400_guides = read_trees(file.path(dir.hvrs["V4_L400_guides"], file_list = list.files(path = dir.hvrs["V4_L400_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V5_L400_guides = read_trees(file.path(dir.hvrs["V5_L400_guides"], file_list = list.files(path = dir.hvrs["V5_L400_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V6_L400_guides = read_trees(file.path(dir.hvrs["V6_L400_guides"], file_list = list.files(path = dir.hvrs["V6_L400_guides"], pattern = "align.tre")), drop = TRUE)[1:N]

V1_L500_guides = read_trees(file.path(dir.hvrs["V1_L500_guides"], file_list = list.files(path = dir.hvrs["V1_L500_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V2_L500_guides = read_trees(file.path(dir.hvrs["V2_L500_guides"], file_list = list.files(path = dir.hvrs["V2_L500_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V3_L500_guides = read_trees(file.path(dir.hvrs["V3_L500_guides"], file_list = list.files(path = dir.hvrs["V3_L500_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V4_L500_guides = read_trees(file.path(dir.hvrs["V4_L500_guides"], file_list = list.files(path = dir.hvrs["V4_L500_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V5_L500_guides = read_trees(file.path(dir.hvrs["V5_L500_guides"], file_list = list.files(path = dir.hvrs["V5_L500_guides"], pattern = "align.tre")), drop = TRUE)[1:N]

print("Reading in files...DONE")

print("Calculating distance of controls...")
VFull_control_path= calculate_path_dist_random_pairs(VFull, VFull_control, cores = CORES)
VFull_guides_path = calculate_path_dist_random_pairs(VFull, VFull_guides, cores = CORES)
print("Calculating distance of controls...DONE")

print("Calculating distance of V1...")
V1_L50_path = calculate_path_dist_random_pairs(VFull, V1_L50, cores = CORES)
V1_L100_path = calculate_path_dist_random_pairs(VFull, V1_L100, cores = CORES)
V1_L200_path = calculate_path_dist_random_pairs(VFull, V1_L200, cores = CORES)
V1_L300_path = calculate_path_dist_random_pairs(VFull, V1_L300, cores = CORES)
V1_L400_path = calculate_path_dist_random_pairs(VFull, V1_L400, cores = CORES)
V1_L500_path = calculate_path_dist_random_pairs(VFull, V1_L500, cores = CORES)

V1_L50_guides_path = calculate_path_dist_random_pairs(VFull, V1_L50_guides, cores = CORES)
V1_L100_guides_path = calculate_path_dist_random_pairs(VFull, V1_L100_guides, cores = CORES)
V1_L200_guides_path = calculate_path_dist_random_pairs(VFull, V1_L200_guides, cores = CORES)
V1_L300_guides_path = calculate_path_dist_random_pairs(VFull, V1_L300_guides, cores = CORES)
V1_L400_guides_path = calculate_path_dist_random_pairs(VFull, V1_L400_guides, cores = CORES)
V1_L500_guides_path = calculate_path_dist_random_pairs(VFull, V1_L500_guides, cores = CORES)

print("Calculating distance of V2...")
V2_L50_path = calculate_path_dist_random_pairs(VFull, V2_L50, cores = CORES)
V2_L100_path = calculate_path_dist_random_pairs(VFull, V2_L100, cores = CORES)
V2_L200_path = calculate_path_dist_random_pairs(VFull, V2_L200, cores = CORES)
V2_L300_path = calculate_path_dist_random_pairs(VFull, V2_L300, cores = CORES)
V2_L400_path = calculate_path_dist_random_pairs(VFull, V2_L400, cores = CORES)
V2_L500_path = calculate_path_dist_random_pairs(VFull, V2_L500, cores = CORES)

V2_L50_guides_path = calculate_path_dist_random_pairs(VFull, V2_L50_guides, cores = CORES)
V2_L100_guides_path = calculate_path_dist_random_pairs(VFull, V2_L100_guides, cores = CORES)
V2_L200_guides_path = calculate_path_dist_random_pairs(VFull, V2_L200_guides, cores = CORES)
V2_L300_guides_path = calculate_path_dist_random_pairs(VFull, V2_L300_guides, cores = CORES)
V2_L400_guides_path = calculate_path_dist_random_pairs(VFull, V2_L400_guides, cores = CORES)
V2_L500_guides_path = calculate_path_dist_random_pairs(VFull, V2_L500_guides, cores = CORES)

print("Calculating distance of V3...")
V3_L50_path = calculate_path_dist_random_pairs(VFull, V3_L50, cores = CORES)
V3_L100_path = calculate_path_dist_random_pairs(VFull, V3_L100, cores = CORES)
V3_L200_path = calculate_path_dist_random_pairs(VFull, V3_L200, cores = CORES)
V3_L300_path = calculate_path_dist_random_pairs(VFull, V3_L300, cores = CORES)
V3_L400_path = calculate_path_dist_random_pairs(VFull, V3_L400, cores = CORES)
V3_L500_path = calculate_path_dist_random_pairs(VFull, V3_L500, cores = CORES)

V3_L50_guides_path = calculate_path_dist_random_pairs(VFull, V3_L50_guides, cores = CORES)
V3_L100_guides_path = calculate_path_dist_random_pairs(VFull, V3_L100_guides, cores = CORES)
V3_L200_guides_path = calculate_path_dist_random_pairs(VFull, V3_L200_guides, cores = CORES)
V3_L300_guides_path = calculate_path_dist_random_pairs(VFull, V3_L300_guides, cores = CORES)
V3_L400_guides_path = calculate_path_dist_random_pairs(VFull, V3_L400_guides, cores = CORES)
V3_L500_guides_path = calculate_path_dist_random_pairs(VFull, V3_L500_guides, cores = CORES)

print("Calculating distance of V4...")
V4_L50_path = calculate_path_dist_random_pairs(VFull, V4_L50, cores = CORES)
V4_L100_path = calculate_path_dist_random_pairs(VFull, V4_L100, cores = CORES)
V4_L200_path = calculate_path_dist_random_pairs(VFull, V4_L200, cores = CORES)
V4_L300_path = calculate_path_dist_random_pairs(VFull, V4_L300, cores = CORES)
V4_L400_path = calculate_path_dist_random_pairs(VFull, V4_L400, cores = CORES)
V4_L500_path = calculate_path_dist_random_pairs(VFull, V4_L500, cores = CORES)

V4_L50_guides_path = calculate_path_dist_random_pairs(VFull, V4_L50_guides, cores = CORES)
V4_L100_guides_path = calculate_path_dist_random_pairs(VFull, V4_L100_guides, cores = CORES)
V4_L200_guides_path = calculate_path_dist_random_pairs(VFull, V4_L200_guides, cores = CORES)
V4_L300_guides_path = calculate_path_dist_random_pairs(VFull, V4_L300_guides, cores = CORES)
V4_L400_guides_path = calculate_path_dist_random_pairs(VFull, V4_L400_guides, cores = CORES)
V4_L500_guides_path = calculate_path_dist_random_pairs(VFull, V4_L500_guides, cores = CORES)

print("Calculating distance of V5...")
V5_L50_path = calculate_path_dist_random_pairs(VFull, V5_L50, cores = CORES)
V5_L100_path = calculate_path_dist_random_pairs(VFull, V5_L100, cores = CORES)
V5_L200_path = calculate_path_dist_random_pairs(VFull, V5_L200, cores = CORES)
V5_L300_path = calculate_path_dist_random_pairs(VFull, V5_L300, cores = CORES)
V5_L400_path = calculate_path_dist_random_pairs(VFull, V5_L400, cores = CORES)
V5_L500_path = calculate_path_dist_random_pairs(VFull, V5_L500, cores = CORES)

V5_L50_guides_path = calculate_path_dist_random_pairs(VFull, V5_L50_guides, cores = CORES)
V5_L100_guides_path = calculate_path_dist_random_pairs(VFull, V5_L100_guides, cores = CORES)
V5_L200_guides_path = calculate_path_dist_random_pairs(VFull, V5_L200_guides, cores = CORES)
V5_L300_guides_path = calculate_path_dist_random_pairs(VFull, V5_L300_guides, cores = CORES)
V5_L400_guides_path = calculate_path_dist_random_pairs(VFull, V5_L400_guides, cores = CORES)
V5_L500_guides_path = calculate_path_dist_random_pairs(VFull, V5_L500_guides, cores = CORES)

print("Calculating distance of V6...")
V6_L50_path = calculate_path_dist_random_pairs(VFull, V6_L50, cores = CORES)
V6_L100_path = calculate_path_dist_random_pairs(VFull, V6_L100, cores = CORES)
V6_L200_path = calculate_path_dist_random_pairs(VFull, V6_L200, cores = CORES)
V6_L300_path = calculate_path_dist_random_pairs(VFull, V6_L300, cores = CORES)
V6_L400_path = calculate_path_dist_random_pairs(VFull, V6_L400, cores = CORES)

V6_L50_guides_path = calculate_path_dist_random_pairs(VFull, V6_L50_guides, cores = CORES)
V6_L100_guides_path = calculate_path_dist_random_pairs(VFull, V6_L100_guides, cores = CORES)
V6_L200_guides_path = calculate_path_dist_random_pairs(VFull, V6_L200_guides, cores = CORES)
V6_L300_guides_path = calculate_path_dist_random_pairs(VFull, V6_L300_guides, cores = CORES)
V6_L400_guides_path = calculate_path_dist_random_pairs(VFull, V6_L400_guides, cores = CORES)

print("Calculating distance of V7...")
V7_L50_path = calculate_path_dist_random_pairs(VFull, V7_L50, cores = CORES)
V7_L100_path = calculate_path_dist_random_pairs(VFull, V7_L100, cores = CORES)
V7_L200_path = calculate_path_dist_random_pairs(VFull, V7_L200, cores = CORES)
V7_L300_path = calculate_path_dist_random_pairs(VFull, V7_L300, cores = CORES)

V7_L50_guides_path = calculate_path_dist_random_pairs(VFull, V7_L50_guides, cores = CORES)
V7_L100_guides_path = calculate_path_dist_random_pairs(VFull, V7_L100_guides, cores = CORES)
V7_L200_guides_path = calculate_path_dist_random_pairs(VFull, V7_L200_guides, cores = CORES)
V7_L300_guides_path = calculate_path_dist_random_pairs(VFull, V7_L300_guides, cores = CORES)

print("Calculating distance of V8...")
V8_L50_path = calculate_path_dist_random_pairs(VFull, V8_L50, cores = CORES)
V8_L100_path = calculate_path_dist_random_pairs(VFull, V8_L100, cores = CORES)
V8_L200_path = calculate_path_dist_random_pairs(VFull, V8_L200, cores = CORES)

V8_L50_guides_path = calculate_path_dist_random_pairs(VFull, V8_L50_guides, cores = CORES)
V8_L100_guides_path = calculate_path_dist_random_pairs(VFull, V8_L100_guides, cores = CORES)
V8_L200_guides_path = calculate_path_dist_random_pairs(VFull, V8_L200_guides, cores = CORES)

setwd(dir.save)
save.image(file = "sim_analysis_bp_path_save.RData", compress = TRUE)

bp_path = 
  data.frame("VFull_control_path" = VFull_control_path,
             "VFull_guides_path" = VFull_guides_path,
             "V1_L50_path" = V1_L50_path,
             "V1_L100_path" = V1_L100_path,
             "V1_L200_path" = V1_L200_path,
             "V1_L300_path" = V1_L300_path,
             "V1_L400_path" = V1_L400_path,
             "V1_L500_path" = V1_L500_path,
             "V1_L50_guides_path" = V1_L50_guides_path,
             "V1_L100_guides_path" = V1_L100_guides_path,
             "V1_L200_guides_path" = V1_L200_guides_path,
             "V1_L300_guides_path" = V1_L300_guides_path,
             "V1_L400_guides_path" = V1_L400_guides_path,
             "V1_L500_guides_path" = V1_L500_guides_path,
             "V2_L50_path" = V2_L50_path,
             "V2_L100_path" = V2_L100_path,
             "V2_L200_path" = V2_L200_path,
             "V2_L300_path" = V2_L300_path,
             "V2_L400_path" = V2_L400_path,
             "V2_L500_path" = V2_L500_path,
             "V2_L50_guides_path" = V2_L50_guides_path,
             "V2_L100_guides_path" = V2_L100_guides_path,
             "V2_L200_guides_path" = V2_L200_guides_path,
             "V2_L300_guides_path" = V2_L300_guides_path,
             "V2_L400_guides_path" = V2_L400_guides_path,
             "V2_L500_guides_path" = V2_L500_guides_path,
             "V3_L50_path" = V3_L50_path,
             "V3_L100_path" = V3_L100_path,
             "V3_L200_path" = V3_L200_path,
             "V3_L300_path" = V3_L300_path,
             "V3_L400_path" = V3_L400_path,
             "V3_L500_path" = V3_L500_path,
             "V3_L50_guides_path" = V3_L50_guides_path,
             "V3_L100_guides_path" = V3_L100_guides_path,
             "V3_L200_guides_path" = V3_L200_guides_path,
             "V3_L300_guides_path" = V3_L300_guides_path,
             "V3_L400_guides_path" = V3_L400_guides_path,
             "V3_L500_guides_path" = V3_L500_guides_path,
             "V4_L50_path" = V4_L50_path,
             "V4_L100_path" = V4_L100_path,
             "V4_L200_path" = V4_L200_path,
             "V4_L300_path" = V4_L300_path,
             "V4_L400_path" = V4_L400_path,
             "V4_L500_path" = V4_L500_path,
             "V4_L50_guides_path" = V4_L50_guides_path,
             "V4_L100_guides_path" = V4_L100_guides_path,
             "V4_L200_guides_path" = V4_L200_guides_path,
             "V4_L300_guides_path" = V4_L300_guides_path,
             "V4_L400_guides_path" = V4_L400_guides_path,
             "V4_L500_guides_path" = V4_L500_guides_path,
             "V5_L50_path" = V5_L50_path,
             "V5_L100_path" = V5_L100_path,
             "V5_L200_path" = V5_L200_path,
             "V5_L300_path" = V5_L300_path,
             "V5_L400_path" = V5_L400_path,
             "V5_L500_path" = V5_L500_path,
             "V5_L50_guides_path" = V5_L50_guides_path,
             "V5_L100_guides_path" = V5_L100_guides_path,
             "V5_L200_guides_path" = V5_L200_guides_path,
             "V5_L300_guides_path" = V5_L300_guides_path,
             "V5_L400_guides_path" = V5_L400_guides_path,
             "V5_L500_guides_path" = V5_L500_guides_path,
             "V6_L50_path" = V6_L50_path,
             "V6_L100_path" = V6_L100_path,
             "V6_L200_path" = V6_L200_path,
             "V6_L300_path" = V6_L300_path,
             "V6_L400_path" = V6_L400_path,
             "V6_L50_guides_path" = V6_L50_guides_path,
             "V6_L100_guides_path" = V6_L100_guides_path,
             "V6_L200_guides_path" = V6_L200_guides_path,
             "V6_L300_guides_path" = V6_L300_guides_path,
             "V6_L400_guides_path" = V6_L400_guides_path,
             "V7_L50_path" = V7_L50_path,
             "V7_L100_path" = V7_L100_path,
             "V7_L200_path" = V7_L200_path,
             "V7_L300_path" = V7_L300_path,
             "V7_L50_guides_path" = V7_L50_guides_path,
             "V7_L100_guides_path" = V7_L100_guides_path,
             "V7_L200_guides_path" = V7_L200_guides_path,
             "V7_L300_guides_path" = V7_L300_guides_path,
             "V8_L50_path" = V8_L50_path,
             "V8_L100_path" = V8_L100_path,
             "V8_L200_path" = V8_L200_path,
             "V8_L50_guides_path" = V8_L50_guides_path,
             "V8_L100_guides_path" = V8_L100_guides_path,
             "V8_L200_guides_path" = V8_L200_guides_path)
setwd(dir.sim.out)
saveRDS(file = "bp_path_random_pairs.rds", object = bp_path)
