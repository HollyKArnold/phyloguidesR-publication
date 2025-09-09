################################################################################
# sim_analysis_ov_nanoclusters_random_pairs_3.R

################################################################################

################################################################################
# SET CONSTANTS
################################################################################
CORES = 100
SAVE_FILE_NAME = "ov_nanoclusters_random_pairs_3.rds"
SUBTREE_SIZE = 3
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
# INPUTS
################################################################################

# Set Home
HOME = "/nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration"
HOME_TMP = "/nfs3/Sharpton_Lab/tmp/projects/arnoldh/2025_hvr_guide_phylogenetic_integration/"

# Make directories
setwd(HOME)
dir.scripts = file.path(HOME, "2025_hvr_guide_phylogenetic_integration/scripts/")
dir.ref = file.path(HOME, "reference_data/")
dir.out = file.path(HOME, "2025_hvr_guide_phylogenetic_integration/data_visulization/")
dir.sim.out = file.path(HOME, "2025_hvr_guide_phylogenetic_integration/out/sim_overlapping/")

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
# 0: READ IN FILES
################################################################################
N = 100 # The number of trees to look at per simulation
setwd(dir.mix)
dir.mix = directories[c("V3V4", "V3V4_guides", "V4", "V4_guides", "V4_V3V4_mix", "V4_V3V4_mix_guides")]
dir.mix

setwd(dir.simulation)
dir.control = directories[c("VFull", "VFull_control", "VFull_guides")]
dir.control

dir.hvrs = c(dir.mix, dir.control)

setwd(dir.mix)
VFull = read_trees(file.path(dir.hvrs["VFull"], file_list = list.files(path = dir.hvrs["VFull"], pattern = ".tre")))[1:N]
VFull_guides = read_trees(file.path(dir.hvrs["VFull_guides"], file_list = list.files(path = dir.hvrs["VFull_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
VFull_control = read_trees(file.path(dir.hvrs["VFull_control"], file_list = list.files(path = dir.hvrs["VFull_control"], pattern = "align.tre")))[1:N]

V4 = read_trees(file.path(dir.hvrs["V4"], file_list = list.files(path = dir.hvrs["V4"], pattern = ".tre")))[1:N]
V3V4 = read_trees(file.path(dir.hvrs["V3V4"], file_list = list.files(path = dir.hvrs["V3V4"], pattern = ".tre")))[1:N]
V4_V3V4_mix = read_trees(file.path(dir.hvrs["V4_V3V4_mix"], file_list = list.files(path = dir.hvrs["V4_V3V4_mix"], pattern = ".tre")))[1:N]

V4_guides = read_trees(file.path(dir.hvrs["V4_guides"], file_list = list.files(path = dir.hvrs["V4_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V3V4_guides = read_trees(file.path(dir.hvrs["V3V4_guides"], file_list = list.files(path = dir.hvrs["V3V4_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V4_V3V4_mix_guides = read_trees(file.path(dir.hvrs["V4_V3V4_mix_guides"], file_list = list.files(path = dir.hvrs["V4_V3V4_mix_guides"], pattern = "align.tre")), drop = TRUE)[1:N]

################################################################################
# 
################################################################################
VFull_control = calculate_shared_subtrees_dist_random_pairs(VFull, VFull_control, cores = CORES, n = SUBTREE_SIZE)
VFull_guides = calculate_shared_subtrees_dist_random_pairs(VFull, VFull_guides, cores = CORES, n = SUBTREE_SIZE)
V4_rf = calculate_shared_subtrees_dist_random_pairs(VFull, V4, cores = CORES, n = SUBTREE_SIZE)
V3V4_rf = calculate_shared_subtrees_dist_random_pairs(VFull, V3V4, cores = CORES, n = SUBTREE_SIZE)
V4_guides_rf = calculate_shared_subtrees_dist_random_pairs(VFull, V4_guides, cores = CORES, n = SUBTREE_SIZE)
V3V4_guides_rf = calculate_shared_subtrees_dist_random_pairs(VFull, V3V4_guides, cores = CORES, n = SUBTREE_SIZE)
V4_V3V4_mix_rf = calculate_shared_subtrees_dist_random_pairs(VFull, V4_V3V4_mix, cores = CORES, n = SUBTREE_SIZE)
V4_V3V4_mix_guides_rf = calculate_shared_subtrees_dist_random_pairs(VFull, V4_V3V4_mix_guides, cores = CORES, n = SUBTREE_SIZE)


dist = data.frame("VFull" = VFull_control,
                  "VFull_guides" = VFull_guides,
                  "V4_rf" = V4_rf,
                  "V3V4_rf" = V3V4_rf,
                  "V4_guides_rf" = V4_guides_rf,
                  "V3V4_guides_rf" = V3V4_guides_rf,
                  "V4_V3V4_mix_rf" = V4_V3V4_mix_rf,
                  "V4_V3V4_mix_guides_rf" = V4_V3V4_mix_guides_rf)

setwd(dir.sim.out)
saveRDS(file = SAVE_FILE_NAME, object = dist)
