
################################################################################
# SET CONSTANTS
################################################################################
N = 100 # Number of trees to consider
CORES = 30 # Number of cores
SUBTREE_SIZE = 2
SAVE_FILE_NAME = paste0(c("nanostructure_nd_", SUBTREE_SIZE, "_random_pairs.rds"), sep = "", collapse = "")

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
dir.sim.out = file.path(HOME, "2025_hvr_guide_phylogenetic_integration/out/sim_nd_mix/")
setwd(HOME_TMP)
dir.simulation = file.path(HOME_TMP, "length_simulations/")
dir.mix = file.path(HOME_TMP, "mix_experiment/")

# Read functions
setwd(dir.scripts)
source("sim_analysis_functions_V2.R")

################################################################################
# INPUT
################################################################################
# Check that all the tree files are present. 
directories <- c(list.dirs(dir.mix, full.names = TRUE, recursive = FALSE),
                 list.dirs(dir.simulation, full.names = TRUE, recursive = FALSE))
names(directories) = basename(directories)
basename(directories)

################################################################################
# 0A: READ IN FILES
################################################################################

# Read in files
setwd(dir.mix)
dir.mix = directories[c("V1V2", "V1V2_guides", 
                        "V1V3", "V1V3_guides", 
                        "V3V4", "V3V4_guides", 
                        "V3V5", "V3V5_guides",
                        "V4", "V4_guides",
                        "V4V5", "V4V5_guides", 
                        "nd_mix", "nd_mix_guides")]

setwd(dir.simulation)
dir.control = directories[c("VFull", "VFull_control", "VFull_guides")]
dir.hvrs = c(dir.mix, dir.control)

# Read in trees
setwd(dir.mix)
VFull = read_trees(file.path(dir.hvrs["VFull"], file_list = list.files(path = dir.hvrs["VFull"], pattern = ".tre")))[1:N]
VFull_guides = read_trees(file.path(dir.hvrs["VFull_guides"], file_list = list.files(path = dir.hvrs["VFull_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
VFull_control = read_trees(file.path(dir.hvrs["VFull_control"], file_list = list.files(path = dir.hvrs["VFull_control"], pattern = "align.tre")))[1:N]

V1V2 = read_trees(file.path(dir.hvrs["V1V2"], file_list = list.files(path = dir.hvrs["V1V2"], pattern = ".tre")))[1:N]
V1V3 = read_trees(file.path(dir.hvrs["V1V3"], file_list = list.files(path = dir.hvrs["V1V3"], pattern = ".tre")))[1:N]
V3V4 = read_trees(file.path(dir.hvrs["V3V4"], file_list = list.files(path = dir.hvrs["V3V4"], pattern = ".tre")))[1:N]
V3V5 = read_trees(file.path(dir.hvrs["V3V5"], file_list = list.files(path = dir.hvrs["V3V5"], pattern = ".tre")))[1:N]
V4 = read_trees(file.path(dir.hvrs["V4"], file_list = list.files(path = dir.hvrs["V4"], pattern = ".tre")))[1:N]
V4V5 = read_trees(file.path(dir.hvrs["V4V5"], file_list = list.files(path = dir.hvrs["V4V5"], pattern = ".tre")))[1:N]
nd_mix = read_trees(file.path(dir.hvrs["nd_mix"], file_list = list.files(path = dir.hvrs["nd_mix"], pattern = ".tre")))[1:N]

V1V2_guides = read_trees(file.path(dir.hvrs["V1V2_guides"], file_list = list.files(path = dir.hvrs["V1V2_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V1V3_guides = read_trees(file.path(dir.hvrs["V1V3_guides"], file_list = list.files(path = dir.hvrs["V1V3_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V3V4_guides = read_trees(file.path(dir.hvrs["V3V4_guides"], file_list = list.files(path = dir.hvrs["V3V4_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V3V5_guides = read_trees(file.path(dir.hvrs["V3V5_guides"], file_list = list.files(path = dir.hvrs["V3V5_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V4_guides = read_trees(file.path(dir.hvrs["V4_guides"], file_list = list.files(path = dir.hvrs["V4_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V4V5_guides = read_trees(file.path(dir.hvrs["V4V5_guides"], file_list = list.files(path = dir.hvrs["V4V5_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
nd_mix_guides = read_trees(file.path(dir.hvrs["nd_mix_guides"], file_list = list.files(path = dir.hvrs["nd_mix_guides"], pattern = "align.tre")), drop = TRUE)[1:N]


################################################################################
# 
################################################################################

VFull_control = calculate_shared_subtrees_dist_random_pairs(VFull, VFull_control, cores = CORES, n = SUBTREE_SIZE)
VFull_guides = calculate_shared_subtrees_dist_random_pairs(VFull, VFull_guides, cores = CORES, n = SUBTREE_SIZE)

V1V2 = calculate_shared_subtrees_dist_random_pairs(VFull, V1V2, cores = CORES, n = SUBTREE_SIZE)
V1V3 = calculate_shared_subtrees_dist_random_pairs(VFull, V1V3, cores = CORES, n = SUBTREE_SIZE)
V3V4 = calculate_shared_subtrees_dist_random_pairs(VFull, V3V4, cores = CORES, n = SUBTREE_SIZE)
V3V5 = calculate_shared_subtrees_dist_random_pairs(VFull, V3V5, cores = CORES, n = SUBTREE_SIZE)
V4 = calculate_shared_subtrees_dist_random_pairs(VFull, V4, cores = CORES, n = SUBTREE_SIZE)
V4V5 = calculate_shared_subtrees_dist_random_pairs(VFull, V4V5, cores = CORES, n = SUBTREE_SIZE)
nd_mix = calculate_shared_subtrees_dist_random_pairs(VFull, nd_mix, cores = CORES, n = SUBTREE_SIZE)

V1V2_guides = calculate_shared_subtrees_dist_random_pairs(VFull, V1V2_guides, cores = CORES, n = SUBTREE_SIZE)
V1V3_guides = calculate_shared_subtrees_dist_random_pairs(VFull, V1V3_guides, cores = CORES, n = SUBTREE_SIZE)
V3V4_guides = calculate_shared_subtrees_dist_random_pairs(VFull, V3V4_guides, cores = CORES, n = SUBTREE_SIZE)
V3V5_guides = calculate_shared_subtrees_dist_random_pairs(VFull, V3V5_guides, cores = CORES, n = SUBTREE_SIZE)
V4_guides = calculate_shared_subtrees_dist_random_pairs(VFull, V4_guides, cores = CORES, n = SUBTREE_SIZE)
V4V5_guides = calculate_shared_subtrees_dist_random_pairs(VFull, V4V5_guides, cores = CORES, n = SUBTREE_SIZE)
nd_mix_guides = calculate_shared_subtrees_dist_random_pairs(VFull, nd_mix_guides, cores = CORES, n = SUBTREE_SIZE)

dist = data.frame("VFull" = VFull_control,
                  "VFull_guides" = VFull_guides,
                  
                  "V1V2" = V1V2,
                  "V1V3" = V1V3,
                  "V3V4" = V3V4,
                  "V3V5" = V3V5,
                  "V4" = V4,
                  "V4V5" = V4V5,
                  "nd_mix" = nd_mix,
                  
                  "V1V2_guides" = V1V2_guides,
                  "V1V3_guides" = V1V3_guides,
                  "V3V4_guides" = V3V4_guides,
                  "V3V5_guides" = V3V5_guides,
                  "V4_guides" = V4_guides,
                  "V4V5_guides" = V4V5_guides,
                  "nd_mix_guides" = nd_mix_guides)

setwd(dir.sim.out)
saveRDS(file = SAVE_FILE_NAME, object = dist)

