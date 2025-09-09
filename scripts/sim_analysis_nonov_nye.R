
################################################################################
# CORES
################################################################################
CORES = 100 # Number of cores to use for the simulation
SAVE_FILENAME = "nonov_nye_random_pairs.rds"
N = 100 # trees to look at per simulation

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
dir.sim.out = file.path(HOME, "2025_hvr_guide_phylogenetic_integration/out/sim_nonoverlapping/")

setwd(HOME_TMP)
dir.simulation = file.path(HOME_TMP, "length_simulations/")
dir.mix = file.path(HOME_TMP, "mix_experiment/")

setwd(dir.scripts)
source("sim_analysis_functions_V2.R")

################################################################################
# SET DIRECTORIES
################################################################################
# Check that all the tree files are present. 
directories <- c(list.dirs(dir.mix, full.names = TRUE, recursive = FALSE),
                 list.dirs(dir.simulation, full.names = TRUE, recursive = FALSE))
names(directories) = basename(directories)
basename(directories)

################################################################################
# 0A: READ IN FILES
################################################################################

setwd(dir.mix)

dir.mix = directories[c("V1V2", "V1V2_guides", 
                        "V3V4", "V3V4_guides", 
                        "V5V7", "V5V7_guides", 
                        "V1V2_V3V4_V5V7_mix", "V1V2_V3V4_V5V7_mix_guides")]
dir.mix

setwd(dir.simulation)
dir.control = directories[c("VFull", "VFull_control", "VFull_guides")]
dir.hvrs = c(dir.mix, dir.control)

setwd(dir.mix)
print("Reading in trees...")
VFull = read_trees(file.path(dir.hvrs["VFull"], file_list = list.files(path = dir.hvrs["VFull"], pattern = ".tre")))[1:N]
VFull_guides = read_trees(file.path(dir.hvrs["VFull_guides"], file_list = list.files(path = dir.hvrs["VFull_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
VFull_control = read_trees(file.path(dir.hvrs["VFull_control"], file_list = list.files(path = dir.hvrs["VFull_control"], pattern = "align.tre")))[1:N]

V1V2 = read_trees(file.path(dir.hvrs["V1V2"], file_list = list.files(path = dir.hvrs["V1V2"], pattern = ".tre")))[1:N]
V3V4 = read_trees(file.path(dir.hvrs["V3V4"], file_list = list.files(path = dir.hvrs["V3V4"], pattern = ".tre")))[1:N]
V5V7 = read_trees(file.path(dir.hvrs["V5V7"], file_list = list.files(path = dir.hvrs["V5V7"], pattern = ".tre")))[1:N]
V1V2_V3V4_V5V7_mix = read_trees(file.path(dir.hvrs["V1V2_V3V4_V5V7_mix"], file_list = list.files(path = dir.hvrs["V1V2_V3V4_V5V7_mix"], pattern = ".tre")))[1:N]

V1V2_guides = read_trees(file.path(dir.hvrs["V1V2_guides"], file_list = list.files(path = dir.hvrs["V1V2_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V3V4_guides = read_trees(file.path(dir.hvrs["V3V4_guides"], file_list = list.files(path = dir.hvrs["V3V4_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V5V7_guides = read_trees(file.path(dir.hvrs["V5V7_guides"], file_list = list.files(path = dir.hvrs["V5V7_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
V1V2_V3V4_V5V7_mix_guides = read_trees(file.path(dir.hvrs["V1V2_V3V4_V5V7_mix_guides"], file_list = list.files(path = dir.hvrs["V1V2_V3V4_V5V7_mix_guides"], pattern = "align.tre")), drop = TRUE)[1:N]
print("Reading in trees...DONE")

################################################################################
# 
################################################################################
VFull_control = calculate_nye_distance_random_pairs(VFull, VFull_control, cores = CORES)
VFull_guides = calculate_nye_distance_random_pairs(VFull, VFull_guides, cores = CORES)

V1V2 = calculate_nye_distance_random_pairs(VFull, V1V2, cores = CORES)
V3V4 = calculate_nye_distance_random_pairs(VFull, V3V4, cores = CORES)
V5V7 = calculate_nye_distance_random_pairs(VFull, V5V7, cores = CORES)

V1V2_guides = calculate_nye_distance_random_pairs(VFull, V1V2_guides, cores = CORES)
V3V4_guides = calculate_nye_distance_random_pairs(VFull, V3V4_guides, cores = CORES)
V5V7_guides = calculate_nye_distance_random_pairs(VFull, V5V7_guides, cores = CORES)

V1V2_V3V4_V5V7_mix = calculate_nye_distance_random_pairs(VFull, V1V2_V3V4_V5V7_mix, cores = CORES)
V1V2_V3V4_V5V7_mix_guides = calculate_nye_distance_random_pairs(VFull, V1V2_V3V4_V5V7_mix_guides, cores = CORES)

nonov_rf = data.frame("VFull" = VFull_control,
                      "VFull_guides" = VFull_guides,
                      "V1V2" = V1V2,
                      "V3V4" = V3V4,
                      "V5V7" = V5V7,
                      "V1V2_guides" = V1V2_guides,
                      "V3V4_guides" = V3V4_guides,
                      "V5V7_guides" = V5V7_guides,
                      "V1V2_V3V4_V5V7_mix" = V1V2_V3V4_V5V7_mix,
                      "V1V2_V3V4_V5V7_mix_guides" = V1V2_V3V4_V5V7_mix_guides)

setwd(dir.sim.out)
saveRDS(file = SAVE_FILENAME, object = nonov_rf)