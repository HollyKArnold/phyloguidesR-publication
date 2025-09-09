################################################################################
#AUTHOR: ARNOLD
#DAY: Jan 24th, 2025
#SCRIPT: sim_analysis_master_V5.R
#DESCRIPTION: 
# SIM 1) Do guides improve tree accuracy across all VRs? 
# SIM 2) Do guides improve tree accuracy across all VRs and different lengths? 
# SIM 3) Do guides enable integration across overlapping VRs? 
# SIM 4) Do guides enable integration across non-overlapping VRs? 
# SIM 5) Do guides enable mixing a diverse set of VRs? 
# EXP) Do guides enable the same results when using true experimental long-read and short-read data?  Addition of experimental data to validate the findings here. 
################################################################################

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
dir.demo = file.path(HOME_TMP, "experimental_demo/")
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

################################################################################
# 0. SIM PREPREP: Check that all files have data in them
################################################################################

# Sanity check that there are no empty files
# size = vector()
# for(i in 1:length(directories)){
#   cur_dir = directories[i]
#   files = list.files(cur_dir, pattern = ".tre")
#   #print(paste0(c("Checking dir: ", basename(cur_dir)), sep = "", collapse = ""))
#   if(length(files) != 100){
#     message = paste0(c("Missing some files in dir: ", cur_dir), sep = "", collapse = "")
#     print(message)
#   }
#   for(j in 1:length(files)){
#     cur_file = files[j]
#     size = c(size, file.size(file.path(cur_dir, cur_file)))
#     if(file.size(file.path(cur_dir, cur_file)) == 677){
#       print(cur_file)
#     }
#   }
# }
# if(min(size) > 1000){
#   print("PASS: All files have something in them")
# }else{
#   print("FAIL: There are some smaller files that could be concerning.")
# }

################################################################################
# 0. SIM PREPREP: Visualize length distributions for each HVR in a barchart.
################################################################################
# setwd(dir.ref)
# hvrs = as_tibble(read.fasta("seqs_V1_LHVR.fasta")) %>%
#   rename("V1" = seq.text) %>%
#   left_join(as_tibble(read.fasta("seqs_V2_LHVR.fasta")), join_by(seq.name)) %>%
#   rename("V2" = seq.text) %>%
#   left_join(as_tibble(read.fasta("seqs_V3_LHVR.fasta")), join_by(seq.name)) %>%
#   rename("V3" = seq.text) %>%
#   left_join(as_tibble(read.fasta("seqs_V4_LHVR.fasta")), join_by(seq.name)) %>%
#   rename("V4" = seq.text) %>%
#   left_join(as_tibble(read.fasta("seqs_V5_LHVR.fasta")), join_by(seq.name)) %>%
#   rename("V5" = seq.text) %>%
#   left_join(as_tibble(read.fasta("seqs_V6_LHVR.fasta")), join_by(seq.name)) %>%
#   rename("V6" = seq.text) %>% 
#   left_join(as_tibble(read.fasta("seqs_V7_LHVR.fasta")), join_by(seq.name)) %>%
#   rename("V7" = seq.text) %>%
#   left_join(as_tibble(read.fasta("seqs_V8_LHVR.fasta")), join_by(seq.name)) %>%
#   rename("V8" = seq.text) %>%
#   left_join(as_tibble(read.fasta("seqs_V9_LHVR.fasta")), join_by(seq.name)) %>%
#   rename("V9" = seq.text) 
# 
# hvrs = 
#   hvrs %>%
#   mutate(V1 = stringr::str_length(V1)) %>%
#   mutate(V2 = stringr::str_length(V2)) %>%
#   mutate(V3 = stringr::str_length(V3)) %>%
#   mutate(V4 = stringr::str_length(V4)) %>%
#   mutate(V5 = stringr::str_length(V5)) %>%
#   mutate(V6 = stringr::str_length(V6)) %>%
#   mutate(V7 = stringr::str_length(V7)) %>%
#   mutate(V8 = stringr::str_length(V8)) %>%
#   mutate(V9 = stringr::str_length(V9)) 
# 
# hvrs = melt(hvrs, variable.name = "HVR", value.name = "Length")
# head(hvrs)
# 
# p = ggplot(hvrs, aes(x = HVR, y = Length, fill = HVR)) +
#   geom_boxplot() +
#   theme_minimal() +
#   labs(title = "Box Plot of Length Across HVR Regions",
#        x = "HVR",
#        y = "Length") +
#   theme(legend.position = "none")
# 
# p
# setwd(dir.out)
#saveRDS(file = "hvr_by_length.rds", hvrs, compress = TRUE) # START PROOF


################################################################################
# SIM 1: How does VR impact ability to topologically reconstruct trees?
################################################################################

#########################
# 1A. RF Distance - Random pairs
#########################
# 
# CORES = 50
# v1_rf = calculate_rf_distance_random_pairs(VFull, V1, cores = CORES)
# v2_rf = calculate_rf_distance_random_pairs(VFull, V2, cores = CORES)
# v3_rf = calculate_rf_distance_random_pairs(VFull, V3, cores = CORES)
# v4_rf = calculate_rf_distance_random_pairs(VFull, V4, cores = CORES)
# v5_rf = calculate_rf_distance_random_pairs(VFull, V5, cores = CORES)
# v6_rf = calculate_rf_distance_random_pairs(VFull, V6, cores = CORES)
# v7_rf = calculate_rf_distance_random_pairs(VFull, V7, cores = CORES)
# v8_rf = calculate_rf_distance_random_pairs(VFull, V8, cores = CORES)
# v9_rf = calculate_rf_distance_random_pairs(VFull, V9, cores = CORES)
# 
# v1_guides_rf = calculate_rf_distance_random_pairs(VFull, V1_guides, cores = CORES)
# v2_guides_rf = calculate_rf_distance_random_pairs(VFull, V2_guides, cores = CORES)
# v3_guides_rf = calculate_rf_distance_random_pairs(VFull, V3_guides, cores = CORES)
# v4_guides_rf = calculate_rf_distance_random_pairs(VFull, V4_guides, cores = CORES)
# v5_guides_rf = calculate_rf_distance_random_pairs(VFull, V5_guides, cores = CORES)
# v6_guides_rf = calculate_rf_distance_random_pairs(VFull, V6_guides, cores = CORES)
# v7_guides_rf = calculate_rf_distance_random_pairs(VFull, V7_guides, cores = CORES)
# v8_guides_rf = calculate_rf_distance_random_pairs(VFull, V8_guides, cores = CORES)
# v9_guides_rf = calculate_rf_distance_random_pairs(VFull, V9_guides, cores = CORES)
# 
# v_control_rf = calculate_rf_distance_random_pairs(VFull, VFull_control, cores = CORES)
# v_guides_rf = calculate_rf_distance_random_pairs(VFull, VFull_guides, cores = CORES)
# 
# # # Make Data Frame Of normalized RF
# hvr_rf_random_pairs =
#   data.frame("V1" = v1_rf,
#              "V2" = v2_rf,
#              "V3" = v3_rf,
#              "V4" = v4_rf,
#              "V5" = v5_rf,
#              "V6" = v6_rf,
#              "V7" = v7_rf,
#              "V8" = v8_rf,
#              "V9" = v9_rf,
#              "V1_guides" = v1_guides_rf,
#              "V2_guides" = v2_guides_rf,
#              "V3_guides" = v3_guides_rf,
#              "V4_guides" = v4_guides_rf,
#              "V5_guides" = v5_guides_rf,
#              "V6_guides" = v6_guides_rf,
#              "V7_guides" = v7_guides_rf,
#              "V8_guides" = v8_guides_rf,
#              "V9_guides" = v9_guides_rf,
#              "VFull" = v_control_rf,
#              "VFull_guides" = v_guides_rf)
# 

setwd(dir.out)
# saveRDS(hvr_rf_random_pairs, file = "hvr_rf_random_pairs.rds", compress = TRUE)
hvr_rf_random_pairs = readRDS("hvr_rf_random_pairs.rds")


#########################
# 1B. Quartet Distance - For random pairwise comparisons
#########################
# VFull_paths = file.path(dir.hvrs["VFull"], file_list = list.files(path = dir.hvrs["VFull"], pattern = ".tre"))
# VFull_control_paths = file.path(dir.hvrs["VFull_control"], file_list = list.files(path = dir.hvrs["VFull_control"], pattern = ".tre"))
# V1_paths = file.path(dir.hvrs["V1_HVR"], file_list = list.files(path = dir.hvrs["V1_HVR"], pattern = ".tre"))
# V2_paths = file.path(dir.hvrs["V2_HVR"], file_list = list.files(path = dir.hvrs["V2_HVR"], pattern = ".tre"))
# V3_paths = file.path(dir.hvrs["V3_HVR"], file_list = list.files(path = dir.hvrs["V3_HVR"], pattern = ".tre"))
# V4_paths = file.path(dir.hvrs["V4_HVR"], file_list = list.files(path = dir.hvrs["V4_HVR"], pattern = ".tre"))
# V5_paths = file.path(dir.hvrs["V5_HVR"], file_list = list.files(path = dir.hvrs["V5_HVR"], pattern = ".tre"))
# V6_paths = file.path(dir.hvrs["V6_HVR"], file_list = list.files(path = dir.hvrs["V6_HVR"], pattern = ".tre"))
# V7_paths = file.path(dir.hvrs["V7_HVR"], file_list = list.files(path = dir.hvrs["V7_HVR"], pattern = ".tre"))
# V8_paths = file.path(dir.hvrs["V8_HVR"], file_list = list.files(path = dir.hvrs["V8_HVR"], pattern = ".tre"))
# V9_paths = file.path(dir.hvrs["V9_HVR"], file_list = list.files(path = dir.hvrs["V9_HVR"], pattern = ".tre"))
# 
# # Get guide sequence paths.
# VFull_guides_paths = file.path(dir.hvrs["VFull_guides"], file_list = list.files(path = dir.hvrs["VFull_guides"], pattern = "align.tre"))
# V1_guides_paths = file.path(dir.hvrs["V1_HVR_guides"], file_list = list.files(path = dir.hvrs["V1_HVR_guides"], pattern = "align.tre"))
# V2_guides_paths = file.path(dir.hvrs["V2_HVR_guides"], file_list = list.files(path = dir.hvrs["V2_HVR_guides"], pattern = "align.tre"))
# V3_guides_paths = file.path(dir.hvrs["V3_HVR_guides"], file_list = list.files(path = dir.hvrs["V3_HVR_guides"], pattern = "align.tre"))
# V4_guides_paths = file.path(dir.hvrs["V4_HVR_guides"], file_list = list.files(path = dir.hvrs["V4_HVR_guides"], pattern = "align.tre"))
# V5_guides_paths = file.path(dir.hvrs["V5_HVR_guides"], file_list = list.files(path = dir.hvrs["V5_HVR_guides"], pattern = "align.tre"))
# V6_guides_paths = file.path(dir.hvrs["V6_HVR_guides"], file_list = list.files(path = dir.hvrs["V6_HVR_guides"], pattern = "align.tre"))
# V7_guides_paths = file.path(dir.hvrs["V7_HVR_guides"], file_list = list.files(path = dir.hvrs["V7_HVR_guides"], pattern = "align.tre"))
# V8_guides_paths = file.path(dir.hvrs["V8_HVR_guides"], file_list = list.files(path = dir.hvrs["V8_HVR_guides"], pattern = "align.tre"))
# V9_guides_paths = file.path(dir.hvrs["V9_HVR_guides"], file_list = list.files(path = dir.hvrs["V9_HVR_guides"], pattern = "align.tre"))
# 
# # Remove guide sequences from trees.
# VFull_guides_paths = write_ref_dropped_trees(file_list = VFull_guides_paths, pattern = "^REF.")
# V1_guides_paths = write_ref_dropped_trees(file_list = V1_guides_paths, pattern = "^REF.")
# V2_guides_paths = write_ref_dropped_trees(file_list = V2_guides_paths, pattern = "^REF.")
# V3_guides_paths = write_ref_dropped_trees(file_list = V3_guides_paths, pattern = "^REF.")
# V4_guides_paths = write_ref_dropped_trees(file_list = V4_guides_paths, pattern = "^REF.")
# V5_guides_paths = write_ref_dropped_trees(file_list = V5_guides_paths, pattern = "^REF.")
# V6_guides_paths = write_ref_dropped_trees(file_list = V6_guides_paths, pattern = "^REF.")
# V7_guides_paths = write_ref_dropped_trees(file_list = V7_guides_paths, pattern = "^REF.")
# V8_guides_paths = write_ref_dropped_trees(file_list = V8_guides_paths, pattern = "^REF.")
# V9_guides_paths = write_ref_dropped_trees(file_list = V9_guides_paths, pattern = "^REF.")
# 
# 
# CORES = 10
# vfull_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = VFull_guides_paths, cores = CORES)
# vfull_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = VFull_control_paths, cores = CORES)
# 
# v1_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V1_paths, cores = CORES)
# v2_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V2_paths, cores = CORES)
# v3_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V3_paths, cores = CORES)
# v4_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V4_paths, cores = CORES)
# v5_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V5_paths, cores = CORES)
# v6_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V6_paths, cores = CORES)
# v7_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V7_paths, cores = CORES)
# v8_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V8_paths, cores = CORES)
# v9_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V9_paths, cores = CORES)
# 
# v1_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V1_guides_paths, cores = CORES)
# v2_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V2_guides_paths, cores = CORES)
# v3_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V3_guides_paths, cores = CORES)
# v4_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V4_guides_paths, cores = CORES)
# v5_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V5_guides_paths, cores = CORES)
# v6_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V6_guides_paths, cores = CORES)
# v7_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V7_guides_paths, cores = CORES)
# v8_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V8_guides_paths, cores = CORES)
# v9_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V9_guides_paths, cores = CORES)
# 
# # Make Data Frame Of normalized RF
# hvr_quartet_random_pairs =
#   data.frame("V1" = v1_quartet,
#              "V2" = v2_quartet,
#              "V3" = v3_quartet,
#              "V4" = v4_quartet,
#              "V5" = v5_quartet,
#              "V6" = v6_quartet,
#              "V7" = v7_quartet,
#              "V8" = v8_quartet,
#              "V9" = v9_quartet,
#              "VFull" = vfull_quartet,
#              "V1_guides" = v1_guides_quartet,
#              "V2_guides" = v2_guides_quartet,
#              "V3_guides" = v3_guides_quartet,
#              "V4_guides" = v4_guides_quartet,
#              "V5_guides" = v5_guides_quartet,
#              "V6_guides" = v6_guides_quartet,
#              "V7_guides" = v7_guides_quartet,
#              "V8_guides" = v8_guides_quartet,
#              "V9_guides" = v9_guides_quartet,
#              "VFull_guides" = vfull_guides_quartet)

setwd(dir.out)
#saveRDS(object = hvr_quartet_random_pairs, "hvr_quartet_random_pairs.rds", compress = TRUE)
hvr_quartet_random_pairs = readRDS("hvr_quartet_random_pairs.rds")
dim(hvr_quartet_random_pairs)
#saveRDS(object = hvr_quartet_random_pairs/QUARTET_MAX, "hvr_quartet_normalized_random_pairs.rds", compress = TRUE)
hvr_quartet_norm_random_pairs = readRDS("hvr_quartet_normalized_random_pairs.rds")
dim(hvr_quartet_norm_random_pairs)


#########################
# 1C. Path Distance -- Calculate between 100 random pairings. 
#########################
# hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_path_random_pairings.R' -r hvr_path_random_pairings -p 100 -q sharpton

setwd(dir.out)
hvr_path_random_pairings = readRDS("hvr_path_random_pairings.rds")
head(hvr_path_random_pairings)



#########################
# 1D. Phylogenetic Information Distance (PID) - Random pairings.
#########################

#hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_pid_random_pairs.R' -r pid_random_pairs -p 100 -q sharpton
setwd(dir.out)
hvr_pid_random_pairings = readRDS("hvr_pid_random_pairs.rds")


#########################
# 1E. MSD --  random pairs
#########################

hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_msd_random_pairs.R' -r msd_random_pairs -p 100 -q sharpton

setwd(dir.out)
hvr_msd_random_pairings = readRDS("hvr_msd_random_pairs.rds")


#########################
# 1Fii. NYE Random Pairings
#########################
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_nye_random_pairs.R' -r nye_random_pairs -p 100 -q sharpton
readRDS("hvr_nye_random_pairs.rds")


#########################
# 1G. JRF - Random pairings
#########################
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_jrf_random_pairs.R' -r jrf_random_pairs -p 60 -q biomed
readRDS("hvr_jrf_random_pairs.rds")


#######################################
# 1H. MCD --Random pairings
#######################################
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_mcd_random_pairs.R' -r mcd_random_pairs -p 100 -q sharpton
readRDS("hvr_mcd_random_pairs.rds")

#######################################
# 1J. Determine subcluster distance - Random pairs
#######################################
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_nanostructure_random_pairs.R' -r  nanoclusters_random_paris -p 100 -q sharpton

s2_VFull = calculate_subtree_size_2(trees1 = VFull, cores = 50)
s2_V1 = calculate_subtree_size_2(trees1 = V1, cores = 50)
s2_V2 = calculate_subtree_size_2(trees1 = V2, cores = 50)
s2_V3 = calculate_subtree_size_2(trees1 = V3, cores = 50)
s2_V4 = calculate_subtree_size_2(trees1 = V4, cores = 50)
s2_V5 = calculate_subtree_size_2(trees1 = V5, cores = 50)
s2_V6 = calculate_subtree_size_2(trees1 = V6, cores = 50)
s2_V7 = calculate_subtree_size_2(trees1 = V7, cores = 50)
s2_V8 = calculate_subtree_size_2(trees1 = V8, cores = 50)
s2_V9 = calculate_subtree_size_2(trees1 = V9, cores = 50)

s2_V1_guides = calculate_subtree_size_2(trees1 = V1_guides, cores = 50)
s2_V2_guides = calculate_subtree_size_2(trees1 = V2_guides, cores = 50)
s2_V3_guides = calculate_subtree_size_2(trees1 = V3_guides, cores = 50)
s2_V4_guides = calculate_subtree_size_2(trees1 = V4_guides, cores = 50)
s2_V5_guides = calculate_subtree_size_2(trees1 = V5_guides, cores = 50)
s2_V6_guides = calculate_subtree_size_2(trees1 = V6_guides, cores = 50)
s2_V7_guides = calculate_subtree_size_2(trees1 = V7_guides, cores = 50)
s2_V8_guides = calculate_subtree_size_2(trees1 = V8_guides, cores = 50)
s2_V9_guides = calculate_subtree_size_2(trees1 = V9_guides, cores = 50)

sizes_2 = data.frame("VFull" = s2_VFull,
                     "V1" = s2_V1,
                     "V2" = s2_V2,
                     "V3" = s2_V3,
                     "V4" = s2_V4,
                     "V5" = s2_V5,
                     "V6" = s2_V6,
                     "V7" = s2_V7,
                     "V8" = s2_V8,
                     "V9" = s2_V9,
                     "V1_guides" = s2_V1_guides,
                     "V2_guides" = s2_V2_guides,
                     "V3_guides" = s2_V3_guides,
                     "V4_guides" = s2_V4_guides,
                     "V5_guides" = s2_V5_guides,
                     "V6_guides" = s2_V6_guides,
                     "V7_guides" = s2_V7_guides,
                     "V8_guides" = s2_V8_guides,
                     "V9_guides" = s2_V9_guides)
setwd(dir.out)
#saveRDS(file = "hvr_sizes_nanocluster_2.rds", sizes_2, compress = TRUE)
dplyr::as_tibble(readRDS("hvr_s2_random_pairs.rds"))
dplyr::as_tibble(readRDS("hvr_s3_random_pairs.rds"))
dplyr::as_tibble(readRDS("hvr_s4_random_pairs.rds"))
dplyr::as_tibble(readRDS("hvr_s5_random_pairs.rds"))
dplyr::as_tibble(readRDS("hvr_s6_random_pairs.rds"))
dplyr::as_tibble(readRDS("hvr_s7_random_pairs.rds"))
dplyr::as_tibble(readRDS("hvr_s8_random_pairs.rds"))
dplyr::as_tibble(readRDS("hvr_s9_random_pairs.rds"))
dplyr::as_tibble(readRDS("hvr_s10_random_pairs.rds"))


#######################################
# 1J.Determine MAST
#######################################

hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_hvr_mast.R' -r  mast -p 100 -q sharpton
dplyr::as_tibble(readRDS("hvr_mast_random_pairs.rds"))

#######################################
# Clean up for simulation 1
#######################################
rm(list = ls())

################################################################################
# SIM 2: How does sequence length impact ability to topologically reconstruct trees?
################################################################################

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

# Check that all the tree files are present. 
directories <- c(list.dirs(dir.mix, full.names = TRUE, recursive = FALSE),
                 list.dirs(dir.simulation, full.names = TRUE, recursive = FALSE))
names(directories) = basename(directories)
basename(directories)

################################################################################
# 0: READ IN FILES
################################################################################
N = 100 # The number of trees to look at per simulation
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

################################################################################
# 0B: Check sequence length and tree sizes are what we expect.
################################################################################


if(all(c(
  check_tree_size(tree_list = V1_L50, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V1_L100, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V1_L200, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V1_L300, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V1_L400, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V1_L500, expected_tree_size = TREE_SIZE),
  
  check_tree_size(tree_list = V1_L50_guides, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V1_L100_guides, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V1_L200_guides, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V1_L300_guides, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V1_L400_guides, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V1_L500_guides, expected_tree_size = TREE_SIZE)))){
  print("PASS: All of V1 length simulations are same size")
}else{
  print("FAIL: Not all of V1 length simulations are same size")
}

if(all(c(
  check_tree_size(tree_list = V2_L50, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V2_L100, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V2_L200, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V2_L300, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V2_L400, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V2_L500, expected_tree_size = TREE_SIZE),
  
  check_tree_size(tree_list = V2_L50_guides, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V2_L100_guides, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V2_L200_guides, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V2_L300_guides, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V2_L400_guides, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V2_L500_guides, expected_tree_size = TREE_SIZE)))){
  print("PASS: All of V2 length simulations are same size")
}else{
  print("FAIL: Not all of V2 length simulations are same size")
}

if(all(c(
  check_tree_size(tree_list = V3_L50, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V3_L100, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V3_L200, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V3_L300, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V3_L400, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V3_L500, expected_tree_size = TREE_SIZE),
  
  check_tree_size(tree_list = V3_L50_guides, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V3_L100_guides, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V3_L200_guides, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V3_L300_guides, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V3_L400_guides, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V3_L500_guides, expected_tree_size = TREE_SIZE)))){
  print("PASS: All of V3 length simulations are same size")
}else{
  print("FAIL: Not all of V3 length simulations are same size")
}

if(all(c(
  check_tree_size(tree_list = V4_L50, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V4_L100, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V4_L200, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V4_L300, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V4_L400, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V4_L500, expected_tree_size = TREE_SIZE),
  
  check_tree_size(tree_list = V4_L50_guides, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V4_L100_guides, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V4_L200_guides, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V4_L300_guides, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V4_L400_guides, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V4_L500_guides, expected_tree_size = TREE_SIZE)))){
  print("PASS: All of V4 length simulations are same size")
}else{
  print("FAIL: Not all of V4 length simulations are same size")
}

if(all(c(
  check_tree_size(tree_list = V5_L50, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V5_L100, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V5_L200, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V5_L300, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V5_L400, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V5_L500, expected_tree_size = TREE_SIZE),
  
  check_tree_size(tree_list = V5_L50_guides, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V5_L100_guides, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V5_L200_guides, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V5_L300_guides, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V5_L400_guides, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V5_L500_guides, expected_tree_size = TREE_SIZE)))){
  print("PASS: All of V5 length simulations are same size")
}else{
  print("FAIL: Not all of V5 length simulations are same size")
}


if(all(c(
  check_tree_size(tree_list = V6_L50, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V6_L100, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V6_L200, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V6_L300, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V6_L400, expected_tree_size = TREE_SIZE),
  
  check_tree_size(tree_list = V6_L50_guides, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V6_L100_guides, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V6_L200_guides, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V6_L300_guides, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V6_L400_guides, expected_tree_size = TREE_SIZE)))){
  print("PASS: All of V6 length simulations are same size")
}else{
  print("FAIL: Not all of V6  length simulations are same size")
}

if(all(c(
  check_tree_size(tree_list = V7_L50, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V7_L100, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V7_L200, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V7_L300, expected_tree_size = TREE_SIZE),
  
  check_tree_size(tree_list = V7_L50_guides, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V7_L100_guides, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V7_L200_guides, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V7_L300_guides, expected_tree_size = TREE_SIZE)))){
  print("PASS: All of V7 length simulations are same size")
}else{
  print("FAIL: Not all of V7 length simulations are same size")
}

if(all(c(
  check_tree_size(tree_list = V8_L50, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V8_L100, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V8_L200, expected_tree_size = TREE_SIZE),
  
  check_tree_size(tree_list = V8_L50_guides, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V8_L100_guides, expected_tree_size = TREE_SIZE),
  check_tree_size(tree_list = V8_L200_guides, expected_tree_size = TREE_SIZE)))){
  print("PASS: All of V8 length simulations are same size")
}else{
  print("FAIL: Not all of V8 length simulations are same size")
}


################################################################################
# 1: Determine how guides x length impact topological reconstruction -- RF
################################################################################
CORES = 20
VFull_control = calculate_rf_distance_random_pairs(VFull, VFull_control, cores = CORES)
VFull_guides = calculate_rf_distance_random_pairs(VFull, VFull_guides, cores = CORES)
V1_L50_rf = calculate_rf_distance_random_pairs(VFull, V1_L50, cores = CORES)
V1_L100_rf = calculate_rf_distance_random_pairs(VFull, V1_L100, cores = CORES)
V1_L200_rf = calculate_rf_distance_random_pairs(VFull, V1_L200, cores = CORES)
V1_L300_rf = calculate_rf_distance_random_pairs(VFull, V1_L300, cores = CORES)
V1_L400_rf = calculate_rf_distance_random_pairs(VFull, V1_L400, cores = CORES)
V1_L500_rf = calculate_rf_distance_random_pairs(VFull, V1_L500, cores = CORES)

V1_L50_guides_rf = calculate_rf_distance_random_pairs(VFull, V1_L50_guides, cores = CORES)
V1_L100_guides_rf = calculate_rf_distance_random_pairs(VFull, V1_L100_guides, cores = CORES)
V1_L200_guides_rf = calculate_rf_distance_random_pairs(VFull, V1_L200_guides, cores = CORES)
V1_L300_guides_rf = calculate_rf_distance_random_pairs(VFull, V1_L300_guides, cores = CORES)
V1_L400_guides_rf = calculate_rf_distance_random_pairs(VFull, V1_L400_guides, cores = CORES)
V1_L500_guides_rf = calculate_rf_distance_random_pairs(VFull, V1_L500_guides, cores = CORES)

V2_L50_rf = calculate_rf_distance_random_pairs(VFull, V2_L50, cores = CORES)
V2_L100_rf = calculate_rf_distance_random_pairs(VFull, V2_L100, cores = CORES)
V2_L200_rf = calculate_rf_distance_random_pairs(VFull, V2_L200, cores = CORES)
V2_L300_rf = calculate_rf_distance_random_pairs(VFull, V2_L300, cores = CORES)
V2_L400_rf = calculate_rf_distance_random_pairs(VFull, V2_L400, cores = CORES)
V2_L500_rf = calculate_rf_distance_random_pairs(VFull, V2_L500, cores = CORES)

V2_L50_guides_rf = calculate_rf_distance_random_pairs(VFull, V2_L50_guides, cores = CORES)
V2_L100_guides_rf = calculate_rf_distance_random_pairs(VFull, V2_L100_guides, cores = CORES)
V2_L200_guides_rf = calculate_rf_distance_random_pairs(VFull, V2_L200_guides, cores = CORES)
V2_L300_guides_rf = calculate_rf_distance_random_pairs(VFull, V2_L300_guides, cores = CORES)
V2_L400_guides_rf = calculate_rf_distance_random_pairs(VFull, V2_L400_guides, cores = CORES)
V2_L500_guides_rf = calculate_rf_distance_random_pairs(VFull, V2_L500_guides, cores = CORES)

V3_L50_rf = calculate_rf_distance_random_pairs(VFull, V3_L50, cores = CORES)
V3_L100_rf = calculate_rf_distance_random_pairs(VFull, V3_L100, cores = CORES)
V3_L200_rf = calculate_rf_distance_random_pairs(VFull, V3_L200, cores = CORES)
V3_L300_rf = calculate_rf_distance_random_pairs(VFull, V3_L300, cores = CORES)
V3_L400_rf = calculate_rf_distance_random_pairs(VFull, V3_L400, cores = CORES)
V3_L500_rf = calculate_rf_distance_random_pairs(VFull, V3_L500, cores = CORES)

V3_L50_guides_rf = calculate_rf_distance_random_pairs(VFull, V3_L50_guides, cores = CORES)
V3_L100_guides_rf = calculate_rf_distance_random_pairs(VFull, V3_L100_guides, cores = CORES)
V3_L200_guides_rf = calculate_rf_distance_random_pairs(VFull, V3_L200_guides, cores = CORES)
V3_L300_guides_rf = calculate_rf_distance_random_pairs(VFull, V3_L300_guides, cores = CORES)
V3_L400_guides_rf = calculate_rf_distance_random_pairs(VFull, V3_L400_guides, cores = CORES)
V3_L500_guides_rf = calculate_rf_distance_random_pairs(VFull, V3_L500_guides, cores = CORES)

V4_L50_rf = calculate_rf_distance_random_pairs(VFull, V4_L50, cores = CORES)
V4_L100_rf = calculate_rf_distance_random_pairs(VFull, V4_L100, cores = CORES)
V4_L200_rf = calculate_rf_distance_random_pairs(VFull, V4_L200, cores = CORES)
V4_L300_rf = calculate_rf_distance_random_pairs(VFull, V4_L300, cores = CORES)
V4_L400_rf = calculate_rf_distance_random_pairs(VFull, V4_L400, cores = CORES)
V4_L500_rf = calculate_rf_distance_random_pairs(VFull, V4_L500, cores = CORES)

V4_L50_guides_rf = calculate_rf_distance_random_pairs(VFull, V4_L50_guides, cores = CORES)
V4_L100_guides_rf = calculate_rf_distance_random_pairs(VFull, V4_L100_guides, cores = CORES)
V4_L200_guides_rf = calculate_rf_distance_random_pairs(VFull, V4_L200_guides, cores = CORES)
V4_L300_guides_rf = calculate_rf_distance_random_pairs(VFull, V4_L300_guides, cores = CORES)
V4_L400_guides_rf = calculate_rf_distance_random_pairs(VFull, V4_L400_guides, cores = CORES)
V4_L500_guides_rf = calculate_rf_distance_random_pairs(VFull, V4_L500_guides, cores = CORES)

V5_L50_rf = calculate_rf_distance_random_pairs(VFull, V5_L50, cores = CORES)
V5_L100_rf = calculate_rf_distance_random_pairs(VFull, V5_L100, cores = CORES)
V5_L200_rf = calculate_rf_distance_random_pairs(VFull, V5_L200, cores = CORES)
V5_L300_rf = calculate_rf_distance_random_pairs(VFull, V5_L300, cores = CORES)
V5_L400_rf = calculate_rf_distance_random_pairs(VFull, V5_L400, cores = CORES)
V5_L500_rf = calculate_rf_distance_random_pairs(VFull, V5_L500, cores = CORES)

V5_L50_guides_rf = calculate_rf_distance_random_pairs(VFull, V5_L50_guides, cores = CORES)
V5_L100_guides_rf = calculate_rf_distance_random_pairs(VFull, V5_L100_guides, cores = CORES)
V5_L200_guides_rf = calculate_rf_distance_random_pairs(VFull, V5_L200_guides, cores = CORES)
V5_L300_guides_rf = calculate_rf_distance_random_pairs(VFull, V5_L300_guides, cores = CORES)
V5_L400_guides_rf = calculate_rf_distance_random_pairs(VFull, V5_L400_guides, cores = CORES)
V5_L500_guides_rf = calculate_rf_distance_random_pairs(VFull, V5_L500_guides, cores = CORES)

V6_L50_rf = calculate_rf_distance_random_pairs(VFull, V6_L50, cores = CORES)
V6_L100_rf = calculate_rf_distance_random_pairs(VFull, V6_L100, cores = CORES)
V6_L200_rf = calculate_rf_distance_random_pairs(VFull, V6_L200, cores = CORES)
V6_L300_rf = calculate_rf_distance_random_pairs(VFull, V6_L300, cores = CORES)
V6_L400_rf = calculate_rf_distance_random_pairs(VFull, V6_L400, cores = CORES)

V6_L50_guides_rf = calculate_rf_distance_random_pairs(VFull, V6_L50_guides, cores = CORES)
V6_L100_guides_rf = calculate_rf_distance_random_pairs(VFull, V6_L100_guides, cores = CORES)
V6_L200_guides_rf = calculate_rf_distance_random_pairs(VFull, V6_L200_guides, cores = CORES)
V6_L300_guides_rf = calculate_rf_distance_random_pairs(VFull, V6_L300_guides, cores = CORES)
V6_L400_guides_rf = calculate_rf_distance_random_pairs(VFull, V6_L400_guides, cores = CORES)

V7_L50_rf = calculate_rf_distance_random_pairs(VFull, V7_L50, cores = CORES)
V7_L100_rf = calculate_rf_distance_random_pairs(VFull, V7_L100, cores = CORES)
V7_L200_rf = calculate_rf_distance_random_pairs(VFull, V7_L200, cores = CORES)
V7_L300_rf = calculate_rf_distance_random_pairs(VFull, V7_L300, cores = CORES)

V7_L50_guides_rf = calculate_rf_distance_random_pairs(VFull, V7_L50_guides, cores = CORES)
V7_L100_guides_rf = calculate_rf_distance_random_pairs(VFull, V7_L100_guides, cores = CORES)
V7_L200_guides_rf = calculate_rf_distance_random_pairs(VFull, V7_L200_guides, cores = CORES)
V7_L300_guides_rf = calculate_rf_distance_random_pairs(VFull, V7_L300_guides, cores = CORES)

V8_L50_rf = calculate_rf_distance_random_pairs(VFull, V8_L50, cores = CORES)
V8_L100_rf = calculate_rf_distance_random_pairs(VFull, V8_L100, cores = CORES)
V8_L200_rf = calculate_rf_distance_random_pairs(VFull, V8_L200, cores = CORES)

V8_L50_guides_rf = calculate_rf_distance_random_pairs(VFull, V8_L50_guides, cores = CORES)
V8_L100_guides_rf = calculate_rf_distance_random_pairs(VFull, V8_L100_guides, cores = CORES)
V8_L200_guides_rf = calculate_rf_distance_random_pairs(VFull, V8_L200_guides, cores = CORES)

bp_rf = 
  data.frame("VFull_control" = VFull_control,
             "VFull_guides" = VFull_guides,
             "V1_L50_rf" = V1_L50_rf,
             "V1_L100_rf" = V1_L100_rf,
             "V1_L200_rf" = V1_L200_rf,
             "V1_L300_rf" = V1_L300_rf,
             "V1_L400_rf" = V1_L400_rf,
             "V1_L500_rf" = V1_L500_rf,
             "V1_L50_guides_rf" = V1_L50_guides_rf,
             "V1_L100_guides_rf" = V1_L100_guides_rf,
             "V1_L200_guides_rf" = V1_L200_guides_rf,
             "V1_L300_guides_rf" = V1_L300_guides_rf,
             "V1_L400_guides_rf" = V1_L400_guides_rf,
             "V1_L500_guides_rf" = V1_L500_guides_rf,
             "V2_L50_rf" = V2_L50_rf,
             "V2_L100_rf" = V2_L100_rf,
             "V2_L200_rf" = V2_L200_rf,
             "V2_L300_rf" = V2_L300_rf,
             "V2_L400_rf" = V2_L400_rf,
             "V2_L500_rf" = V2_L500_rf,
             "V2_L50_guides_rf" = V2_L50_guides_rf,
             "V2_L100_guides_rf" = V2_L100_guides_rf,
             "V2_L200_guides_rf" = V2_L200_guides_rf,
             "V2_L300_guides_rf" = V2_L300_guides_rf,
             "V2_L400_guides_rf" = V2_L400_guides_rf,
             "V2_L500_guides_rf" = V2_L500_guides_rf,
             "V3_L50_rf" = V3_L50_rf,
             "V3_L100_rf" = V3_L100_rf,
             "V3_L200_rf" = V3_L200_rf,
             "V3_L300_rf" = V3_L300_rf,
             "V3_L400_rf" = V3_L400_rf,
             "V3_L500_rf" = V3_L500_rf,
             "V3_L50_guides_rf" = V3_L50_guides_rf,
             "V3_L100_guides_rf" = V3_L100_guides_rf,
             "V3_L200_guides_rf" = V3_L200_guides_rf,
             "V3_L300_guides_rf" = V3_L300_guides_rf,
             "V3_L400_guides_rf" = V3_L400_guides_rf,
             "V3_L500_guides_rf" = V3_L500_guides_rf,
             "V4_L50_rf" = V4_L50_rf,
             "V4_L100_rf" = V4_L100_rf,
             "V4_L200_rf" = V4_L200_rf,
             "V4_L300_rf" = V4_L300_rf,
             "V4_L400_rf" = V4_L400_rf,
             "V4_L500_rf" = V4_L500_rf,
             "V4_L50_guides_rf" = V4_L50_guides_rf,
             "V4_L100_guides_rf" = V4_L100_guides_rf,
             "V4_L200_guides_rf" = V4_L200_guides_rf,
             "V4_L300_guides_rf" = V4_L300_guides_rf,
             "V4_L400_guides_rf" = V4_L400_guides_rf,
             "V4_L500_guides_rf" = V4_L500_guides_rf,
             "V5_L50_rf" = V5_L50_rf,
             "V5_L100_rf" = V5_L100_rf,
             "V5_L200_rf" = V5_L200_rf,
             "V5_L300_rf" = V5_L300_rf,
             "V5_L400_rf" = V5_L400_rf,
             "V5_L500_rf" = V5_L500_rf,
             "V5_L50_guides_rf" = V5_L50_guides_rf,
             "V5_L100_guides_rf" = V5_L100_guides_rf,
             "V5_L200_guides_rf" = V5_L200_guides_rf,
             "V5_L300_guides_rf" = V5_L300_guides_rf,
             "V5_L400_guides_rf" = V5_L400_guides_rf,
             "V5_L500_guides_rf" = V5_L500_guides_rf,
             "V6_L50_rf" = V6_L50_rf,
             "V6_L100_rf" = V6_L100_rf,
             "V6_L200_rf" = V6_L200_rf,
             "V6_L300_rf" = V6_L300_rf,
             "V6_L400_rf" = V6_L400_rf,
             "V6_L50_guides_rf" = V6_L50_guides_rf,
             "V6_L100_guides_rf" = V6_L100_guides_rf,
             "V6_L200_guides_rf" = V6_L200_guides_rf,
             "V6_L300_guides_rf" = V6_L300_guides_rf,
             "V6_L400_guides_rf" = V6_L400_guides_rf,
             "V7_L50_rf" = V7_L50_rf,
             "V7_L100_rf" = V7_L100_rf,
             "V7_L200_rf" = V7_L200_rf,
             "V7_L300_rf" = V7_L300_rf,
             "V7_L50_guides_rf" = V7_L50_guides_rf,
             "V7_L100_guides_rf" = V7_L100_guides_rf,
             "V7_L200_guides_rf" = V7_L200_guides_rf,
             "V7_L300_guides_rf" = V7_L300_guides_rf,
             "V8_L50_rf" = V8_L50_rf,
             "V8_L100_rf" = V8_L100_rf,
             "V8_L200_rf" = V8_L200_rf,
             "V8_L50_guides_rf" = V8_L50_guides_rf,
             "V8_L100_guides_rf" = V8_L100_guides_rf,
             "V8_L200_guides_rf" = V8_L200_guides_rf)
setwd(dir.sim.out)
#saveRDS(file = "bp_rf_random_pairs.rds", object = bp_rf)
#bp_rf = readRDS("bp_rf_random_pairs.rds")
#bp_rf
dplyr::as_tibble(readRDS("bp_rf_random_pairs.rds"))

################################################################################
# 1B: Determine how guides x length impact topological reconstruction -- JRF
################################################################################

hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_bp_jrf_random_pairs.R' -r bp_jrf_random_pairs -p 100 -q sharpton
dplyr::as_tibble(readRDS("bp_jrf_random_pairs.rds"))

################################################################################
# 1C: Determine how guides x length impact topological reconstruction -- PID
################################################################################
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_bp_pid_random_pairs.R' -r bp_pid_random_pairs -p 100 -q sharpton
dplyr::as_tibble(readRDS("bp_pid_random_pairs.rds"))

################################################################################
# 1D: Determine how guides x length impact topological reconstruction -- MSD
################################################################################

hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_bp_msd_random_pairs.R' -r bp_msd_random_pairs -p 30 -q biomed
dplyr::as_tibble(readRDS("bp_msd_random_pairs.rds"))

################################################################################
# 1E: Determine how guides x length impact topological reconstruction -- QD
################################################################################

VFull_paths = file.path(dir.hvrs["VFull"], file_list = list.files(path = dir.hvrs["VFull"], pattern = ".tre"))[1:N]
VFull_guides_paths = file.path(dir.hvrs["VFull_guides"], file_list = list.files(path = dir.hvrs["VFull_guides"], pattern = "align.tre"), drop = TRUE)[1:N]
VFull_control_paths = file.path(dir.hvrs["VFull_control"], file_list = list.files(path = dir.hvrs["VFull_control"], pattern = "align.tre"))[1:N]

V1_L50_paths = file.path(dir.hvrs["V1_L50"], file_list = list.files(path = dir.hvrs["V1_L50"], pattern = ".tre"))[1:N]
V2_L50_paths = file.path(dir.hvrs["V2_L50"], file_list = list.files(path = dir.hvrs["V2_L50"], pattern = ".tre"))[1:N]
V3_L50_paths = file.path(dir.hvrs["V3_L50"], file_list = list.files(path = dir.hvrs["V3_L50"], pattern = ".tre"))[1:N]
V4_L50_paths = file.path(dir.hvrs["V4_L50"], file_list = list.files(path = dir.hvrs["V4_L50"], pattern = ".tre"))[1:N]
V5_L50_paths = file.path(dir.hvrs["V5_L50"], file_list = list.files(path = dir.hvrs["V5_L50"], pattern = ".tre"))[1:N]
V6_L50_paths = file.path(dir.hvrs["V6_L50"], file_list = list.files(path = dir.hvrs["V6_L50"], pattern = ".tre"))[1:N]
V7_L50_paths = file.path(dir.hvrs["V7_L50"], file_list = list.files(path = dir.hvrs["V7_L50"], pattern = ".tre"))[1:N]
V8_L50_paths = file.path(dir.hvrs["V8_L50"], file_list = list.files(path = dir.hvrs["V8_L50"], pattern = ".tre"))[1:N]

V1_L100_paths = file.path(dir.hvrs["V1_L100"], file_list = list.files(path = dir.hvrs["V1_L100"], pattern = ".tre"))[1:N]
V2_L100_paths = file.path(dir.hvrs["V2_L100"], file_list = list.files(path = dir.hvrs["V2_L100"], pattern = ".tre"))[1:N]
V3_L100_paths = file.path(dir.hvrs["V3_L100"], file_list = list.files(path = dir.hvrs["V3_L100"], pattern = ".tre"))[1:N]
V4_L100_paths = file.path(dir.hvrs["V4_L100"], file_list = list.files(path = dir.hvrs["V4_L100"], pattern = ".tre"))[1:N]
V5_L100_paths = file.path(dir.hvrs["V5_L100"], file_list = list.files(path = dir.hvrs["V5_L100"], pattern = ".tre"))[1:N]
V6_L100_paths = file.path(dir.hvrs["V6_L100"], file_list = list.files(path = dir.hvrs["V6_L100"], pattern = ".tre"))[1:N]
V7_L100_paths = file.path(dir.hvrs["V7_L100"], file_list = list.files(path = dir.hvrs["V7_L100"], pattern = ".tre"))[1:N]
V8_L100_paths = file.path(dir.hvrs["V8_L100"], file_list = list.files(path = dir.hvrs["V8_L100"], pattern = ".tre"))[1:N]

V1_L200_paths = file.path(dir.hvrs["V1_L200"], file_list = list.files(path = dir.hvrs["V1_L200"], pattern = ".tre"))[1:N]
V2_L200_paths = file.path(dir.hvrs["V2_L200"], file_list = list.files(path = dir.hvrs["V2_L200"], pattern = ".tre"))[1:N]
V3_L200_paths = file.path(dir.hvrs["V3_L200"], file_list = list.files(path = dir.hvrs["V3_L200"], pattern = ".tre"))[1:N]
V4_L200_paths = file.path(dir.hvrs["V4_L200"], file_list = list.files(path = dir.hvrs["V4_L200"], pattern = ".tre"))[1:N]
V5_L200_paths = file.path(dir.hvrs["V5_L200"], file_list = list.files(path = dir.hvrs["V5_L200"], pattern = ".tre"))[1:N]
V6_L200_paths = file.path(dir.hvrs["V6_L200"], file_list = list.files(path = dir.hvrs["V6_L200"], pattern = ".tre"))[1:N]
V7_L200_paths = file.path(dir.hvrs["V7_L200"], file_list = list.files(path = dir.hvrs["V7_L200"], pattern = ".tre"))[1:N]
V8_L200_paths = file.path(dir.hvrs["V8_L200"], file_list = list.files(path = dir.hvrs["V8_L200"], pattern = ".tre"))[1:N]

V1_L300_paths = file.path(dir.hvrs["V1_L100"], file_list = list.files(path = dir.hvrs["V1_L100"], pattern = ".tre"))[1:N]
V2_L300_paths = file.path(dir.hvrs["V2_L300"], file_list = list.files(path = dir.hvrs["V2_L300"], pattern = ".tre"))[1:N]
V3_L300_paths = file.path(dir.hvrs["V3_L300"], file_list = list.files(path = dir.hvrs["V3_L300"], pattern = ".tre"))[1:N]
V4_L300_paths = file.path(dir.hvrs["V4_L300"], file_list = list.files(path = dir.hvrs["V4_L300"], pattern = ".tre"))[1:N]
V5_L300_paths = file.path(dir.hvrs["V5_L300"], file_list = list.files(path = dir.hvrs["V5_L300"], pattern = ".tre"))[1:N]
V6_L300_paths = file.path(dir.hvrs["V6_L300"], file_list = list.files(path = dir.hvrs["V6_L300"], pattern = ".tre"))[1:N]
V7_L300_paths = file.path(dir.hvrs["V7_L300"], file_list = list.files(path = dir.hvrs["V7_L300"], pattern = ".tre"))[1:N]

V1_L400_paths = file.path(dir.hvrs["V1_L100"], file_list = list.files(path = dir.hvrs["V1_L100"], pattern = ".tre"))[1:N]
V2_L400_paths = file.path(dir.hvrs["V2_L400"], file_list = list.files(path = dir.hvrs["V2_L400"], pattern = ".tre"))[1:N]
V3_L400_paths = file.path(dir.hvrs["V3_L400"], file_list = list.files(path = dir.hvrs["V3_L400"], pattern = ".tre"))[1:N]
V4_L400_paths = file.path(dir.hvrs["V4_L400"], file_list = list.files(path = dir.hvrs["V4_L400"], pattern = ".tre"))[1:N]
V5_L400_paths = file.path(dir.hvrs["V5_L400"], file_list = list.files(path = dir.hvrs["V5_L400"], pattern = ".tre"))[1:N]
V6_L400_paths = file.path(dir.hvrs["V6_L400"], file_list = list.files(path = dir.hvrs["V6_L400"], pattern = ".tre"))[1:N]

V1_L500_paths = file.path(dir.hvrs["V1_L100"], file_list = list.files(path = dir.hvrs["V1_L100"], pattern = ".tre"))[1:N]
V2_L500_paths = file.path(dir.hvrs["V2_L500"], file_list = list.files(path = dir.hvrs["V2_L500"], pattern = ".tre"))[1:N]
V3_L500_paths = file.path(dir.hvrs["V3_L500"], file_list = list.files(path = dir.hvrs["V3_L500"], pattern = ".tre"))[1:N]
V4_L500_paths = file.path(dir.hvrs["V4_L500"], file_list = list.files(path = dir.hvrs["V4_L500"], pattern = ".tre"))[1:N]
V5_L500_paths = file.path(dir.hvrs["V5_L500"], file_list = list.files(path = dir.hvrs["V5_L500"], pattern = ".tre"))[1:N]

# # Get guide sequence paths.
VFull_guides_paths = file.path(dir.hvrs["VFull_guides"], file_list = list.files(path = dir.hvrs["VFull_guides"], pattern = "align.tre"))[1:N]
V1_L50_guides_paths = file.path(dir.hvrs["V1_L50_guides"], file_list = list.files(path = dir.hvrs["V1_L50_guides"], pattern = "align.tre"))[1:N]
V2_L50_guides_paths = file.path(dir.hvrs["V2_L50_guides"], file_list = list.files(path = dir.hvrs["V2_L50_guides"], pattern = "align.tre"))[1:N]
V3_L50_guides_paths = file.path(dir.hvrs["V3_L50_guides"], file_list = list.files(path = dir.hvrs["V3_L50_guides"], pattern = "align.tre"))[1:N]
V4_L50_guides_paths = file.path(dir.hvrs["V4_L50_guides"], file_list = list.files(path = dir.hvrs["V4_L50_guides"], pattern = "align.tre"))[1:N]
V5_L50_guides_paths = file.path(dir.hvrs["V5_L50_guides"], file_list = list.files(path = dir.hvrs["V5_L50_guides"], pattern = "align.tre"))[1:N]
V6_L50_guides_paths = file.path(dir.hvrs["V6_L50_guides"], file_list = list.files(path = dir.hvrs["V6_L50_guides"], pattern = "align.tre"))[1:N]
V7_L50_guides_paths = file.path(dir.hvrs["V7_L50_guides"], file_list = list.files(path = dir.hvrs["V7_L50_guides"], pattern = "align.tre"))[1:N]
V8_L50_guides_paths = file.path(dir.hvrs["V8_L50_guides"], file_list = list.files(path = dir.hvrs["V8_L50_guides"], pattern = "align.tre"))[1:N]

V1_L100_guides_paths = file.path(dir.hvrs["V1_L100_guides"], file_list = list.files(path = dir.hvrs["V1_L100_guides"], pattern = "align.tre"))[1:N]
V2_L100_guides_paths = file.path(dir.hvrs["V2_L100_guides"], file_list = list.files(path = dir.hvrs["V2_L100_guides"], pattern = "align.tre"))[1:N]
V3_L100_guides_paths = file.path(dir.hvrs["V3_L100_guides"], file_list = list.files(path = dir.hvrs["V3_L100_guides"], pattern = "align.tre"))[1:N]
V4_L100_guides_paths = file.path(dir.hvrs["V4_L100_guides"], file_list = list.files(path = dir.hvrs["V4_L100_guides"], pattern = "align.tre"))[1:N]
V5_L100_guides_paths = file.path(dir.hvrs["V5_L100_guides"], file_list = list.files(path = dir.hvrs["V5_L100_guides"], pattern = "align.tre"))[1:N]
V6_L100_guides_paths = file.path(dir.hvrs["V6_L100_guides"], file_list = list.files(path = dir.hvrs["V6_L100_guides"], pattern = "align.tre"))[1:N]
V7_L100_guides_paths = file.path(dir.hvrs["V7_L100_guides"], file_list = list.files(path = dir.hvrs["V7_L100_guides"], pattern = "align.tre"))[1:N]
V8_L100_guides_paths = file.path(dir.hvrs["V8_L100_guides"], file_list = list.files(path = dir.hvrs["V8_L100_guides"], pattern = "align.tre"))[1:N]

V1_L200_guides_paths = file.path(dir.hvrs["V1_L200_guides"], file_list = list.files(path = dir.hvrs["V1_L200_guides"], pattern = "align.tre"))[1:N]
V2_L200_guides_paths = file.path(dir.hvrs["V2_L200_guides"], file_list = list.files(path = dir.hvrs["V2_L200_guides"], pattern = "align.tre"))[1:N]
V3_L200_guides_paths = file.path(dir.hvrs["V3_L200_guides"], file_list = list.files(path = dir.hvrs["V3_L200_guides"], pattern = "align.tre"))[1:N]
V4_L200_guides_paths = file.path(dir.hvrs["V4_L200_guides"], file_list = list.files(path = dir.hvrs["V4_L200_guides"], pattern = "align.tre"))[1:N]
V5_L200_guides_paths = file.path(dir.hvrs["V5_L200_guides"], file_list = list.files(path = dir.hvrs["V5_L200_guides"], pattern = "align.tre"))[1:N]
V6_L200_guides_paths = file.path(dir.hvrs["V6_L200_guides"], file_list = list.files(path = dir.hvrs["V6_L200_guides"], pattern = "align.tre"))[1:N]
V7_L200_guides_paths = file.path(dir.hvrs["V7_L200_guides"], file_list = list.files(path = dir.hvrs["V7_L200_guides"], pattern = "align.tre"))[1:N]
V8_L200_guides_paths = file.path(dir.hvrs["V8_L200_guides"], file_list = list.files(path = dir.hvrs["V8_L200_guides"], pattern = "align.tre"))[1:N]

V1_L300_guides_paths = file.path(dir.hvrs["V1_L300_guides"], file_list = list.files(path = dir.hvrs["V1_L300_guides"], pattern = "align.tre"))[1:N]
V2_L300_guides_paths = file.path(dir.hvrs["V2_L300_guides"], file_list = list.files(path = dir.hvrs["V2_L300_guides"], pattern = "align.tre"))[1:N]
V3_L300_guides_paths = file.path(dir.hvrs["V3_L300_guides"], file_list = list.files(path = dir.hvrs["V3_L300_guides"], pattern = "align.tre"))[1:N]
V4_L300_guides_paths = file.path(dir.hvrs["V4_L300_guides"], file_list = list.files(path = dir.hvrs["V4_L300_guides"], pattern = "align.tre"))[1:N]
V5_L300_guides_paths = file.path(dir.hvrs["V5_L300_guides"], file_list = list.files(path = dir.hvrs["V5_L300_guides"], pattern = "align.tre"))[1:N]
V6_L300_guides_paths = file.path(dir.hvrs["V6_L300_guides"], file_list = list.files(path = dir.hvrs["V6_L300_guides"], pattern = "align.tre"))[1:N]
V7_L300_guides_paths = file.path(dir.hvrs["V7_L300_guides"], file_list = list.files(path = dir.hvrs["V7_L300_guides"], pattern = "align.tre"))[1:N]

V1_L400_guides_paths = file.path(dir.hvrs["V1_L400_guides"], file_list = list.files(path = dir.hvrs["V1_L400_guides"], pattern = "align.tre"))[1:N]
V2_L400_guides_paths = file.path(dir.hvrs["V2_L400_guides"], file_list = list.files(path = dir.hvrs["V2_L400_guides"], pattern = "align.tre"))[1:N]
V3_L400_guides_paths = file.path(dir.hvrs["V3_L400_guides"], file_list = list.files(path = dir.hvrs["V3_L400_guides"], pattern = "align.tre"))[1:N]
V4_L400_guides_paths = file.path(dir.hvrs["V4_L400_guides"], file_list = list.files(path = dir.hvrs["V4_L400_guides"], pattern = "align.tre"))[1:N]
V5_L400_guides_paths = file.path(dir.hvrs["V5_L400_guides"], file_list = list.files(path = dir.hvrs["V5_L400_guides"], pattern = "align.tre"))[1:N]
V6_L400_guides_paths = file.path(dir.hvrs["V6_L400_guides"], file_list = list.files(path = dir.hvrs["V6_L400_guides"], pattern = "align.tre"))[1:N]

V1_L500_guides_paths = file.path(dir.hvrs["V1_L500_guides"], file_list = list.files(path = dir.hvrs["V1_L500_guides"], pattern = "align.tre"))[1:N]
V2_L500_guides_paths = file.path(dir.hvrs["V2_L500_guides"], file_list = list.files(path = dir.hvrs["V2_L500_guides"], pattern = "align.tre"))[1:N]
V3_L500_guides_paths = file.path(dir.hvrs["V3_L500_guides"], file_list = list.files(path = dir.hvrs["V3_L500_guides"], pattern = "align.tre"))[1:N]
V4_L500_guides_paths = file.path(dir.hvrs["V4_L500_guides"], file_list = list.files(path = dir.hvrs["V4_L500_guides"], pattern = "align.tre"))[1:N]
V5_L500_guides_paths = file.path(dir.hvrs["V5_L500_guides"], file_list = list.files(path = dir.hvrs["V5_L500_guides"], pattern = "align.tre"))[1:N]

# The reference was already written, so just read it in
VFull_guides_paths = file.path(dir.hvrs["VFull_guides"], file_list = list.files(path = dir.hvrs["VFull_guides"], pattern = "shuffle._refdropped.tre"))[1:N]
VFull_guides_paths

# Make reference dropped trees
V1_L50_guides_paths = write_ref_dropped_trees(file_list = V1_L50_guides_paths, pattern = "^REF.")
V2_L50_guides_paths = write_ref_dropped_trees(file_list = V2_L50_guides_paths, pattern = "^REF.")
V3_L50_guides_paths = write_ref_dropped_trees(file_list = V3_L50_guides_paths, pattern = "^REF.")
V4_L50_guides_paths = write_ref_dropped_trees(file_list = V4_L50_guides_paths, pattern = "^REF.")
V5_L50_guides_paths = write_ref_dropped_trees(file_list = V5_L50_guides_paths, pattern = "^REF.")
V6_L50_guides_paths = write_ref_dropped_trees(file_list = V6_L50_guides_paths, pattern = "^REF.")
V7_L50_guides_paths = write_ref_dropped_trees(file_list = V7_L50_guides_paths, pattern = "^REF.")
V8_L50_guides_paths = write_ref_dropped_trees(file_list = V8_L50_guides_paths, pattern = "^REF.")


V1_L100_guides_paths = write_ref_dropped_trees(file_list = V1_L100_guides_paths, pattern = "^REF.")
V2_L100_guides_paths = write_ref_dropped_trees(file_list = V2_L100_guides_paths, pattern = "^REF.")
V3_L100_guides_paths = write_ref_dropped_trees(file_list = V3_L100_guides_paths, pattern = "^REF.")
V4_L100_guides_paths = write_ref_dropped_trees(file_list = V4_L100_guides_paths, pattern = "^REF.")
V5_L100_guides_paths = write_ref_dropped_trees(file_list = V5_L100_guides_paths, pattern = "^REF.")
V6_L100_guides_paths = write_ref_dropped_trees(file_list = V6_L100_guides_paths, pattern = "^REF.")
V7_L100_guides_paths = write_ref_dropped_trees(file_list = V7_L100_guides_paths, pattern = "^REF.")
V8_L100_guides_paths = write_ref_dropped_trees(file_list = V8_L100_guides_paths, pattern = "^REF.")

V1_L200_guides_paths = write_ref_dropped_trees(file_list = V1_L200_guides_paths, pattern = "^REF.")
V2_L200_guides_paths = write_ref_dropped_trees(file_list = V2_L200_guides_paths, pattern = "^REF.")
V3_L200_guides_paths = write_ref_dropped_trees(file_list = V3_L200_guides_paths, pattern = "^REF.")
V4_L200_guides_paths = write_ref_dropped_trees(file_list = V4_L200_guides_paths, pattern = "^REF.")
V5_L200_guides_paths = write_ref_dropped_trees(file_list = V5_L200_guides_paths, pattern = "^REF.")
V6_L200_guides_paths = write_ref_dropped_trees(file_list = V6_L200_guides_paths, pattern = "^REF.")
V7_L200_guides_paths = write_ref_dropped_trees(file_list = V7_L200_guides_paths, pattern = "^REF.")
V8_L200_guides_paths = write_ref_dropped_trees(file_list = V8_L200_guides_paths, pattern = "^REF.")

V1_L300_guides_paths = write_ref_dropped_trees(file_list = V1_L300_guides_paths, pattern = "^REF.")
V2_L300_guides_paths = write_ref_dropped_trees(file_list = V2_L300_guides_paths, pattern = "^REF.")
V3_L300_guides_paths = write_ref_dropped_trees(file_list = V3_L300_guides_paths, pattern = "^REF.")
V4_L300_guides_paths = write_ref_dropped_trees(file_list = V4_L300_guides_paths, pattern = "^REF.")
V5_L300_guides_paths = write_ref_dropped_trees(file_list = V5_L300_guides_paths, pattern = "^REF.")
V6_L300_guides_paths = write_ref_dropped_trees(file_list = V6_L300_guides_paths, pattern = "^REF.")
V7_L300_guides_paths = write_ref_dropped_trees(file_list = V7_L300_guides_paths, pattern = "^REF.")

V1_L400_guides_paths = write_ref_dropped_trees(file_list = V1_L400_guides_paths, pattern = "^REF.")
V2_L400_guides_paths = write_ref_dropped_trees(file_list = V2_L400_guides_paths, pattern = "^REF.")
V3_L400_guides_paths = write_ref_dropped_trees(file_list = V3_L400_guides_paths, pattern = "^REF.")
V4_L400_guides_paths = write_ref_dropped_trees(file_list = V4_L400_guides_paths, pattern = "^REF.")
V5_L400_guides_paths = write_ref_dropped_trees(file_list = V5_L400_guides_paths, pattern = "^REF.")
V6_L400_guides_paths = write_ref_dropped_trees(file_list = V6_L400_guides_paths, pattern = "^REF.")

V1_L500_guides_paths = write_ref_dropped_trees(file_list = V1_L500_guides_paths, pattern = "^REF.")
V2_L500_guides_paths = write_ref_dropped_trees(file_list = V2_L500_guides_paths, pattern = "^REF.")
V3_L500_guides_paths = write_ref_dropped_trees(file_list = V3_L500_guides_paths, pattern = "^REF.")
V4_L500_guides_paths = write_ref_dropped_trees(file_list = V4_L500_guides_paths, pattern = "^REF.")
V5_L500_guides_paths = write_ref_dropped_trees(file_list = V5_L500_guides_paths, pattern = "^REF.")

CORES = 10
vfull_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = VFull_guides_paths, cores = CORES)
vfull_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = VFull_control_paths, cores = CORES)

V1_L50_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V1_L50_guides_paths, cores = CORES)
V2_L50_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V2_L50_guides_paths, cores = CORES)
V3_L50_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V3_L50_guides_paths, cores = CORES)
V4_L50_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V4_L50_guides_paths, cores = CORES)
V5_L50_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V5_L50_guides_paths, cores = CORES)
V6_L50_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V6_L50_guides_paths, cores = CORES)
V7_L50_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V7_L50_guides_paths, cores = CORES)
V8_L50_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V8_L50_guides_paths, cores = CORES)

V1_L50_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V1_L50_paths, cores = CORES)
V2_L50_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V2_L50_paths, cores = CORES)
V3_L50_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V3_L50_paths, cores = CORES)
V4_L50_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V4_L50_paths, cores = CORES)
V5_L50_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V5_L50_paths, cores = CORES)
V6_L50_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V6_L50_paths, cores = CORES)
V7_L50_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V7_L50_paths, cores = CORES)
V8_L50_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V8_L50_paths, cores = CORES)

V1_L100_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V1_L100_guides_paths, cores = CORES)
V2_L100_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V2_L100_guides_paths, cores = CORES)
V3_L100_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V3_L100_guides_paths, cores = CORES)
V4_L100_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V4_L100_guides_paths, cores = CORES)
V5_L100_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V5_L100_guides_paths, cores = CORES)
V6_L100_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V6_L100_guides_paths, cores = CORES)
V7_L100_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V7_L100_guides_paths, cores = CORES)
V8_L100_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V8_L100_guides_paths, cores = CORES)

V1_L100_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V1_L100_paths, cores = CORES)
V2_L100_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V2_L100_paths, cores = CORES)
V3_L100_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V3_L100_paths, cores = CORES)
V4_L100_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V4_L100_paths, cores = CORES)
V5_L100_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V5_L100_paths, cores = CORES)
V6_L100_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V6_L100_paths, cores = CORES)
V7_L100_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V7_L100_paths, cores = CORES)
V8_L100_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V8_L100_paths, cores = CORES)

V1_L200_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V1_L200_guides_paths, cores = CORES)
V2_L200_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V2_L200_guides_paths, cores = CORES)
V3_L200_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V3_L200_guides_paths, cores = CORES)
V4_L200_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V4_L200_guides_paths, cores = CORES)
V5_L200_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V5_L200_guides_paths, cores = CORES)
V6_L200_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V6_L200_guides_paths, cores = CORES)
V7_L200_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V7_L200_guides_paths, cores = CORES)
V8_L200_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V8_L200_guides_paths, cores = CORES)

V1_L200_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V1_L200_paths, cores = CORES)
V2_L200_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V2_L200_paths, cores = CORES)
V3_L200_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V3_L200_paths, cores = CORES)
V4_L200_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V4_L200_paths, cores = CORES)
V5_L200_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V5_L200_paths, cores = CORES)
V6_L200_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V6_L200_paths, cores = CORES)
V7_L200_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V7_L200_paths, cores = CORES)
V8_L200_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V8_L200_paths, cores = CORES)

V1_L300_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V1_L300_guides_paths, cores = CORES)
V2_L300_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V2_L300_guides_paths, cores = CORES)
V3_L300_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V3_L300_guides_paths, cores = CORES)
V4_L300_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V4_L300_guides_paths, cores = CORES)
V5_L300_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V5_L300_guides_paths, cores = CORES)
V6_L300_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V6_L300_guides_paths, cores = CORES)
V7_L300_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V7_L300_guides_paths, cores = CORES)

V1_L300_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V1_L300_paths, cores = CORES)
V2_L300_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V2_L300_paths, cores = CORES)
V3_L300_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V3_L300_paths, cores = CORES)
V4_L300_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V4_L300_paths, cores = CORES)
V5_L300_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V5_L300_paths, cores = CORES)
V6_L300_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V6_L300_paths, cores = CORES)
V7_L300_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V7_L300_paths, cores = CORES)

V1_L400_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V1_L400_guides_paths, cores = CORES)
V2_L400_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V2_L400_guides_paths, cores = CORES)
V3_L400_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V3_L400_guides_paths, cores = CORES)
V4_L400_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V4_L400_guides_paths, cores = CORES)
V5_L400_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V5_L400_guides_paths, cores = CORES)
V6_L400_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V6_L400_guides_paths, cores = CORES)

V1_L400_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V1_L400_paths, cores = CORES)
V2_L400_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V2_L400_paths, cores = CORES)
V3_L400_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V3_L400_paths, cores = CORES)
V4_L400_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V4_L400_paths, cores = CORES)
V5_L400_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V5_L400_paths, cores = CORES)
V6_L400_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V6_L400_paths, cores = CORES)

V1_L500_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V1_L500_guides_paths, cores = CORES)
V2_L500_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V2_L500_guides_paths, cores = CORES)
V3_L500_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V3_L500_guides_paths, cores = CORES)
V4_L500_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V4_L500_guides_paths, cores = CORES)
V5_L500_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V5_L500_guides_paths, cores = CORES)

V1_L500_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V1_L500_paths, cores = CORES)
V2_L500_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V2_L500_paths, cores = CORES)
V3_L500_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V3_L500_paths, cores = CORES)
V4_L500_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V4_L500_paths, cores = CORES)
V5_L500_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V5_L500_paths, cores = CORES)

bp_quartet =
  data.frame("VFull_control" = vfull_quartet,
             "VFull_guides" = vfull_guides_quartet,
             "V1_L50_quartet" = V1_L50_quartet,
             "V1_L100_quartet" = V1_L100_quartet,
             "V1_L200_quartet" = V1_L200_quartet,
             "V1_L300_quartet" = V1_L300_quartet,
             "V1_L400_quartet" = V1_L400_quartet,
             "V1_L500_quartet" = V1_L500_quartet,
             "V1_L50_guides_quartet" = V1_L50_guides_quartet,
             "V1_L100_guides_quartet" = V1_L100_guides_quartet,
             "V1_L200_guides_quartet" = V1_L200_guides_quartet,
             "V1_L300_guides_quartet" = V1_L300_guides_quartet,
             "V1_L400_guides_quartet" = V1_L400_guides_quartet,
             "V1_L500_guides_quartet" = V1_L500_guides_quartet,
             "V2_L50_quartet" = V2_L50_quartet,
             "V2_L100_quartet" = V2_L100_quartet,
             "V2_L200_quartet" = V2_L200_quartet,
             "V2_L300_quartet" = V2_L300_quartet,
             "V2_L400_quartet" = V2_L400_quartet,
             "V2_L500_quartet" = V2_L500_quartet,
             "V2_L50_guides_quartet" = V2_L50_guides_quartet,
             "V2_L100_guides_quartet" = V2_L100_guides_quartet,
             "V2_L200_guides_quartet" = V2_L200_guides_quartet,
             "V2_L300_guides_quartet" = V2_L300_guides_quartet,
             "V2_L400_guides_quartet" = V2_L400_guides_quartet,
             "V2_L500_guides_quartet" = V2_L500_guides_quartet,
             "V3_L50_quartet" = V3_L50_quartet,
             "V3_L100_quartet" = V3_L100_quartet,
             "V3_L200_quartet" = V3_L200_quartet,
             "V3_L300_quartet" = V3_L300_quartet,
             "V3_L400_quartet" = V3_L400_quartet,
             "V3_L500_quartet" = V3_L500_quartet,
             "V3_L50_guides_quartet" = V3_L50_guides_quartet,
             "V3_L100_guides_quartet" = V3_L100_guides_quartet,
             "V3_L200_guides_quartet" = V3_L200_guides_quartet,
             "V3_L300_guides_quartet" = V3_L300_guides_quartet,
             "V3_L400_guides_quartet" = V3_L400_guides_quartet,
             "V3_L500_guides_quartet" = V3_L500_guides_quartet,
             "V4_L50_quartet" = V4_L50_quartet,
             "V4_L100_quartet" = V4_L100_quartet,
             "V4_L200_quartet" = V4_L200_quartet,
             "V4_L300_quartet" = V4_L300_quartet,
             "V4_L400_quartet" = V4_L400_quartet,
             "V4_L500_quartet" = V4_L500_quartet,
             "V4_L50_guides_quartet" = V4_L50_guides_quartet,
             "V4_L100_guides_quartet" = V4_L100_guides_quartet,
             "V4_L200_guides_quartet" = V4_L200_guides_quartet,
             "V4_L300_guides_quartet" = V4_L300_guides_quartet,
             "V4_L400_guides_quartet" = V4_L400_guides_quartet,
             "V4_L500_guides_quartet" = V4_L500_guides_quartet,
             "V5_L50_quartet" = V5_L50_quartet,
             "V5_L100_quartet" = V5_L100_quartet,
             "V5_L200_quartet" = V5_L200_quartet,
             "V5_L300_quartet" = V5_L300_quartet,
             "V5_L400_quartet" = V5_L400_quartet,
             "V5_L500_quartet" = V5_L500_quartet,
             "V5_L50_guides_quartet" = V5_L50_guides_quartet,
             "V5_L100_guides_quartet" = V5_L100_guides_quartet,
             "V5_L200_guides_quartet" = V5_L200_guides_quartet,
             "V5_L300_guides_quartet" = V5_L300_guides_quartet,
             "V5_L400_guides_quartet" = V5_L400_guides_quartet,
             "V5_L500_guides_quartet" = V5_L500_guides_quartet,
             "V6_L50_quartet" = V6_L50_quartet,
             "V6_L100_quartet" = V6_L100_quartet,
             "V6_L200_quartet" = V6_L200_quartet,
             "V6_L300_quartet" = V6_L300_quartet,
             "V6_L400_quartet" = V6_L400_quartet,
             "V6_L50_guides_quartet" = V6_L50_guides_quartet,
             "V6_L100_guides_quartet" = V6_L100_guides_quartet,
             "V6_L200_guides_quartet" = V6_L200_guides_quartet,
             "V6_L300_guides_quartet" = V6_L300_guides_quartet,
             "V6_L400_guides_quartet" = V6_L400_guides_quartet,
             "V7_L50_quartet" = V7_L50_quartet,
             "V7_L100_quartet" = V7_L100_quartet,
             "V7_L200_quartet" = V7_L200_quartet,
             "V7_L300_quartet" = V7_L300_quartet,
             "V7_L50_guides_quartet" = V7_L50_guides_quartet,
             "V7_L100_guides_quartet" = V7_L100_guides_quartet,
             "V7_L200_guides_quartet" = V7_L200_guides_quartet,
             "V7_L300_guides_quartet" = V7_L300_guides_quartet,
             "V8_L50_quartet" = V8_L50_quartet,
             "V8_L100_quartet" = V8_L100_quartet,
             "V8_L200_quartet" = V8_L200_quartet,
             "V8_L50_guides_quartet" = V8_L50_guides_quartet,
             "V8_L100_guides_quartet" = V8_L100_guides_quartet,
             "V8_L200_guides_quartet" = V8_L200_guides_quartet)
bp_quartet
setwd(dir.sim.out)
#saveRDS(file = "bp_quartet_random_pairs.rds", object = bp_quartet, compress = TRUE)
bp_quartet = as_tibble(readRDS("bp_quartet_random_pairs.rds"))
bp_quartet
QUARTET_MAX = choose(TREE_SIZE, 4)
#saveRDS(file = "bp_quartet_normalized_random_pairs.rds", object = bp_quartet/QUARTET_MAX, compress = TRUE)
bp_quartet_normalized = as_tibble(readRDS("bp_quartet_normalized.rds"))
bp_quartet_normalized
dplyr::as_tibble(readRDS("bp_quartet_normalized_random_pairs.rds"))
dplyr::as_tibble(readRDS("bp_quartet_random_pairs.rds"))

################################################################################
# 1F: Determine how guides x length impact topological reconstruction -- PATH
################################################################################
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_bp_path_random_pairs.R' -r bp_path_random_pairs -p 50 -q biomed
dplyr::as_tibble(readRDS("bp_path_random_pairs.rds"))


################################################################################
# 1G: Determine how guides x length impact topological reconstruction -- NYE
################################################################################
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_bp_nye_random_pairs.R' -r bp_nye_random_pairs -p 100 -q sharpton
dplyr::as_tibble(readRDS("bp_nye_random_pairs.rds"))

################################################################################
# 1H: Determine how guides x length impact topological reconstruction -- MCD
################################################################################
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_bp_mcd_random_pairs.R' -r bp_mcd_random_pairs -p 100 -q sharpton
dplyr::as_tibble(readRDS("bp_mcd_random_pairs.rds"))

################################################################################
# 1H: Determine how guides x length impact topological reconstruction -- Nanostructure
################################################################################

hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_bp_nanostructure_random_pairs_2.R' -r  bp_nanoclusters_random_pairs_2 -p 30 -q all.q
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_bp_nanostructure_random_pairs_3.R' -r  bp_nanoclusters_random_pairs_3 -p 30 -q all.q
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_bp_nanostructure_random_pairs_4.R' -r  bp_nanoclusters_random_pairs_4 -p 100 -q sharpton
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_bp_nanostructure_random_pairs_5.R' -r  bp_nanoclusters_random_pairs_5 -p 100 -q sharpton
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_bp_nanostructure_random_pairs_6.R' -r  bp_nanoclusters_random_pairs_6 -p 100 -q sharpton
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_bp_nanostructure_random_pairs_7.R' -r  bp_nanoclusters_random_pairs_7 -p 100 -q sharpton
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_bp_nanostructure_random_pairs_8.R' -r  bp_nanoclusters_random_pairs_8 -p 100 -q sharpton
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_bp_nanostructure_random_pairs_9.R' -r  bp_nanoclusters_random_pairs_9 -p 30 -q all.q
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_bp_nanostructure_random_pairs_10.R' -r  bp_nanoclusters_random_pairs_10 -p 30 -q all.q

as_tibble(readRDS("nanostructure_bp_2_random_pairs.rds"))
as_tibble(readRDS("nanostructure_bp_3_random_pairs.rds"))
as_tibble(readRDS("nanostructure_bp_4_random_pairs.rds"))
as_tibble(readRDS("nanostructure_bp_5_random_pairs.rds"))
as_tibble(readRDS("nanostructure_bp_6_random_pairs.rds"))
as_tibble(readRDS("nanostructure_bp_7_random_pairs.rds"))
as_tibble(readRDS("nanostructure_bp_8_random_pairs.rds"))
as_tibble(readRDS("nanostructure_bp_9_random_pairs.rds"))
as_tibble(readRDS("nanostructure_bp_10_random_pairs.rds"))

#######################################
# 1J.Determine MAST
#######################################

hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_bp_mast.R' -r  bp_mast -p 30 -q all.q
as_tibble(readRDS("bp_mast.rds"))

#######################################
# Clean up
#######################################
rm(list = ls())


################################################################################
# SIM 3: Combination of overlapping HVRs 
################################################################################

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
dir.sim.out = file.path(HOME, "2025_hvr_guide_phylogenetic_integration/out/sim_overlapping/")

setwd(HOME_TMP)
dir.simulation = file.path(HOME_TMP, "length_simulations/")
dir.mix = file.path(HOME_TMP, "mix_experiment/")

setwd(dir.scripts)
source("sim_analysis_functions_V2.R")

################################################################################
# SET CONSTANTS
################################################################################
TREE_SIZE = 6846

# Check that all the tree files are present. 
directories <- c(list.dirs(dir.mix, full.names = TRUE, recursive = FALSE),
                 list.dirs(dir.simulation, full.names = TRUE, recursive = FALSE))
names(directories) = basename(directories)
basename(directories)

################################################################################
# 0A: READ IN FILES
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
# 0B: Check sequence length and tree sizes are what we expect.
################################################################################
if(all(c(check_tree_size(tree_list = VFull, expected_tree_size = TREE_SIZE),
         check_tree_size(tree_list = VFull_guides, expected_tree_size = TREE_SIZE),
         check_tree_size(tree_list = VFull_control, expected_tree_size = TREE_SIZE),
         check_tree_size(tree_list = V4, expected_tree_size = TREE_SIZE),
         check_tree_size(tree_list = V3V4, expected_tree_size = TREE_SIZE),
         check_tree_size(tree_list = V4_V3V4_mix, expected_tree_size = TREE_SIZE),
         check_tree_size(tree_list = V4_guides, expected_tree_size = TREE_SIZE),
         check_tree_size(tree_list = V3V4_guides, expected_tree_size = TREE_SIZE),
         check_tree_size(tree_list = V4_V3V4_mix_guides, expected_tree_size = TREE_SIZE)))){
  print("PASS: All of overlapping simulations are same size")
}else{
  print("FAIL: Not all overlapping simulations are same size")
}


################################################################################
# RF distance sim 3
################################################################################
CORES = 50
VFull_control = calculate_rf_distance_random_pairs(VFull, VFull_control, cores = CORES)
VFull_guides = calculate_rf_distance_random_pairs(VFull, VFull_guides, cores = CORES)
V4_rf = calculate_rf_distance_random_pairs(VFull, V4, cores = CORES)
V3V4_rf = calculate_rf_distance_random_pairs(VFull, V3V4, cores = CORES)
V4_guides_rf = calculate_rf_distance_random_pairs(VFull, V4_guides, cores = CORES)
V3V4_guides_rf = calculate_rf_distance_random_pairs(VFull, V3V4_guides, cores = CORES)
V4_V3V4_mix_rf = calculate_rf_distance_random_pairs(VFull, V4_V3V4_mix, cores = CORES)
V4_V3V4_mix_guides_rf = calculate_rf_distance_random_pairs(VFull, V4_V3V4_mix_guides, cores = CORES)

ov_rf = data.frame("VFull" = VFull_control,
                   "VFull_guides" = VFull_guides,
                   "V4_rf" = V4_rf,
                   "V3V4_rf" = V3V4_rf,
                   "V4_guides_rf" = V4_guides_rf,
                   "V3V4_guides_rf" = V3V4_guides_rf,
                   "V4_V3V4_mix_rf" = V4_V3V4_mix_rf,
                   "V4_V3V4_mix_guides_rf" = V4_V3V4_mix_guides_rf)

setwd(dir.sim.out)
#saveRDS(file = "ov_rf_random_pairs.rds", object = ov_rf)
ov_rf = readRDS("ov_rf_random_pairs.rds")
ov_rf
as_tibble(readRDS("ov_rf_random_pairs.rds"))

################################################################################
# JRF distance Sim 3
################################################################################
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_ov_jrf.R' -r ov_jrf -p 30 -q all.q
as_tibble(readRDS("ov_jrf_random_pairs.rds"))

################################################################################
# Nye distance sim 3
################################################################################
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_ov_nye.R' -r ov_nye -p 30 -q all.q
as_tibble(readRDS("ov_nye_random_pairs.rds"))

################################################################################
# mcd distance sim 3
################################################################################
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_ov_msd.R' -r ov_msd -p 30 -q all.q
as_tibble(readRDS("ov_msd_random_pairs.rds"))

################################################################################
# 1E: pid distance sim 3
################################################################################
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_ov_pid.R' -r ov_pid -p 30 -q all.q
as_tibble(readRDS("ov_pid_random_pairs.rds"))

################################################################################
# 1F: mcd distance sim 3
################################################################################
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_ov_mcd.R' -r ov_mcd -p 30 -q all.q
as_tibble(readRDS("ov_mcd_random_pairs.rds"))


################################################################################
# 1F: path distance sim 3
################################################################################
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_ov_path.R' -r ov_path -p 30 -q all.q
as_tibble(readRDS("ov_path_random_pairs.rds"))


################################################################################
# 1G: quartet distance sim 3
################################################################################
VFull_paths = file.path(dir.hvrs["VFull"], file_list = list.files(path = dir.hvrs["VFull"], pattern = ".tre"))[1:N]
VFull_guides_paths = file.path(dir.hvrs["VFull_guides"], file_list = list.files(path = dir.hvrs["VFull_guides"], pattern = "align.tre"), drop = TRUE)[1:N]
VFull_control_paths = file.path(dir.hvrs["VFull_control"], file_list = list.files(path = dir.hvrs["VFull_control"], pattern = "align.tre"))[1:N]

V4_paths = file.path(dir.hvrs["V4"], file_list = list.files(path = dir.hvrs["V4"], pattern = ".tre"))[1:N]
V3V4_paths = file.path(dir.hvrs["V3V4"], file_list = list.files(path = dir.hvrs["V3V4"], pattern = ".tre"))[1:N]
V4_V3V4_mix_paths = file.path(dir.hvrs["V4_V3V4_mix"], file_list = list.files(path = dir.hvrs["V4_V3V4_mix"], pattern = ".tre"))[1:N]

V4_guides_paths = file.path(dir.hvrs["V4_guides"], file_list = list.files(path = dir.hvrs["V4_guides"], pattern = ".tre"))[1:N]
V3V4_guides_paths = file.path(dir.hvrs["V3V4_guides"], file_list = list.files(path = dir.hvrs["V3V4_guides"], pattern = ".tre"))[1:N]
V4_V3V4_mix_guides_paths = file.path(dir.hvrs["V4_V3V4_mix_guides"], file_list = list.files(path = dir.hvrs["V4_V3V4_mix_guides"], pattern = ".tre"))[1:N]

V4_guides_paths = write_ref_dropped_trees(file_list = V4_guides_paths, pattern = "^REF.")
V3V4_guides_paths = write_ref_dropped_trees(file_list = V3V4_guides_paths, pattern = "^REF.")
V4_V3V4_mix_guides_paths = write_ref_dropped_trees(file_list = V4_V3V4_mix_guides_paths, pattern = "^REF.")


CORES = 10
vfull_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = VFull_control_paths, cores = CORES)
vfull_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = VFull_guides_paths, cores = CORES)

V4_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V4_paths, cores = CORES)
V3V4_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V3V4_paths, cores = CORES)
V4_V3V4_mix_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V4_V3V4_mix_paths, cores = CORES)

V4_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V4_guides_paths, cores = CORES)
V3V4_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V3V4_guides_paths, cores = CORES)
V4_V3V4_mix_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V4_V3V4_mix_guides_paths, cores = CORES)

ov_quartet =
  data.frame("VFull_control" = vfull_quartet,
             "VFull_guides" = vfull_guides_quartet,
             "V4_quartet" = V4_quartet,
             "V3V4_quartet" = V3V4_quartet,
             "V4_V3V4_mix_quartet" = V4_V3V4_mix_quartet,
             "V4_guides_quartet" = V4_guides_quartet,
             "V3V4_guides_quartet" = V3V4_guides_quartet,
             "V4_V3V4_mix_guides_quartet" = V4_V3V4_mix_guides_quartet) 
ov_quartet
setwd(dir.sim.out)
saveRDS(file = "ov_quartet_random_pairs.rds", object = ov_quartet, compress = TRUE)
ov_quartet = as_tibble(readRDS("ov_quartet_random_pairs.rds"))
ov_quartet
QUARTET_MAX = choose(TREE_SIZE, 4)
saveRDS(file = "ov_quartet_normalized_random_pairs.rds", object = ov_quartet/QUARTET_MAX, compress = TRUE)
ov_quartet_normalized = as_tibble(readRDS("ov_quartet_normalized_random_pairs.rds"))
ov_quartet_normalized
as_tibble(readRDS("ov_quartet_random_pairs.rds"))
as_tibble(readRDS("ov_quartet_normalized_random_pairs.rds"))



################################################################################
# 1H: nanosimilarity distance sim 3
################################################################################

hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_ov_nanostructure_random_pairs_2.R' -r  ov_nanoclusters_random_pairs_2 -p 50 -q sharpton
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_ov_nanostructure_random_pairs_3.R' -r  ov_nanoclusters_random_pairs_3 -p 100 -q sharpton
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_ov_nanostructure_random_pairs_4.R' -r  ov_nanoclusters_random_pairs_4 -p 100 -q sharpton
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_ov_nanostructure_random_pairs_5.R' -r  ov_nanoclusters_random_pairs_5 -p 100 -q sharpton
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_ov_nanostructure_random_pairs_6.R' -r  ov_nanoclusters_random_pairs_6 -p 100 -q sharpton
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_ov_nanostructure_random_pairs_7.R' -r  ov_nanoclusters_random_pairs_7 -p 100 -q sharpton
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_ov_nanostructure_random_pairs_8.R' -r  ov_nanoclusters_random_pairs_8 -p 100 -q sharpton
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_ov_nanostructure_random_pairs_9.R' -r  ov_nanoclusters_random_pairs_9 -p 100 -q sharpton
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_ov_nanostructure_random_pairs_10.R' -r ov_nanoclusters_random_pairs_10 -p 100 -q sharpton

as_tibble(readRDS("ov_nanoclusters_random_pairs_2.rds"))
as_tibble(readRDS("ov_nanoclusters_random_pairs_3.rds"))
as_tibble(readRDS("ov_nanoclusters_random_pairs_4.rds"))
as_tibble(readRDS("ov_nanoclusters_random_pairs_5.rds"))
as_tibble(readRDS("ov_nanoclusters_random_pairs_6.rds"))
as_tibble(readRDS("ov_nanoclusters_random_pairs_7.rds"))
as_tibble(readRDS("ov_nanoclusters_random_pairs_8.rds"))
as_tibble(readRDS("ov_nanoclusters_random_pairs_9.rds"))
as_tibble(readRDS("ov_nanoclusters_random_pairs_10.rds"))

#######################################
# 1J.Determine MAST sim 3
#######################################

hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_ov_mast.R' -r  ov_mast -p 30 -q all.q

as_tibble(readRDS("ov_mast_random_pairs.rds"))

#######################################
# Clean up
#######################################
rm(list = ls())


################################################################################
# SIM 4: Combination of nonoverlapping HVRs 
################################################################################

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
# SET CONSTANTS
################################################################################
TREE_SIZE = 6846

# Check that all the tree files are present. 
directories <- c(list.dirs(dir.mix, full.names = TRUE, recursive = FALSE),
                 list.dirs(dir.simulation, full.names = TRUE, recursive = FALSE))
names(directories) = basename(directories)
basename(directories)

################################################################################
# 0A: READ IN FILES
################################################################################
N = 100 # The number of trees to look at per simulation
setwd(dir.mix)

dir.mix = directories[c("V1V2", "V1V2_guides", 
                        "V3V4", "V3V4_guides", 
                        "V5V7", "V5V7_guides", 
                        "V1V2_V3V4_V5V7_mix", "V1V2_V3V4_V5V7_mix_guides")]
dir.mix

setwd(dir.simulation)
dir.control = directories[c("VFull", "VFull_control", "VFull_guides")]
dir.control

dir.hvrs = c(dir.mix, dir.control)


setwd(dir.mix)
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

################################################################################
# 0B: Check sequence length and tree sizes are what we expect.
################################################################################
if(all(c(check_tree_size(tree_list = VFull, expected_tree_size = TREE_SIZE),
         check_tree_size(tree_list = VFull_guides, expected_tree_size = TREE_SIZE),
         check_tree_size(tree_list = VFull_control, expected_tree_size = TREE_SIZE),
         check_tree_size(tree_list = V1V2, expected_tree_size = TREE_SIZE),
         check_tree_size(tree_list = V3V4, expected_tree_size = TREE_SIZE),
         check_tree_size(tree_list = V5V7, expected_tree_size = TREE_SIZE),
         check_tree_size(tree_list = V1V2_guides, expected_tree_size = TREE_SIZE),
         check_tree_size(tree_list = V3V4_guides, expected_tree_size = TREE_SIZE),
         check_tree_size(tree_list = V5V7_guides, expected_tree_size = TREE_SIZE),
         check_tree_size(tree_list = V1V2_V3V4_V5V7_mix, expected_tree_size = TREE_SIZE),
         check_tree_size(tree_list = V1V2_V3V4_V5V7_mix_guides, expected_tree_size = TREE_SIZE)))){
  print("PASS: All of overlapping simulations are same size")
}else{
  print("FAIL: Not all overlapping simulations are same size")
}

################################################################################
# 1A RF
################################################################################
CORES = 50
VFull_control_rf = calculate_rf_distance_random_pairs(VFull, VFull_control, cores = CORES)
VFull_guides_rf = calculate_rf_distance_random_pairs(VFull, VFull_guides, cores = CORES)
V1V2_rf = calculate_rf_distance_random_pairs(VFull, V1V2, cores = CORES)
V3V4_rf = calculate_rf_distance_random_pairs(VFull, V3V4, cores = CORES)
V5V7_rf = calculate_rf_distance_random_pairs(VFull, V5V7, cores = CORES)

V1V2_guides_rf = calculate_rf_distance_random_pairs(VFull, V1V2_guides, cores = CORES)
V3V4_guides_rf = calculate_rf_distance_random_pairs(VFull, V3V4_guides, cores = CORES)
V5V7_guides_rf = calculate_rf_distance_random_pairs(VFull, V5V7_guides, cores = CORES)

V1V2_V3V4_V5V7_mix_rf = calculate_rf_distance_random_pairs(VFull, V1V2_V3V4_V5V7_mix, cores = CORES)
V1V2_V3V4_V5V7_mix_guides_rf = calculate_rf_distance_random_pairs(VFull, V1V2_V3V4_V5V7_mix_guides, cores = CORES)

nonov_rf = data.frame("VFull" = VFull_control_rf,
                      "VFull_guides" = VFull_guides_rf,
                      "V1V2" = V1V2_rf,
                      "V3V4" = V3V4_rf,
                      "V5V7" = V5V7_rf,
                      
                      "V1V2_guides" = V1V2_guides_rf,
                      "V3V4_guides" = V3V4_guides_rf,
                      "V5V7_guides" = V5V7_guides_rf,
                      "V1V2_V3V4_V5V7_mix" = V1V2_V3V4_V5V7_mix_rf,
                      "V1V2_V3V4_V5V7_mix_guides" = V1V2_V3V4_V5V7_mix_guides_rf)
as_tibble(nonov_rf)
setwd(dir.sim.out)
saveRDS(file = "nonov_rf_random_pairs.rds", object = nonov_rf)
nonov_rf = readRDS("nonov_rf_random_pairs.rds")
nonov_rf

as_tibble(readRDS("nonov_rf_random_pairs.rds"))


################################################################################
# 1B: Determine how guides impact topological reconstruction nonoverlapping -- JRF
################################################################################
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_nonov_jrf.R' -r nonov_jrf -p 100 -q sharpton 
as_tibble(readRDS("nonov_jrf_random_pairs.rds"))

################################################################################
# 1C: Determine how guidesimpact topological reconstruction nonoverlapping-- NYE
################################################################################

hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_nonov_nye.R' -r nonov_nye -p 100 -q sharpton
as_tibble(readRDS("nonov_nye_random_pairs.rds"))

################################################################################
# 1D: Determine how guides x length impact topological reconstruction -- MSD
################################################################################
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_nonov_msd.R' -r nonov_msd -p 30 -q all.q
as_tibble(readRDS("nonov_msd_random_pairs.rds"))

################################################################################
# 1E: Determine how guides x length impact topological reconstruction -- PID
################################################################################
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_nonov_pid.R' -r nonov_pid -p 100 -q sharpton
as_tibble(readRDS("nonov_pid_random_pairs.rds"))

################################################################################
# 1F: Determine how guides x length impact topological reconstruction -- MCD
################################################################################
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_nonov_mcd.R' -r nonov_mcd -p 30 -q all.q
as_tibble(readRDS("nonov_mcd_random_pairs.rds"))


################################################################################
# 1F: Determine how guides x length impact topological reconstruction -- PATH
################################################################################
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_nonov_path.R' -r nonov_path -p 30 -q biomed
as_tibble(readRDS("nonov_path_random_pairs.rds"))

################################################################################
# 1G: Determine how guides x length impact topological reconstruction -- QD
################################################################################
VFull_paths = file.path(dir.hvrs["VFull"], file_list = list.files(path = dir.hvrs["VFull"], pattern = ".tre"))[1:N]
VFull_guides_paths = file.path(dir.hvrs["VFull_guides"], file_list = list.files(path = dir.hvrs["VFull_guides"], pattern = "align.tre"))[1:N]
VFull_control_paths = file.path(dir.hvrs["VFull_control"], file_list = list.files(path = dir.hvrs["VFull_control"], pattern = "align.tre"))[1:N]

V1V2_paths = file.path(dir.hvrs["V1V2"], file_list = list.files(path = dir.hvrs["V1V2"], pattern = ".tre"))[1:N]
V3V4_paths = file.path(dir.hvrs["V3V4"], file_list = list.files(path = dir.hvrs["V3V4"], pattern = ".tre"))[1:N]
V5V7_paths = file.path(dir.hvrs["V5V7"], file_list = list.files(path = dir.hvrs["V5V7"], pattern = ".tre"))[1:N]
V1V2_V3V4_V5V7_mix_paths = file.path(dir.hvrs["V1V2_V3V4_V5V7_mix"], file_list = list.files(path = dir.hvrs["V1V2_V3V4_V5V7_mix"], pattern = ".tre"))[1:N]

V1V2_guides_paths = file.path(dir.hvrs["V1V2_guides"], file_list = list.files(path = dir.hvrs["V1V2_guides"], pattern = ".tre"))[1:N]
V3V4_guides_paths = file.path(dir.hvrs["V3V4_guides"], file_list = list.files(path = dir.hvrs["V3V4_guides"], pattern = ".tre"))[1:N]
V5V7_guides_paths = file.path(dir.hvrs["V5V7"], file_list = list.files(path = dir.hvrs["V5V7"], pattern = ".tre"))[1:N]
V1V2_V3V4_V5V7_mix_guides_paths = file.path(dir.hvrs["V1V2_V3V4_V5V7_mix_guides"], file_list = list.files(path = dir.hvrs["V1V2_V3V4_V5V7_mix_guides"], pattern = ".tre"))[1:N]

V1V2_guides_paths =  write_ref_dropped_trees(file_list = V1V2_guides_paths, pattern = "^REF.")
V3V4_guides_paths =  write_ref_dropped_trees(file_list = V3V4_guides_paths, pattern = "^REF.")
V5V7_guides_paths =  write_ref_dropped_trees(file_list = V5V7_guides_paths, pattern = "^REF.")
VFull_guides_paths =  write_ref_dropped_trees(file_list = VFull_guides_paths, pattern = "^REF.")
V1V2_V3V4_V5V7_mix_guides_paths =  write_ref_dropped_trees(file_list = V1V2_V3V4_V5V7_mix_guides_paths, pattern = "^REF.")

CORES = 10
vfull_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = VFull_control_paths, cores = CORES)
vfull_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = VFull_guides_paths, cores = CORES)

V1V2_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V1V2_paths, cores = CORES)
V3V4_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V3V4_paths, cores = CORES)
V5V7_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V5V7_paths, cores = CORES)
V1V2_V3V4_V5V7_mix_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V1V2_V3V4_V5V7_mix_paths, cores = CORES)

V1V2_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V1V2_guides_paths, cores = CORES)
V3V4_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V3V4_guides_paths, cores = CORES)
V5V7_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V3V4_guides_paths, cores = CORES)
V1V2_V3V4_V5V7_mix_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V1V2_V3V4_V5V7_mix_guides_paths, cores = CORES)


nonov_quartet =
  data.frame("VFull_control" = vfull_quartet,
             "VFull_guides" = vfull_guides_quartet,
             "V1V2_quartet" = V1V2_quartet,
             "V3V4_quartet" = V3V4_quartet,
             "V5V7_quartet" = V5V7_quartet,
             "V1V2_V3V4_V5V7_mix_quartet" = V1V2_V3V4_V5V7_mix_quartet,
             
             "V1V2_guides_quartet" = V1V2_guides_quartet,
             "V3V4_guides_quartet" = V3V4_guides_quartet,
             "V5V7_guides_quartet" = V5V7_guides_quartet,
             "V1V2_V3V4_V5V7_mix_guides_quartet" = V1V2_V3V4_V5V7_mix_guides_quartet) 
nonov_quartet
setwd(dir.sim.out)
saveRDS(file = "nonov_quartet_random_pairs.rds", object = nonov_quartet, compress = TRUE)
nonov_quartet = as_tibble(readRDS("nonov_quartet_random_pairs.rds"))
nonov_quartet
QUARTET_MAX = choose(TREE_SIZE, 4)
saveRDS(file = "nonov_quartet_normalized_random_pairs.rds", object = nonov_quartet/QUARTET_MAX, compress = TRUE)
nonov_quartet_normalized = as_tibble(readRDS("nonov_quartet_normalized_random_pairs.rds"))
nonov_quartet_normalized

as_tibble(readRDS("nonov_quartet_normalized_random_pairs.rds"))
as_tibble(readRDS("nonov_quartet_random_pairs.rds"))


################################################################################
# 1H: Determine how guides x length impact topological reconstruction -- Nanostructure
################################################################################
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_nonov_nanostructure_random_pairs_2.R' -r  nonov_nanoclusters_random_pairs_2 -p 100 -q sharpton
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_nonov_nanostructure_random_pairs_3.R' -r  nonov_nanoclusters_random_pairs_3 -p 100 -q sharpton
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_nonov_nanostructure_random_pairs_4.R' -r  nonov_nanoclusters_random_pairs_4 -p 100 -q sharpton
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_nonov_nanostructure_random_pairs_5.R' -r  nonov_nanoclusters_random_pairs_5 -p 100 -q sharpton
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_nonov_nanostructure_random_pairs_6.R' -r  nonov_nanoclusters_random_pairs_6 -p 100 -q sharpton
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_nonov_nanostructure_random_pairs_7.R' -r  nonov_nanoclusters_random_pairs_7 -p 100 -q sharpton
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_nonov_nanostructure_random_pairs_8.R' -r  nonov_nanoclusters_random_pairs_8 -p 100 -q sharpton
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_nonov_nanostructure_random_pairs_9.R' -r  nonov_nanoclusters_random_pairs_9 -p 100 -q sharpton
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_nonov_nanostructure_random_pairs_10.R' -r  nonov_nanoclusters_random_pairs_10 -p 100 -q sharpton
as_tibble(readRDS("nonov_nanostructure_2_random_pairs.rds"))
as_tibble(readRDS("nonov_nanostructure_3_random_pairs.rds"))
as_tibble(readRDS("nonov_nanostructure_4_random_pairs.rds"))
as_tibble(readRDS("nonov_nanostructure_5_random_pairs.rds"))
as_tibble(readRDS("nonov_nanostructure_6_random_pairs.rds"))
as_tibble(readRDS("nonov_nanostructure_7_random_pairs.rds"))
as_tibble(readRDS("nonov_nanostructure_8_random_pairs.rds"))
as_tibble(readRDS("nonov_nanostructure_9_random_pairs.rds"))
as_tibble(readRDS("nonov_nanostructure_10_random_pairs.rds"))

#######################################
# 1J.Determine MAST
#######################################

hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_nonov_mast.R' -r  nonov_mast -p 30 -q biomed
as_tibble(readRDS("nonov_mast_random_pairs.rds"))


#######################################
# Clean up
#######################################
rm(list = ls())



################################################################################
# SIM 5: ND Mix
################################################################################

################################################################################
# SET DIRECTORIES
################################################################################
# Set Home
HOME = "/nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration"
HOME_TMP = "/nfs3/Sharpton_Lab/tmp/projects/arnoldh/2025_hvr_guide_phylogenetic_integration/"
HOME_ND = "/nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2022NeurodegenerativeDisease/2022NeurodegenerativeDisease/metadata/metaanalysis_study_summary/"


# Make directories
setwd(HOME)

dir.scripts = file.path(HOME, "2025_hvr_guide_phylogenetic_integration/scripts/")
dir.ref = file.path(HOME, "reference_data/")
dir.out = file.path(HOME, "2025_hvr_guide_phylogenetic_integration/data_visulization/")
dir.sim.out = file.path(HOME, "2025_hvr_guide_phylogenetic_integration/out/sim_nd_mix/")

setwd(HOME_TMP)
dir.simulation = file.path(HOME_TMP, "length_simulations/")
dir.mix = file.path(HOME_TMP, "mix_experiment/")

setwd(dir.scripts)
source("sim_analysis_functions_V2.R")

################################################################################
# SET CONSTANTS
################################################################################
TREE_SIZE = 6846

# Check that all the tree files are present. 
directories <- c(list.dirs(dir.mix, full.names = TRUE, recursive = FALSE),
                 list.dirs(dir.simulation, full.names = TRUE, recursive = FALSE))
names(directories) = basename(directories)
basename(directories)

################################################################################
# 0A: READ IN FILES
################################################################################
N = 100 # The number of trees to look at per simulation


# See what HVRs are in our ND study
setwd(HOME_ND)
nd_meta = readxl::read_excel(path = "study_summary26.xlsx", sheet = 2, skip = 2)
table(nd_meta$HVR)

setwd(dir.mix)
dir.mix = directories[c("V1V2", "V1V2_guides", 
                        "V1V3", "V1V3_guides", 
                        "V3V4", "V3V4_guides", 
                        "V3V5", "V3V5_guides",
                        "V4", "V4_guides",
                        "V4V5", "V4V5_guides", 
                        "nd_mix", "nd_mix_guides")]
dir.mix

setwd(dir.simulation)
dir.control = directories[c("VFull", "VFull_control", "VFull_guides")]
dir.control
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
# 0B: Check sequence length and tree sizes are what we expect.
################################################################################
if(all(c(check_tree_size(tree_list = VFull, expected_tree_size = TREE_SIZE),
         check_tree_size(tree_list = VFull_guides, expected_tree_size = TREE_SIZE),
         check_tree_size(tree_list = VFull_control, expected_tree_size = TREE_SIZE),
         
         check_tree_size(tree_list = V1V2, expected_tree_size = TREE_SIZE),
         check_tree_size(tree_list = V1V3, expected_tree_size = TREE_SIZE),
         check_tree_size(tree_list = V3V4, expected_tree_size = TREE_SIZE),
         check_tree_size(tree_list = V3V5, expected_tree_size = TREE_SIZE),
         check_tree_size(tree_list = V4, expected_tree_size = TREE_SIZE),
         check_tree_size(tree_list = V4V5, expected_tree_size = TREE_SIZE),
         check_tree_size(tree_list = nd_mix, expected_tree_size = TREE_SIZE),
         
         check_tree_size(tree_list = V1V2_guides, expected_tree_size = TREE_SIZE),
         check_tree_size(tree_list = V1V3_guides, expected_tree_size = TREE_SIZE),
         check_tree_size(tree_list = V3V4_guides, expected_tree_size = TREE_SIZE),
         check_tree_size(tree_list = V3V5_guides, expected_tree_size = TREE_SIZE),
         check_tree_size(tree_list = V4_guides, expected_tree_size = TREE_SIZE),
         check_tree_size(tree_list = V4V5_guides, expected_tree_size = TREE_SIZE),
         check_tree_size(tree_list = nd_mix_guides, expected_tree_size = TREE_SIZE)))){
  print("PASS: All of overlapping simulations are same size")
}else{
  print("FAIL: Not all overlapping simulations are same size")
}

################################################################################
# 1A: Determine how guides x length impact topological reconstruction -- RF
################################################################################
CORES = 50
VFull_control = calculate_rf_distance_random_pairs(VFull, VFull_control, cores = CORES)
VFull_guides = calculate_rf_distance_random_pairs(VFull, VFull_guides, cores = CORES)

V1V2_rf = calculate_rf_distance_random_pairs(VFull, V1V2, cores = CORES)
V1V3_rf = calculate_rf_distance_random_pairs(VFull, V1V3, cores = CORES)
V3V4_rf = calculate_rf_distance_random_pairs(VFull, V3V4, cores = CORES)
V3V5_rf = calculate_rf_distance_random_pairs(VFull, V3V5, cores = CORES)
V4_rf = calculate_rf_distance_random_pairs(VFull, V4, cores = CORES)
V4V5_rf = calculate_rf_distance_random_pairs(VFull, V4V5, cores = CORES)
nd_mix_rf = calculate_rf_distance_random_pairs(VFull, nd_mix, cores = CORES)

V1V2_guides_rf = calculate_rf_distance_random_pairs(VFull, V1V2_guides, cores = CORES)
V1V3_guides_rf = calculate_rf_distance_random_pairs(VFull, V1V3_guides, cores = CORES)
V3V4_guides_rf = calculate_rf_distance_random_pairs(VFull, V3V4_guides, cores = CORES)
V3V5_guides_rf = calculate_rf_distance_random_pairs(VFull, V3V5_guides, cores = CORES)
V4_guides_rf = calculate_rf_distance_random_pairs(VFull, V4_guides, cores = CORES)
V4V5_guides_rf = calculate_rf_distance_random_pairs(VFull, V4V5_guides, cores = CORES)
nd_mix_guides_rf = calculate_rf_distance_random_pairs(VFull, nd_mix_guides, cores = CORES)

nd_rf = data.frame("VFull" = VFull_control,
                   "VFull_guides" = VFull_guides,
                   
                   "V1V2_rf" = V1V2_rf,
                   "V1V3_rf" = V1V3_rf,
                   "V3V4_rf" = V3V4_rf,
                   "V3V5_rf" = V3V5_rf,
                   "V4_rf" = V4_rf,
                   "V4V5_rf" = V4V5_rf,
                   "nd_mix_rf" = nd_mix_rf,
                   
                   "V1V2_guides_rf" = V1V2_guides_rf,
                   "V1V3_guides_rf" = V1V3_guides_rf,
                   "V3V4_guides_rf" = V3V4_guides_rf,
                   "V3V5_guides_rf" = V3V5_guides_rf,
                   "V4_guides_rf" = V4_guides_rf,
                   "V4V5_guides_rf" = V4V5_guides_rf,
                   "nd_mix_guides_rf" = nd_mix_guides_rf)

setwd(dir.sim.out)
#saveRDS(file = "nd_rf_random_pairs.rds", object = nd_rf)
nd_rf = readRDS("nd_rf_random_pairs.rds")
nd_rf
as_tibble(readRDS("nd_rf_random_pairs.rds"))

################################################################################
# 1B: Determine how guides impact topological reconstruction nd mix -- JRF
################################################################################
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_nd_jrf.R' -r nd_jrf -p 30 -q all.q
as_tibble(readRDS("sim_analysis_nd_jrf.rds"))

################################################################################
# 1C: Determine how guidesimpact topological reconstruction nonoverlapping-- NYE
################################################################################

hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_nd_nye.R' -r nd_nye -p 30 -q biomed
as_tibble(readRDS("sim_analysis_nd_nye.rds"))

################################################################################
# 1D: Determine how guides x length impact topological reconstruction -- MSD
################################################################################
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_nd_msd.R' -r nd_msd -p 30 -q biomed
as_tibble(readRDS("sim_analysis_nd_msd.rds"))

################################################################################
# 1E: Determine how guides x length impact topological reconstruction -- PID
################################################################################
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_nd_pid.R' -r nd_pid -p 30 -q all.q
as_tibble(readRDS("sim_analysis_nd_pid.rds"))

################################################################################
# 1F: Determine how guides x length impact topological reconstruction -- MCD
################################################################################
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_nd_mcd.R' -r nd_mcd -p 30 -q all.q
as_tibble(readRDS("sim_analysis_nd_mcd.rds"))


################################################################################
# 1F: Determine how guides x length impact topological reconstruction -- PATH
################################################################################
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_nd_path.R' -r nd_path -p 30 -q all.q
as_tibble(readRDS("sim_analysis_nd_path.rds"))

################################################################################
# 1G: Determine how guides x length impact topological reconstruction -- QD
################################################################################
VFull_paths = file.path(dir.hvrs["VFull"], file_list = list.files(path = dir.hvrs["VFull"], pattern = ".tre"))[1:N]
VFull_guides_paths = file.path(dir.hvrs["VFull_guides"], file_list = list.files(path = dir.hvrs["VFull_guides"], pattern = "align.tre"))[1:N]
VFull_control_paths = file.path(dir.hvrs["VFull_control"], file_list = list.files(path = dir.hvrs["VFull_control"], pattern = "align.tre"))[1:N]

V1V2_paths = file.path(dir.hvrs["V1V2"], file_list = list.files(path = dir.hvrs["V1V2"], pattern = ".tre"))[1:N]
V1V3_paths = file.path(dir.hvrs["V1V3"], file_list = list.files(path = dir.hvrs["V1V3"], pattern = ".tre"))[1:N]
V3V4_paths = file.path(dir.hvrs["V3V4"], file_list = list.files(path = dir.hvrs["V3V4"], pattern = ".tre"))[1:N]
V3V5_paths = file.path(dir.hvrs["V3V5"], file_list = list.files(path = dir.hvrs["V3V5"], pattern = ".tre"))[1:N]
V4_paths = file.path(dir.hvrs["V4"], file_list = list.files(path = dir.hvrs["V4"], pattern = ".tre"))[1:N]
V4V5_paths = file.path(dir.hvrs["V4V5"], file_list = list.files(path = dir.hvrs["V4V5"], pattern = ".tre"))[1:N]
nd_mix_paths = file.path(dir.hvrs["nd_mix"], file_list = list.files(path = dir.hvrs["nd_mix"], pattern = ".tre"))[1:N]

V1V2_guides_paths = file.path(dir.hvrs["V1V2_guides"], file_list = list.files(path = dir.hvrs["V1V2_guides"], pattern = ".tre"))[1:N]
V1V3_guides_paths = file.path(dir.hvrs["V1V3_guides"], file_list = list.files(path = dir.hvrs["V1V3_guides"], pattern = ".tre"))[1:N]
V3V4_guides_paths = file.path(dir.hvrs["V3V4_guides"], file_list = list.files(path = dir.hvrs["V3V4_guides"], pattern = ".tre"))[1:N]
V3V5_guides_paths = file.path(dir.hvrs["V3V5_guides"], file_list = list.files(path = dir.hvrs["V3V5_guides"], pattern = ".tre"))[1:N]
V4_guides_paths = file.path(dir.hvrs["V4_guides"], file_list = list.files(path = dir.hvrs["V4_guides"], pattern = ".tre"))[1:N]
V4V5_guides_paths = file.path(dir.hvrs["V4V5_guides"], file_list = list.files(path = dir.hvrs["V4V5_guides"], pattern = ".tre"))[1:N]
nd_mix_guides_paths = file.path(dir.hvrs["nd_mix_guides"], file_list = list.files(path = dir.hvrs["nd_mix_guides"], pattern = ".tre"))[1:N]

VFull_guides_paths =  write_ref_dropped_trees(file_list = VFull_guides_paths, pattern = "^REF.")
V1V2_guides_paths =  write_ref_dropped_trees(file_list = V1V2_guides_paths, pattern = "^REF.")
V1V3_guides_paths =  write_ref_dropped_trees(file_list = V1V3_guides_paths, pattern = "^REF.")
V3V4_guides_paths =  write_ref_dropped_trees(file_list = V3V4_guides_paths, pattern = "^REF.")
V3V5_guides_paths =  write_ref_dropped_trees(file_list = V3V5_guides_paths, pattern = "^REF.")
V4_guides_paths =  write_ref_dropped_trees(file_list = V4_guides_paths, pattern = "^REF.")
V4V5_guides_paths =  write_ref_dropped_trees(file_list = V4V5_guides_paths, pattern = "^REF.")
nd_mix_guides_paths =  write_ref_dropped_trees(file_list = nd_mix_guides_paths, pattern = "^REF.")

CORES = 10
vfull_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = VFull_control_paths, cores = CORES)
vfull_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = VFull_guides_paths, cores = CORES)

V1V2_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V1V2_paths, cores = CORES)
V1V3_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V1V3_paths, cores = CORES)
V3V4_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V3V4_paths, cores = CORES)
V3V5_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V3V5_paths, cores = CORES)
V4_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V4_paths, cores = CORES)
V4V5_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V4V5_paths, cores = CORES)
nd_mix_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = nd_mix_paths, cores = CORES)

V1V2_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V1V2_guides_paths, cores = CORES)
V1V3_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V1V3_guides_paths, cores = CORES)
V3V4_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V3V4_guides_paths, cores = CORES)
V3V5_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V3V5_guides_paths, cores = CORES)
V4_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V4_guides_paths, cores = CORES)
V4V5_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = V4V5_guides_paths, cores = CORES)
nd_mix_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = VFull_paths, paths.trees2 = nd_mix_guides_paths, cores = CORES)


dist = data.frame("VFull" = vfull_quartet,
                  "VFull_guides" = vfull_guides_quartet,
                  
                  "V1V2_quartet" = V1V2_quartet,
                  "V1V3_quartet" = V1V3_quartet,
                  "V3V4_quartet" = V3V4_quartet,
                  "V3V5_quartet" = V3V5_quartet,
                  "V4_quartet" = V4_quartet,
                  "V4V5_quartet" = V4V5_quartet,
                  "nd_mix_quartet" = nd_mix_quartet,
                  
                  "V1V2_guides_quartet" = V1V2_guides_quartet,
                  "V1V3_guides_quartet" = V1V3_guides_quartet,
                  "V3V4_guides_quartet" = V3V4_guides_quartet,
                  "V3V5_guides_quartet" = V3V5_guides_quartet,
                  "V4_guides_quartet" = V4_guides_quartet,
                  "V4V5_guides_quartet" = V4V5_guides_quartet,
                  "nd_mix_guides_quartet" = nd_mix_guides_quartet)


setwd(dir.sim.out)
saveRDS(file = "nd_quartet_random_pairs.rds", object = dist, compress = TRUE)
QUARTET_MAX = choose(TREE_SIZE, 4)
saveRDS(file = "nd_quartet_normalized_random_pairs.rds", object = dist/QUARTET_MAX, compress = TRUE)

as_tibble(readRDS("nd_quartet_normalized_random_pairs.rds"))
as_tibble(readRDS("nd_quartet_random_pairs.rds"))


################################################################################
# 1H: Determine how guides x length impact topological reconstruction -- Nanostructure
################################################################################
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_nd_nanostructure_random_pairs_2.R' -r  nd_nanoclusters_random_pairs_2 -p 30 -q all.q
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_nd_nanostructure_random_pairs_3.R' -r  nd_nanoclusters_random_pairs_3 -p 30 -q all.q
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_nd_nanostructure_random_pairs_4.R' -r  nd_nanoclusters_random_pairs_4 -p 30 -q all.q
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_nd_nanostructure_random_pairs_5.R' -r  nd_nanoclusters_random_pairs_5 -p 30 -q all.q
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_nd_nanostructure_random_pairs_6.R' -r  nd_nanoclusters_random_pairs_6 -p 30 -q all.q
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_nd_nanostructure_random_pairs_7.R' -r  nd_nanoclusters_random_pairs_7 -p 30 -q all.q
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_nd_nanostructure_random_pairs_8.R' -r  nd_nanoclusters_random_pairs_8 -p 30 -q all.q
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_nd_nanostructure_random_pairs_9.R' -r  nd_nanoclusters_random_pairs_9 -p 30 -q all.q
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_nd_nanostructure_random_pairs_10.R' -r  nd_nanoclusters_random_pairs_10 -p 30 -q all.q
as_tibble(readRDS("nanostructure_nd_2_random_pairs.rds"))
as_tibble(readRDS("nanostructure_nd_3_random_pairs.rds"))
as_tibble(readRDS("nanostructure_nd_4_random_pairs.rds"))
as_tibble(readRDS("nanostructure_nd_5_random_pairs.rds"))
as_tibble(readRDS("nanostructure_nd_6_random_pairs.rds"))
as_tibble(readRDS("nanostructure_nd_7_random_pairs.rds"))
as_tibble(readRDS("nanostructure_nd_8_random_pairs.rds"))
as_tibble(readRDS("nanostructure_nd_9_random_pairs.rds"))
as_tibble(readRDS("nanostructure_nd_10_random_pairs.rds"))


#######################################
# 1J.Determine MAST
#######################################

hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/sim_analysis_nd_mast.R' -r  nd_mast -p 100 -q sharpton
as_tibble(readRDS("nd_mast_random_pairs.rds"))


#######################################
# Clean up
#######################################
rm(list = ls())

################################################################################
################################################################################
# EXPERIMENTAL DATA
# Aug 8th 2025
################################################################################
################################################################################
rm(list = ls())

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
library(ggbreak)
library(grid)

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
dir.sim.out = file.path(HOME, "2025_hvr_guide_phylogenetic_integration/out/experimental_demo/")

setwd(HOME_TMP)
dir.demo = file.path(HOME_TMP, "experimental_demo/")
dir.demo
setwd(dir.scripts)
source("sim_analysis_functions_V2.R")

################################################################################
# SET CONSTANTS
################################################################################
TREE_SIZE = 598
RF_MAX = 2*(TREE_SIZE -3)
QUARTET_MAX = choose(TREE_SIZE, 4)

# Check that all the tree files are present. 
directories <- c(list.dirs(dir.demo, full.names = TRUE, recursive = FALSE))
names(directories) = basename(directories)
basename(directories)

################################################################################
# EXP DEMO PREPREP: READ INPUT FILES
################################################################################
setwd(dir.demo)
dir.demos = directories[c("short", "long", "long_control", "short_guides_optimized")] 
dir.demos
long = read_trees(file.path(dir.demos["long"], file_list = list.files(path = dir.demos["long"], pattern = ".tre")))
long

long_control = read_trees(file.path(dir.demos["long_control"], file_list = list.files(path = dir.demos["long_control"], pattern = ".tre")))
long_control 

short = read_trees(file.path(dir.demos["short"], file_list = list.files(path = dir.demos["short"], pattern = ".tre")))
short

short_guides_opt = read_trees(file.path(dir.demos["short_guides_optimized"], file_list = list.files(path = dir.demos["short_guides_optimized"], pattern = "align.tre")), drop = TRUE)
short_guides_opt


#########################
# EXP DEMO RF Distance - Random pairs
#########################
CORES = 10 
short_rf = calculate_rf_distance_random_pairs(long, short, cores = CORES)
long_control = calculate_rf_distance_random_pairs(long, long_control, cores = CORES)
short_guides_opt_rf = calculate_rf_distance_random_pairs(long, short_guides_opt, cores = CORES)
exp_random_pairs = 
  data.frame("short" = short_rf,
             "control" = long_control,
             "short_guides_opt" = short_guides_opt_rf)
# Convert to long format for ggplot
exp_random_pairs_long <- exp_random_pairs %>%
  #select(!control) %>%
  pivot_longer(cols = everything(),
               names_to = "Condition",
               values_to = "Value")

# Create box plot
p = ggplot(exp_random_pairs_long, aes(x = Condition, y = Value, fill = Condition)) +
  geom_boxplot() +
  labs(title = "Random Pair RF Distributions",
       x = "Condition",
       y = "Value") +
  theme(legend.position = "none")+ 
  scale_y_break(c(0.01, 0.51))
p

exp_random_pairs
setwd(dir.out)
saveRDS(exp_random_pairs, file = "experimental_demo_rf_random_pairs.rds", compress = TRUE)
setwd("../../data_visulization/")
as_tibble(readRDS("experimental_demo_rf_random_pairs.rds"))
#########################
# EXP DEMO Quartet Distance - For All possible pairwise comparisons
#########################
short_paths = file.path(dir.demos["short"], file_list = list.files(path = dir.demos["short"], pattern = ".tre"))
long_paths = file.path(dir.demos["long"], file_list = list.files(path = dir.demos["long"], pattern = ".tre"))
long_control_paths = file.path(dir.demos["long_control"], file_list = list.files(path = dir.demos["long_control"], pattern = ".tre"))
short_guides_paths = file.path(dir.demos["short_guides_optimized"], file_list = list.files(path = dir.demos["short_guides_optimized"], pattern = "align.tre"))
short_guides_paths = write_ref_dropped_trees(file_list = short_guides_paths, pattern = "^REF.")

short_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = long_paths, paths.trees2 = short_paths, cores = CORES)
short_guides_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = long_paths, paths.trees2 = short_guides_paths, cores = CORES)
long_control_quartet = calculate_quartet_distances_random_pairs(paths.trees1 = long_paths, paths.trees2 = long_control_paths, cores = CORES)
CORES = 10

exp_demo_quartet =
  data.frame(
    "short" = short_quartet,
    "short_guides" = short_guides_quartet,
    "control" = long_control_quartet)
head(exp_demo_quartet)

setwd(dir.out)
saveRDS(object = exp_demo_quartet, "experimental_demo_quartet.rds", compress = TRUE)
saveRDS(object = exp_demo_quartet/QUARTET_MAX, "experimental_demo_quartet_normalized.rds", compress = TRUE)
as_tibble(readRDS("experimental_demo_quartet.rds"))
as_tibble(readRDS("experimental_demo_quartet_normalized.rds"))

#########################
# EXP DEMO Path Distance -- Between 100 random pairings.
#########################
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/experimental_demo_path_random_pairings.R' -r experimental_demo_path_random_pairings -p 10 -q sharpton
setwd(dir.out)
path = readRDS("experimental_demo_path_random_pairings.rds")
head(path)

#########################
# EXPERIMENTAL DEMO Phylogenetic Information Distance (PID) - Random pairings.
#########################
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/experimental_demo_pid_random_pairs.R' -r experimental_demo_random_pairs -p 10 -q sharpton
setwd(dir.out)
pid = readRDS("experimental_demo_pid_random_pairs.rds")
head(pid)


#########################
# EXPERIMENTAL DEMO MSD -- Random pairings.
#########################
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/experimental_demo_msd_random_pairs.R' -r msd_random_pairs -p 10 -q sharpton
setwd(dir.out)
experimental_demo_msd_random_pairings = readRDS("experimental_demo_msd_random_pairs.rds")
head(experimental_demo_msd_random_pairings)

#########################
# EXPERIMENTAL DEMO NYE Random Pairings
#########################
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/experimental_demo_nye_random_paris.R' -r experimental_demo_nye_random_pairs -p 10 -q sharpton
setwd(dir.out)
experimental_demo_nye_random_pairings = readRDS("experimental_demo_nye_random_pairs.rds")
head(experimental_demo_nye_random_pairings)

#########################
# 1G. JRF - Random pairings
#########################
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/experimental_demo_jrf_random_pairs.R' -r experimental_demo_jrf_random_pairs -p 10 -q biomed
setwd(dir.out)
experimental_demo_jrf_random_pairings = readRDS("experimental_demo_jrf_random_pairs.rds")
experimental_demo_jrf_random_pairings

#######################################
# 1H. MCD --Random pairings
#######################################
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/experimental_demo_mcd_random_pairs.R' -r mcd_random_pairs -p 10 -q sharpton
setwd(dir.out)
experimental_demo_mcd_random_pairings = readRDS("experimental_demo_mcd_random_pairs.rds")
head(experimental_demo_mcd_random_pairings)


#######################################
# 1J. Determine subcluster distance - Random pairs S2 - S10
#######################################
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/experimental_demo_nanostructure_random_pairs.R' -r  experimental_demo_nanostructure_random_pairs.R -p 10 -q sharpton
setwd(dir.out)
experimental_demo_nonostructure_random_pairings = readRDS("experimental_demo_subtrees_random_pairs.rds")
head(experimental_demo_nonostructure_random_pairings)
dim(experimental_demo_nonostructure_random_pairings)

#######################################
# 1J.Determine MAST
#######################################
hpcman queue submit 'Rscript /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/experimental_demo_mast_random_pairs.R' -r  mast -p 10 -q sharpton
setwd(dir.out)
experimental_demo_mast_random_pairings = readRDS("experimental_demo_mast_random_pairs.rds")
head(experimental_demo_mast_random_pairings)
hist(experimental_demo_mast_random_pairings$short_guides)

#######################################
# Clean up
#######################################
rm(list = ls())



