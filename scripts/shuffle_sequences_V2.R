
#AUTHOR: ARNOLD
#DAY: Jan 13th, 2025
# shuffle_sequences_V2.R
# The purpose of this script is to take each of the input files and shuffle 
# their sequences

# LIBRARIES
library(phylotools)
library(dplyr)
.libPaths(c(.libPaths(), "/home/micro/arnoldho/R/x86_64-pc-linux-gnu-library/4.1"))
library("PIMMA")
library(foreach)
library(doParallel)

# DIRECTORIES
HOME = "/nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration"
HOME_TMP = "/nfs3/Sharpton_Lab/tmp/projects/arnoldh/2025_hvr_guide_phylogenetic_integration/"

setwd(HOME)
dir.scripts = file.path(HOME, "2025_hvr_guide_phylogenetic_integration/scripts/")
dir.ref = file.path(HOME, "reference_data/")

setwd(HOME_TMP)
dir.simulation = file.path(HOME_TMP, "length_simulations/")
dir.mix = file.path(HOME_TMP, "mix_experiment/")
dir.demo = file.path(HOME_TMP, "experimental_demo/")
TEMPLATE_SEED = file.path(dir.ref, "ref.align")

if(!file.exists(TEMPLATE_SEED)){
  print("FAIL: Cannot find template seed file.")
}else{print("PASS: Can find the seed file.")}

# Load scripts
setwd(dir.scripts)
source("shuffle_sequences_functions_V2.R")

# Shuffle sequences for the full length alignment
setwd(dir.simulation)

# Full length
# write_shuffles(in_fasta = file.path(dir.ref, "ref.align"), N = 100, out_fasta_prefix = "VFull_", out_dir = file.path(dir.simulation, "V_Full_Length"))
# write_shuffles(in_fasta = file.path(dir.ref, "ref.align"), N = 100, out_fasta_prefix = "VFull_control_", out_dir = file.path(dir.simulation, "VFull_control"))

# # # 
# # # V1 sims, no guide
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V1_L50.align"), N = 100, out_fasta_prefix = "V1_L50_", out_dir = file.path(dir.simulation, "V1_L50"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V1_L100.align"), N = 100, out_fasta_prefix = "V1_L100_", out_dir = file.path(dir.simulation, "V1_L100"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V1_L200.align"), N = 100, out_fasta_prefix = "V1_L200_", out_dir = file.path(dir.simulation, "V1_L200"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V1_L300.align"), N = 100, out_fasta_prefix = "V1_L300_", out_dir = file.path(dir.simulation, "V1_L300"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V1_L400.align"), N = 100, out_fasta_prefix = "V1_L400_", out_dir = file.path(dir.simulation, "V1_L400"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V1_L500.align"), N = 100, out_fasta_prefix = "V1_L500_", out_dir = file.path(dir.simulation, "V1_L500"))
# #
# # # V2 sims, no guide
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V2_L50.align"), N = 100, out_fasta_prefix = "V2_L50_", out_dir = file.path(dir.simulation, "V2_L50"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V2_L100.align"), N = 100, out_fasta_prefix = "V2_L100_", out_dir = file.path(dir.simulation, "V2_L100"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V2_L200.align"), N = 100, out_fasta_prefix = "V2_L200_", out_dir = file.path(dir.simulation, "V2_L200"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V2_L300.align"), N = 100, out_fasta_prefix = "V2_L300_", out_dir = file.path(dir.simulation, "V2_L300"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V2_L400.align"), N = 100, out_fasta_prefix = "V2_L400_", out_dir = file.path(dir.simulation, "V2_L400"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V2_L500.align"), N = 100, out_fasta_prefix = "V2_L500_", out_dir = file.path(dir.simulation, "V2_L500"))
# #
# # # V3 sims, no guide
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V3_L50.align"), N = 100, out_fasta_prefix = "V3_L50_", out_dir = file.path(dir.simulation, "V3_L50"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V3_L100.align"), N = 100, out_fasta_prefix = "V3_L100_", out_dir = file.path(dir.simulation, "V3_L100"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V3_L200.align"), N = 100, out_fasta_prefix = "V3_L200_", out_dir = file.path(dir.simulation, "V3_L200"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V3_L300.align"), N = 100, out_fasta_prefix = "V3_L300_", out_dir = file.path(dir.simulation, "V3_L300"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V3_L400.align"), N = 100, out_fasta_prefix = "V3_L400_", out_dir = file.path(dir.simulation, "V3_L400"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V3_L500.align"), N = 100, out_fasta_prefix = "V3_L500_", out_dir = file.path(dir.simulation, "V3_L500"))
# 
# # # V4 sims, no guide
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V4_L50.align"), N = 100, out_fasta_prefix = "V4_L50_", out_dir = file.path(dir.simulation, "V4_L50"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V4_L100.align"), N = 100, out_fasta_prefix = "V4_L100_", out_dir = file.path(dir.simulation, "V4_L100"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V4_L200.align"), N = 100, out_fasta_prefix = "V4_L200_", out_dir = file.path(dir.simulation, "V4_L200"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V4_L300.align"), N = 100, out_fasta_prefix = "V4_L300_", out_dir = file.path(dir.simulation, "V4_L300"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V4_L400.align"), N = 100, out_fasta_prefix = "V4_L400_", out_dir = file.path(dir.simulation, "V4_L400"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V4_L500.align"), N = 100, out_fasta_prefix = "V4_L500_", out_dir = file.path(dir.simulation, "V4_L500"))
# #
# # # V5 sims, no guide
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V5_L50.align"), N = 100, out_fasta_prefix = "V5_L50_", out_dir = file.path(dir.simulation, "V5_L50"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V5_L100.align"), N = 100, out_fasta_prefix = "V5_L100_", out_dir = file.path(dir.simulation, "V5_L100"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V5_L200.align"), N = 100, out_fasta_prefix = "V5_L200_", out_dir = file.path(dir.simulation, "V5_L200"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V5_L300.align"), N = 100, out_fasta_prefix = "V5_L300_", out_dir = file.path(dir.simulation, "V5_L300"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V5_L400.align"), N = 100, out_fasta_prefix = "V5_L400_", out_dir = file.path(dir.simulation, "V5_L400"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V5_L500.align"), N = 100, out_fasta_prefix = "V5_L500_", out_dir = file.path(dir.simulation, "V5_L500"))
# #
# # # V6 sims, no guide
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V6_L50.align"), N = 100, out_fasta_prefix = "V6_L50_", out_dir = file.path(dir.simulation, "V6_L50"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V6_L100.align"), N = 100, out_fasta_prefix = "V6_L100_", out_dir = file.path(dir.simulation, "V6_L100"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V6_L200.align"), N = 100, out_fasta_prefix = "V6_L200_", out_dir = file.path(dir.simulation, "V6_L200"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V6_L300.align"), N = 100, out_fasta_prefix = "V6_L300_", out_dir = file.path(dir.simulation, "V6_L300"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V6_L400.align"), N = 100, out_fasta_prefix = "V6_L400_", out_dir = file.path(dir.simulation, "V6_L400"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V6_L500.align"), N = 100, out_fasta_prefix = "V6_L500_", out_dir = file.path(dir.simulation, "V6_L500"))
# #
# # # V7 sims, no guide
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V7_L50.align"), N = 100, out_fasta_prefix = "V7_L50_", out_dir = file.path(dir.simulation, "V7_L50"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V7_L100.align"), N = 100, out_fasta_prefix = "V7_L100_", out_dir = file.path(dir.simulation, "V7_L100"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V7_L200.align"), N = 100, out_fasta_prefix = "V7_L200_", out_dir = file.path(dir.simulation, "V7_L200"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V7_L300.align"), N = 100, out_fasta_prefix = "V7_L300_", out_dir = file.path(dir.simulation, "V7_L300"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V7_L400.align"), N = 100, out_fasta_prefix = "V7_L400_", out_dir = file.path(dir.simulation, "V7_L400"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V7_L500.align"), N = 100, out_fasta_prefix = "V7_L500_", out_dir = file.path(dir.simulation, "V7_L500"))
# #
# # # V8 sims, no guide
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V8_L50.align"), N = 100, out_fasta_prefix = "V8_L50_", out_dir = file.path(dir.simulation, "V8_L50"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V8_L100.align"), N = 100, out_fasta_prefix = "V8_L100_", out_dir = file.path(dir.simulation, "V8_L100"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V8_L200.align"), N = 100, out_fasta_prefix = "V8_L200_", out_dir = file.path(dir.simulation, "V8_L200"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V8_L300.align"), N = 100, out_fasta_prefix = "V8_L300_", out_dir = file.path(dir.simulation, "V8_L300"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V8_L400.align"), N = 100, out_fasta_prefix = "V8_L400_", out_dir = file.path(dir.simulation, "V8_L400"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V8_L500.align"), N = 100, out_fasta_prefix = "V8_L500_", out_dir = file.path(dir.simulation, "V8_L500"))
# #
# # # V9 sims, no guide
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V9_L50.align"), N = 100, out_fasta_prefix = "V9_L50_", out_dir = file.path(dir.simulation, "V9_L50"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V9_L100.align"), N = 100, out_fasta_prefix = "V9_L100_", out_dir = file.path(dir.simulation, "V9_L100"))

# V1 sims, with guide
#write_shuffles(in_fasta = file.path(dir.ref, "seqs_V1_L50_guides.align"), N = 100, out_fasta_prefix = "V1_L50_guides_", out_dir = file.path(dir.simulation, "V1_L50_guides"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V1_L100_guides.align"), N = 100, out_fasta_prefix = "V1_L100_guides_", out_dir = file.path(dir.simulation, "V1_L100_guides"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V1_L200_guides.align"), N = 100, out_fasta_prefix = "V1_L200_guides_", out_dir = file.path(dir.simulation, "V1_L200_guides"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V1_L300_guides.align"), N = 100, out_fasta_prefix = "V1_L300_guides_", out_dir = file.path(dir.simulation, "V1_L300_guides"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V1_L400_guides.align"), N = 100, out_fasta_prefix = "V1_L400_guides_", out_dir = file.path(dir.simulation, "V1_L400_guides"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V1_L500_guides.align"), N = 100, out_fasta_prefix = "V1_L500_guides_", out_dir = file.path(dir.simulation, "V1_L500_guides"))
# 
# # V2 sims, with guide
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V2_L50_guides.align"), N = 100, out_fasta_prefix = "V2_L50_guides_", out_dir = file.path(dir.simulation, "V2_L50_guides"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V2_L100_guides.align"), N = 100, out_fasta_prefix = "V2_L100_guides_", out_dir = file.path(dir.simulation, "V2_L100_guides"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V2_L200_guides.align"), N = 100, out_fasta_prefix = "V2_L200_guides_", out_dir = file.path(dir.simulation, "V2_L200_guides"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V2_L300_guides.align"), N = 100, out_fasta_prefix = "V2_L300_guides_", out_dir = file.path(dir.simulation, "V2_L300_guides"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V2_L400_guides.align"), N = 100, out_fasta_prefix = "V2_L400_guides_", out_dir = file.path(dir.simulation, "V2_L400_guides"))
# write_shuffles(in_fasta = file.path(dir.ref, "seqs_V2_L500_guides.align"), N = 100, out_fasta_prefix = "V2_L500_guides_", out_dir = file.path(dir.simulation, "V2_L500_guides"))

# V3 sims, with guide
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V3_L50_guides.align"), N = 100, out_fasta_prefix = "V3_L50_guides_", out_dir = file.path(dir.simulation, "V3_L50_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V3_L100_guides.align"), N = 100, out_fasta_prefix = "V3_L100_guides_", out_dir = file.path(dir.simulation, "V3_L100_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V3_L200_guides.align"), N = 100, out_fasta_prefix = "V3_L200_guides_", out_dir = file.path(dir.simulation, "V3_L200_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V3_L300_guides.align"), N = 100, out_fasta_prefix = "V3_L300_guides_", out_dir = file.path(dir.simulation, "V3_L300_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V3_L400_guides.align"), N = 100, out_fasta_prefix = "V3_L400_guides_", out_dir = file.path(dir.simulation, "V3_L400_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V3_L500_guides.align"), N = 100, out_fasta_prefix = "V3_L500_guides_", out_dir = file.path(dir.simulation, "V3_L500_guides"))

# V4 sims, with guide
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V4_L50_guides.align"), N = 100, out_fasta_prefix = "V4_L50_guides_", out_dir = file.path(dir.simulation, "V4_L50_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V4_L100_guides.align"), N = 100, out_fasta_prefix = "V4_L100_guides_", out_dir = file.path(dir.simulation, "V4_L100_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V4_L200_guides.align"), N = 100, out_fasta_prefix = "V4_L200_guides_", out_dir = file.path(dir.simulation, "V4_L200_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V4_L300_guides.align"), N = 100, out_fasta_prefix = "V4_L300_guides_", out_dir = file.path(dir.simulation, "V4_L300_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V4_L400_guides.align"), N = 100, out_fasta_prefix = "V4_L400_guides_", out_dir = file.path(dir.simulation, "V4_L400_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V4_L500_guides.align"), N = 100, out_fasta_prefix = "V4_L500_guides_", out_dir = file.path(dir.simulation, "V4_L500_guides"))

# V5 sims, with guide
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V5_L50_guides.align"), N = 100, out_fasta_prefix = "V5_L50_guides_", out_dir = file.path(dir.simulation, "V5_L50_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V5_L100_guides.align"), N = 100, out_fasta_prefix = "V5_L100_guides_", out_dir = file.path(dir.simulation, "V5_L100_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V5_L200_guides.align"), N = 100, out_fasta_prefix = "V5_L200_guides_", out_dir = file.path(dir.simulation, "V5_L200_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V5_L300_guides.align"), N = 100, out_fasta_prefix = "V5_L300_guides_", out_dir = file.path(dir.simulation, "V5_L300_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V5_L400_guides.align"), N = 100, out_fasta_prefix = "V5_L400_guides_", out_dir = file.path(dir.simulation, "V5_L400_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V5_L500_guides.align"), N = 100, out_fasta_prefix = "V5_L500_guides_", out_dir = file.path(dir.simulation, "V5_L500_guides"))

# V6 sims, with guide
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V6_L50_guides.align"), N = 100, out_fasta_prefix = "V6_L50_guides_", out_dir = file.path(dir.simulation, "V6_L50_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V6_L100_guides.align"), N = 100, out_fasta_prefix = "V6_L100_guides_", out_dir = file.path(dir.simulation, "V6_L100_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V6_L200_guides.align"), N = 100, out_fasta_prefix = "V6_L200_guides_", out_dir = file.path(dir.simulation, "V6_L200_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V6_L300_guides.align"), N = 100, out_fasta_prefix = "V6_L300_guides_", out_dir = file.path(dir.simulation, "V6_L300_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V6_L400_guides.align"), N = 100, out_fasta_prefix = "V6_L400_guides_", out_dir = file.path(dir.simulation, "V6_L400_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V6_L500_guides.align"), N = 100, out_fasta_prefix = "V6_L500_guides_", out_dir = file.path(dir.simulation, "V6_L500_guides"))

# V7 sims, with guide
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V7_L50_guides.align"), N = 100, out_fasta_prefix = "V7_L50_guides_", out_dir = file.path(dir.simulation, "V7_L50_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V7_L100_guides.align"), N = 100, out_fasta_prefix = "V7_L100_guides_", out_dir = file.path(dir.simulation, "V7_L100_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V7_L200_guides.align"), N = 100, out_fasta_prefix = "V7_L200_guides_", out_dir = file.path(dir.simulation, "V7_L200_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V7_L300_guides.align"), N = 100, out_fasta_prefix = "V7_L300_guides_", out_dir = file.path(dir.simulation, "V7_L300_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V7_L400_guides.align"), N = 100, out_fasta_prefix = "V7_L400_guides_", out_dir = file.path(dir.simulation, "V7_L400_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V7_L500_guides.align"), N = 100, out_fasta_prefix = "V7_L500_guides_", out_dir = file.path(dir.simulation, "V7_L500_guides"))

# V8 sims, with guide
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V8_L50_guides.align"), N = 100, out_fasta_prefix = "V8_L50_guides_", out_dir = file.path(dir.simulation, "V8_L50_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V8_L100_guides.align"), N = 100, out_fasta_prefix = "V8_L100_guides_", out_dir = file.path(dir.simulation, "V8_L100_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V8_L200_guides.align"), N = 100, out_fasta_prefix = "V8_L200_guides_", out_dir = file.path(dir.simulation, "V8_L200_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V8_L300_guides.align"), N = 100, out_fasta_prefix = "V8_L300_guides_", out_dir = file.path(dir.simulation, "V8_L300_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V8_L400_guides.align"), N = 100, out_fasta_prefix = "V8_L400_guides_", out_dir = file.path(dir.simulation, "V8_L400_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V8_L500_guides.align"), N = 100, out_fasta_prefix = "V8_L500_guides_", out_dir = file.path(dir.simulation, "V8_L500_guides"))

# V9 sims, with guide
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V9_L50_guides.align"), N = 100, out_fasta_prefix = "V9_L50_guides_", out_dir = file.path(dir.simulation, "V9_L50_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V9_L100_guides.align"), N = 100, out_fasta_prefix = "V9_L100_guides_", out_dir = file.path(dir.simulation, "V9_L100_guides"))

# HVR Sims without guide
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V1_LHVR.align"), N = 100, out_fasta_prefix = "V1_HVR_", out_dir = file.path(dir.simulation, "V1_HVR"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V2_LHVR.align"), N = 100, out_fasta_prefix = "V2_HVR_", out_dir = file.path(dir.simulation, "V2_HVR"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V3_LHVR.align"), N = 100, out_fasta_prefix = "V3_HVR_", out_dir = file.path(dir.simulation, "V3_HVR"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V4_LHVR.align"), N = 100, out_fasta_prefix = "V4_HVR_", out_dir = file.path(dir.simulation, "V4_HVR"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V5_LHVR.align"), N = 100, out_fasta_prefix = "V5_HVR_", out_dir = file.path(dir.simulation, "V5_HVR"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V6_LHVR.align"), N = 100, out_fasta_prefix = "V6_HVR_", out_dir = file.path(dir.simulation, "V6_HVR"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V7_LHVR.align"), N = 100, out_fasta_prefix = "V7_HVR_", out_dir = file.path(dir.simulation, "V7_HVR"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V8_LHVR.align"), N = 100, out_fasta_prefix = "V8_HVR_", out_dir = file.path(dir.simulation, "V8_HVR"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V9_LHVR.align"), N = 100, out_fasta_prefix = "V9_HVR_", out_dir = file.path(dir.simulation, "V9_HVR"))

# HVR Sims with guide
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V1_LHVR_guides.align"), N = 100, out_fasta_prefix = "V1_HVR_guides_", out_dir = file.path(dir.simulation, "V1_HVR_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V2_LHVR_guides.align"), N = 100, out_fasta_prefix = "V2_HVR_guides_", out_dir = file.path(dir.simulation, "V2_HVR_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V3_LHVR_guides.align"), N = 100, out_fasta_prefix = "V3_HVR_guides_", out_dir = file.path(dir.simulation, "V3_HVR_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V4_LHVR_guides.align"), N = 100, out_fasta_prefix = "V4_HVR_guides_", out_dir = file.path(dir.simulation, "V4_HVR_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V5_LHVR_guides.align"), N = 100, out_fasta_prefix = "V5_HVR_guides_", out_dir = file.path(dir.simulation, "V5_HVR_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V6_LHVR_guides.align"), N = 100, out_fasta_prefix = "V6_HVR_guides_", out_dir = file.path(dir.simulation, "V6_HVR_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V7_LHVR_guides.align"), N = 100, out_fasta_prefix = "V7_HVR_guides_", out_dir = file.path(dir.simulation, "V7_HVR_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V8_LHVR_guides.align"), N = 100, out_fasta_prefix = "V8_HVR_guides_", out_dir = file.path(dir.simulation, "V8_HVR_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V9_LHVR_guides.align"), N = 100, out_fasta_prefix = "V9_HVR_guides_", out_dir = file.path(dir.simulation, "V9_HVR_guides"))

# Mixed VR regions with and without guides for single VR regions
setwd(dir.mix)
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V1V2_LHVR.align"), N = 100, out_fasta_prefix = "V1V2_", out_dir = file.path(dir.mix, "V1V2"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V1V2_LHVR_guides.align"), N = 100, out_fasta_prefix = "V1V2_guides_", out_dir = file.path(dir.mix, "V1V2_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V1V3_LHVR.align"), N = 100, out_fasta_prefix = "V1V3_", out_dir = file.path(dir.mix, "V1V3"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V1V3_LHVR_guides.align"), N = 100, out_fasta_prefix = "V1V3_guides_", out_dir = file.path(dir.mix, "V1V3_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V3V4_LHVR.align"), N = 100, out_fasta_prefix = "V3V4_", out_dir = file.path(dir.mix, "V3V4"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V3V4_LHVR_guides.align"), N = 100, out_fasta_prefix = "V3V4_guides_", out_dir = file.path(dir.mix, "V3V4_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V4_LHVR.align"), N = 100, out_fasta_prefix = "V4_", out_dir = file.path(dir.mix, "V4"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V4_LHVR_guides.align"), N = 100, out_fasta_prefix = "V4_guides_", out_dir = file.path(dir.mix, "V4_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V4V5_LHVR.align"), N = 100, out_fasta_prefix = "V4V5_", out_dir = file.path(dir.mix, "V4V5"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V4V5_LHVR_guides.align"), N = 100, out_fasta_prefix = "V4V5_guides_", out_dir = file.path(dir.mix, "V4V5_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V5V7_LHVR.align"), N = 100, out_fasta_prefix = "V5V7_", out_dir = file.path(dir.mix, "V5V7"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V5V7_LHVR_guides.align"), N = 100, out_fasta_prefix = "V5V7_guides_", out_dir = file.path(dir.mix, "V5V7_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V7V9_LHVR.align"), N = 100, out_fasta_prefix = "V7V9_", out_dir = file.path(dir.mix, "V7V9"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V7V9_LHVR_guides.align"), N = 100, out_fasta_prefix = "V7V9_guides_", out_dir = file.path(dir.mix, "V7V9_guides"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V3V5_LHVR.align"), N = 100, out_fasta_prefix = "V3V5_", out_dir = file.path(dir.mix, "V3V5"))
write_shuffles(in_fasta = file.path(dir.ref, "seqs_V3V5_LHVR_guides.align"), N = 100, out_fasta_prefix = "V3V5_guides_", out_dir = file.path(dir.mix, "V3V5_guides"))

# Alignments for mixed regions
# # The V4 and V3 V4 Mix Overlapping regions
fastas = c(file.path(dir.ref, "seqs_V3V4_LHVR.align"),
           file.path(dir.ref, "seqs_V4_LHVR.align"))
names(fastas) = c("V3V4", "V4")
fastas
write_shuffles_mix(in_fasta = fastas, N = 100, out_fasta_prefix = "V4_V3V4_mix_", out_dir = file.path(dir.mix, "V4_V3V4_mix/"), guides = FALSE)


# The V1V2, V3V4, and V5V7 Non-overlapping experiement
fastas = c(file.path(dir.ref, "seqs_V1V2_LHVR.align"),
           file.path(dir.ref, "seqs_V3V4_LHVR.align"),
           file.path(dir.ref, "seqs_V5V7_LHVR.align"))
names(fastas) = c("V1V2", "V3V4", "V5V7")
fastas
write_shuffles_mix(in_fasta = fastas, N = 100, out_fasta_prefix = "V1V2_V3V4_V5V7_mix_", out_dir = file.path(dir.mix, "V1V2_V3V4_V5V7_mix/"), guides = FALSE)

# 
# Show how we can combined a mix set of sequences 
fastas = c(file.path(dir.ref, "seqs_V1V2_LHVR.align"),
           file.path(dir.ref, "seqs_V1V3_LHVR.align"),
           file.path(dir.ref, "seqs_V3V4_LHVR.align"),
           file.path(dir.ref, "seqs_V3V5_LHVR.align"),
           file.path(dir.ref, "seqs_V4_LHVR.align"),
           file.path(dir.ref, "seqs_V4V5_LHVR.align"))
names(fastas) = c("V1V2", "V1V3", "V3V4", "V3V5", "V4", "V4V5")
fastas
write_shuffles_mix(in_fasta = fastas, N = 100, out_fasta_prefix = "nd_mix_", out_dir = file.path(dir.mix, "nd_mix/"), guides = FALSE)


fastas = c(file.path(dir.ref, "seqs_V1_LHVR.align"),
           file.path(dir.ref, "seqs_V2_LHVR.align"),
           file.path(dir.ref, "seqs_V3_LHVR.align"),
           file.path(dir.ref, "seqs_V4_LHVR.align"),
           file.path(dir.ref, "seqs_V5_LHVR.align"),
           file.path(dir.ref, "seqs_V6_LHVR.align"),
           file.path(dir.ref, "seqs_V7_LHVR.align"),
           file.path(dir.ref, "seqs_V8_LHVR.align"),
           file.path(dir.ref, "seqs_V9_LHVR.align"))
names(fastas) = c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9")
fastas
write_shuffles_mix(in_fasta = fastas, N = 100, out_fasta_prefix = "V1V9_mix_", out_dir = file.path(dir.mix, "V1V9_mix/"), guides = FALSE)



# The V4 and V3 V4 Mix Overlapping regions with guides
fastas = c(file.path(dir.ref, "seqs_V3V4_LHVR_guides.align"),
           file.path(dir.ref, "seqs_V4_LHVR_guides.align"))
names(fastas) = c("V3V4", "V4")
fastas
write_shuffles_mix(in_fasta = fastas, N = 100, out_fasta_prefix = "V4_V3V4_mix_guides_", out_dir = file.path(dir.mix, "V4_V3V4_mix_guides/"), guides = TRUE)


# The V1V2, V3V4, and V5V7 Non-overlapping experiment
fastas = c(file.path(dir.ref, "seqs_V1V2_LHVR_guides.align"),
           file.path(dir.ref, "seqs_V3V4_LHVR_guides.align"),
           file.path(dir.ref, "seqs_V5V7_LHVR_guides.align"))
names(fastas) = c("V1V2", "V3V4", "V5V7")
fastas
write_shuffles_mix(in_fasta = fastas, N = 100, out_fasta_prefix = "V1V2_V3V4_V5V7_mix_guides_", out_dir = file.path(dir.mix, "V1V2_V3V4_V5V7_mix_guides/"), guides = TRUE)


fastas = c(file.path(dir.ref, "seqs_V1V2_LHVR_guides.align"),
           file.path(dir.ref, "seqs_V1V3_LHVR_guides.align"),
           file.path(dir.ref, "seqs_V3V4_LHVR_guides.align"),
           file.path(dir.ref, "seqs_V3V5_LHVR_guides.align"),
           file.path(dir.ref, "seqs_V4_LHVR_guides.align"),
           file.path(dir.ref, "seqs_V4V5_LHVR_guides.align"))
names(fastas) = c("V1V2", "V1V3", "V3V4", "V3V5", "V4", "V4V5")
fastas
write_shuffles_mix(in_fasta = fastas, N = 100, out_fasta_prefix = "nd_mix_guides_", out_dir = file.path(dir.mix, "nd_mix_guides/"), guides = TRUE)


fastas = c(file.path(dir.ref, "seqs_V1_LHVR_guides.align"),
           file.path(dir.ref, "seqs_V2_LHVR_guides.align"),
           file.path(dir.ref, "seqs_V3_LHVR_guides.align"),
           file.path(dir.ref, "seqs_V4_LHVR_guides.align"),
           file.path(dir.ref, "seqs_V5_LHVR_guides.align"),
           file.path(dir.ref, "seqs_V6_LHVR_guides.align"),
           file.path(dir.ref, "seqs_V7_LHVR_guides.align"),
           file.path(dir.ref, "seqs_V8_LHVR_guides.align"),
           file.path(dir.ref, "seqs_V9_LHVR_guides.align"))
names(fastas) = c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9")
fastas
write_shuffles_mix(in_fasta = fastas, N = 100, out_fasta_prefix = "V1V9_mix_guides_", out_dir = file.path(dir.mix, "V1V9_mix_guides/"), guides = TRUE)
write_shuffles(in_fasta = file.path(dir.ref, "seqs_vfull_guides.align"), N = 100, out_fasta_prefix = "VFull_guides_", out_dir = file.path(dir.simulation, "VFull_guides/"))

# Write shuffled data for experimental data
setwd(dir.demo)
write_shuffles(in_fasta = file.path(dir.ref, "experimental_demo_long.align"), N = 100, out_fasta_prefix = "long_", out_dir = file.path(dir.demo, "long"))
write_shuffles(in_fasta = file.path(dir.ref, "experimental_demo_short.align"), N = 100, out_fasta_prefix = "short_", out_dir = file.path(dir.demo, "short"))
write_shuffles(in_fasta = file.path(dir.ref, "experimental_demo_short_guides_few.align"), N = 100, out_fasta_prefix = "short_guides_few_", out_dir = file.path(dir.demo, "short_guides_few"))
write_shuffles(in_fasta = file.path(dir.ref, "experimental_demo_long.align"), N = 100, out_fasta_prefix = "long_control_", out_dir = file.path(dir.demo, "long_control"))
write_shuffles(in_fasta = file.path(dir.ref, "experimental_demo_short_guides_many.align"), N = 100, out_fasta_prefix = "short_guides_many_", out_dir = file.path(dir.demo, "short_guides_many"))
write_shuffles(in_fasta = file.path(dir.ref, "experimental_demo_short_guides_blast_99.align"), N = 100, out_fasta_prefix = "short_guides_blast_99_", out_dir = file.path(dir.demo, "short_guides_blast_99"))
write_shuffles(in_fasta = file.path(dir.ref, "experimental_demo_short_guides_blast_99.align"), N = 100, out_fasta_prefix = "short_guides_optimized_", out_dir = file.path(dir.demo, "short_guides_optimized"))

