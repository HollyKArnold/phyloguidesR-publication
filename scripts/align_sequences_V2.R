
#AUTHOR: ARNOLD
#DAY: Jan 13th, 2025
# align_sequences_V2.R
# The purpose of this file is to align sequences with and without guides.

# LIBRARIES
.libPaths(c(.libPaths(), "/home/micro/arnoldho/R/x86_64-pc-linux-gnu-library/4.1"))
library(phylotools)
library(dplyr)
library("PIMMA")

# DIRECTORIES
HOME = "/nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration"
setwd(HOME)
dir.scripts = file.path(HOME, "2025_hvr_guide_phylogenetic_integration/scripts/")
dir.ref = file.path(HOME, "reference_data/")

TEMPLATE_SEED = file.path(dir.ref, "ref.align")
if(!file.exists(TEMPLATE_SEED)){
  print("FAIL: Cannot find template seed file.")
}else{print("PASS: Can find the seed file.")}

# Load scripts
setwd(dir.scripts)
source("align_sequences_functions_V2.R")

setwd(dir.ref)

# Align V1 - V9 HVR
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V1_LHVR.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V1_LHVR_guides.fasta")

mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V2_LHVR.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V2_LHVR_guides.fasta")

mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V3_LHVR.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V3_LHVR_guides.fasta")

mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V4_LHVR.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V4_LHVR_guides.fasta")

mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V5_LHVR.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V5_LHVR_guides.fasta")

mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V6_LHVR.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V6_LHVR_guides.fasta")

mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V7_LHVR.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V7_LHVR_guides.fasta")

mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V8_LHVR.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V8_LHVR_guides.fasta")

mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V9_LHVR.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V9_LHVR_guides.fasta")


# Align V1 across length with and without guides
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V1_L50.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V1_L50_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V1_L100.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V1_L100_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V1_L200.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V1_L200_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V1_L300.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V1_L300_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V1_L400.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V1_L400_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V1_L500.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V1_L500_guides.fasta")

# Align V2 across length with and without guides
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V2_L50.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V2_L50_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V2_L100.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V2_L100_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V2_L200.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V2_L200_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V2_L300.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V2_L300_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V2_L400.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V2_L400_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V2_L500.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V2_L500_guides.fasta")


# Align V3 across length with and without guides
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V3_L50.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V3_L50_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V3_L100.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V3_L100_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V3_L200.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V3_L200_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V3_L300.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V3_L300_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V3_L400.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V3_L400_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V3_L500.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V3_L500_guides.fasta")

# Align V4 across length with and without guides
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V4_L50.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V4_L50_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V4_L100.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V4_L100_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V4_L200.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V4_L200_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V4_L300.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V4_L300_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V4_L400.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V4_L400_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V4_L500.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V4_L500_guides.fasta")

# Align V5 across length with and without guides
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V5_L50.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V5_L50_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V5_L100.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V5_L100_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V5_L200.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V5_L200_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V5_L300.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V5_L300_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V5_L400.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V5_L400_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V5_L500.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V5_L500_guides.fasta")

# Align V6 across length with and without guides
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V6_L50.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V6_L50_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V6_L100.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V6_L100_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V6_L200.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V6_L200_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V6_L300.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V6_L300_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V6_L400.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V6_L400_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V6_L500.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V6_L500_guides.fasta")

# Align V7 across length with and without guides
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V7_L50.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V7_L50_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V7_L100.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V7_L100_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V7_L200.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V7_L200_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V7_L300.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V7_L300_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V7_L400.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V7_L400_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V7_L500.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V7_L500_guides.fasta")

# Align V8 across length with and without guides
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V8_L50.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V8_L50_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V8_L100.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V8_L100_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V8_L200.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V8_L200_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V8_L300.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V8_L300_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V8_L400.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V8_L400_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V8_L500.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V8_L500_guides.fasta")

# Align V9 across length with and without guides
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V9_L50.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V9_L50_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V9_L100.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V9_L100_guides.fasta")

# Align overlapping, disjoint, and VRs used for mixed analysis.
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V1V2_LHVR.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V1V2_LHVR_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V1V3_LHVR.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V1V3_LHVR_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V3V4_LHVR.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V3V4_LHVR_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V3V5_LHVR.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V3V5_LHVR_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V4V5_LHVR.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V4V5_LHVR_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V5V7_LHVR.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V5V7_LHVR_guides.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V7V9_LHVR.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_V7V9_LHVR_guides.fasta")

# Control: Full length sequences with themselves.
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "seqs_vfull_guides.fasta")

# Experimental data: Align sequences with and without guides. 
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "experimental_demo_short.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "experimental_demo_short_guides_few.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "experimental_demo_short_guides_many.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "experimental_demo_long.fasta")
mothur_align_seqs(template = TEMPLATE_SEED, output.directory = dir.ref, candidate = "experimental_demo_short_guides_blast_99.fasta")
