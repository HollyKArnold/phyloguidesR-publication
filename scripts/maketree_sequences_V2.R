#AUTHOR: ARNOLD
#DAY: Jan 13th, 2025
# maketree_sequences_functions.R
# Makes trees for each input aligments.

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
source("maketree_sequences_functions_V2.R")

# Get paths to all simulated files for which to align.
setwd(dir.mix)
directories <- list.dirs(dir.mix, full.names = TRUE, recursive = FALSE)
directories
basename(directories)
cmds = vector()
for(i in 1:length(directories)){
  print(basename(directories[i]))
  cmds = c(make_fasttrees(in_dir = directories[i], file_prefix = paste0(c(basename(directories[i]), "_"), sep = "", collapse = ""), out_dir = directories[i], N = 100), cmds)
}

# setwd(dir.scripts)
# sink("tree_commands_mix_V2.sh")
# cat(cmds, sep = "\n")
# sink()

setwd(dir.simulation)
directories <- list.dirs(dir.simulation, full.names = TRUE, recursive = FALSE)
directories
basename(directories)
cmds = vector()
for(i in 1:length(directories)){
  print(basename(directories[i]))
  cmds = c(make_fasttrees(in_dir = directories[i], file_prefix = paste0(c(basename(directories[i]), "_"), sep = "", collapse = ""), out_dir = directories[i], N = 100), cmds)
}

setwd(dir.scripts)
sink("tree_commands_simulation_V2.sh")
cat(cmds, sep = "\n")
sink()

# Get paths to all simulated files for which to align.
setwd(dir.demo)
directories <- list.dirs(dir.demo, full.names = TRUE, recursive = FALSE)
directories
basename(directories)
cmds = vector()
for(i in 1:length(directories)){
  print(basename(directories[i]))
  cmds = c(cmds, make_fasttrees(in_dir = directories[i], file_prefix = paste0(c(basename(directories[i]), "_"), sep = "", collapse = ""), out_dir = directories[i], N = 100))
}

#setwd(dir.scripts)
#sink("tree_commands_experimental_demo_V2.sh")
#cat(cmds, sep = "\n")
#sink()

# Increase memory for the large guide file
setwd(dir.demo)
directories <- list.dirs(dir.demo, full.names = TRUE, recursive = FALSE)[5]
directories
basename(directories)
cmds = vector()
for(i in 1:length(directories)){
  print(basename(directories[i]))
  cmds = c(cmds, make_fasttrees_large(in_dir = directories[i], file_prefix = paste0(c(basename(directories[i]), "_"), sep = "", collapse = ""), out_dir = directories[i], N = 100))
}

setwd(dir.scripts)
sink("tree_commands_experimental_demo_large_V2.sh")
cat(cmds, sep = "\n")
sink()


# Make trees with optimized guides.
setwd(dir.demo)
directories = list.dirs(dir.demo, full.names = TRUE, recursive = FALSE)
directories <- list.dirs(dir.demo, full.names = TRUE, recursive = FALSE)[grepl("short_guides_optimized", directories)]
directories
basename(directories)
cmds = vector()
for(i in 1:length(directories)){
  print(basename(directories[i]))
  cmds = c(cmds, make_fasttrees_ruby(in_dir = directories[i], file_prefix = paste0(c(basename(directories[i]), "_"), sep = "", collapse = ""), out_dir = directories[i], N = 100))
}

setwd(dir.scripts)
sink("tree_commands_experimental_demo_guides_optimized_V2.sh")
cat(cmds, sep = "\n")
sink()

