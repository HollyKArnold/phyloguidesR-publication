# phyloguidesR-publication
Hello! This git contains all scripts and data to reproduce figures associated with "Guide sequences build better trees and harmonize regional datasets". 

This paper introduces phyloguidesR - a work flow that researchers can use to build more accurate phylogenetic trees from short-read data, as well as to harmonize datasets across amplicon regions. 

The companion PhyloguidesR tutorial site is located [here](https://github.com/HollyKArnold/phyloguidesR). 

## Scripts
If you would like to recreate analyses for the publication or are a reviewer, we recommend starting in the scripts folder, and looking at them in the following order. Indentations denote the case when one script calles a second. 

1. Scripts for generating and validating guide sequences. Guide sequences are first curated, then extrapolated VRs are validated to ensure they match expected patterns of length, and entropy.

- `curate_silva_seed_master.sh`: This script downloads sequences and curates them for use as guides.
  - `curate_silva_seed_remove_eukaryotes.pl`: Removes eukaryotic sequences from guides as they are out of scope.
- `calculate_entropy_V1.R`: Calculate entropy along the 16S gene (Figure 2B)
- `guide_sequence_validation_V1`: Plot expected VR length vs. guide sequence VR length.
  
2. Scripts for generating short-read sequences

- `simulate_sequences_seed_master.sh`: Simulate short-read sequences.
  - `simulate_sequences_constant_legnth.pl`: Simulate short-read sequences of different lengths.
  - `simulate_sequences_hvr.pl`: Simulate short-read sequences by variable region. 

3. Some scripts for generating short-read alignments and trees. Short-read sequences are combined with or without guide sequences. If used in a mixed simulation, reads are mixed from VRs in equal proportions. 

- `align_sequences_V2.R`: Aligns short-read sequences for hypotheses H1 - H4.
  - `align_sequences_functions_V2.R`: Functions for `align_sequences_V2.R`.
- `align_sequences_mix_V1.R`: Aligns short-read sequences for hypothesis H5. 
  - `align_sequences_mix_functions_V1.R`: Functions for `align_sequences_mix_V1.R`.
- `shuffle_sequences_V2.R`: A script to produce 100 randomly shuffled alignments per simulation.
  `shuffled_ssequences_functions_V2.R`: Functions for `shuffle_sequences_V2.R`.
- `maketree_sequences_V2.R`: Makes trees from input alignments. 
  - `make_sequences_functions_V2.R`: Functions for `maketree_sequences_V2.R`. 
  - The following files were used to build trees for all simulations: `tree_commands_mix_V2.sh`, `tree_commands_simulation_V2.sh`, `tree_commands_experimental_demo_V2.sh`, `tree_commands_experimental_demo_large_V2.sh`, `tree_commands_experimental_demo_guides_optimized_V2.sh`
  
4. Next, short-read trees were compared to long-read trees using global, local, and MAST metrics.
- `sim_analysis_master_V5.R`: Compares trees made from short-reads to long read trees for all hypotheses (4A - 4E, and experimental data) . Produces plot Figure 2A VR lengths.
  - `sim_analysis_functions_V2.R`: Functions for `sim_analysis_master_V5.R`. 

4A. Hypothesis 1: Calculate distance between short-read trees and full-length control trees
  
  - `sim_analysis_jrf_random_pairs.R` Calculate JRF distance between short-read trees and full-length control trees.
  - `sim_analysis_nye_random_pairs.R` Calculate Nye distance between short-read trees and full-length control trees.
  - `sim_analysis_msd_random_pairs.R`: Calculate MSD distance between short-read trees and full-length control trees.
  - `sim_analysis_pid_random_pairs.R`: Calculate PID distance between short-read trees and full-length control trees.
  - `sim_analysis_mcd_random_pairs.R`: Calculate MCD distance between short-read trees and full-legnth control trees.
  - `sim_analysis_path_random_pairings.R`: Calculate Path distance between short-read trees and full-legnth control trees.
  - `sim_analysis_nanostructure_random_pairs.R`: Calculates cladal similarity of size 2, 3, ..., 10 of short-read and long-read trees.
  - `sim_analysis_hvr_mast.R`: Calculates MAST of short-read and long-read trees. 

4B. Hypothesis 2: Do guides improve accuracy across short-read trees across different read lengths? 

  - `sim_analysis_bp_jrf_random_pairs.R`: Calculate JRF distance between short-read and full-length control trees across different short-read lengths.
  - `sim_analysis_bp_nye_random_pairs.R`: Calculate Nye distance between short-read and full-length control trees across different short-read lengths.
  - `sim_analysis_bp_msd_random_pairs.R`: Calculate MSD distance between short-read and full-length control trees across different short-read lengths.
  - `sim_analysis_bp_pid_random_pairs.R`: Calculate PID distance between short-read and full-length control trees across different short-read lengths.
  - `sim_analysis_bp_mcd_random_pairs.R`: Calculate MCD distance between short-read and full-length control trees across different short-read lengths.
  - `sim_analysis_bp_path_random_pairs.R`: Calculate PATH distance between short-read and full-length control trees across different short-read lengths.
  - `sim_analysis_bp_nanostructure_random_pairs_[2, ... 10].R`: Calculate similarity of sister pairs, triplets, ..., decuplets between long-read and short-read trees across VR and short-read lengths. 
  - `sim_analysis_bp_mast.R`: Calculate MAST distance between short-read and full-length control trees across different short-read lengths.
  
4C. Hypothesis 3: Do guides allow for integration of overlapping VRs? 

  - `sim_analysis_ov_jrf.R`: Calculate JRF distances between short-read trees compared to long-read trees for mixed overlapping regions and single region controls. 
  - `sim_analysis_ov_nye.R`: Calculate Nye distances between short-read trees compared to long-read trees for mixed overlapping regions and single region controls. 
  - `sim_analysis_ov_msd.R`: Calculate MSD distances between short-read trees compared to long-read trees for mixed overlapping regions and single region controls. 
  - `sim_analysis_ov_pid.R`: Calculate PID distances between short-read trees compared to long-read trees for mixed overlapping regions and single region controls. 
  - `sim_analysis_ov_mcd.R`: Calculate MCD distances between short-read trees compared to long-read trees for mixed overlapping regions and single region controls. 
  - `sim_analysis_ov_path.R`: Calculate Path distances between short-read trees compared to long-read trees for mixed overlapping regions and single region controls. 
  - `sim_analysis_ov_nanostructure_random_pairs_[2, ..., 10].R`: Calculate shared sister pairs, triplets, ..., decuplets between short-read trees compared to long-read trees for mixed overlapping regions and single region controls. 
  - `sim_analysis_ov_mast.R`: Calculate Mast distances between short-read trees compared to long-read trees for mixed overlapping regions and single region controls.


4D. Hypothesis 4: Do guides allow for integration of non-overlapping VRs? 

  - `sim_analysis_nonov_jrf.R`: Calculate JRF distances between short-read trees compared to long-read trees for mixed non-overlapping regions and single region controls. 
  - `sim_analysis_nonov_nye.R`:  Calculate Nye distances between short-read trees compared to long-read trees for mixed non-overlapping regions and single region controls. 
  - `sim_analysis_nonov_msd.R`:  Calculate MSD distances between short-read trees compared to long-read trees for mixed non-overlapping regions and single region controls. 
  - `sim_analysis_nonov_pid.R`:  Calculate PID distances between short-read trees compared to long-read trees for mixed non-overlapping regions and single region controls. 
  - `sim_analysis_nonov_mcd.R`: Calculate MCD distances between short-read trees compared to long-read trees for mixed non-overlapping regions and single region controls. 
  - `sim_analysis_nonov_path.R`:  Calculate Path distances between short-read trees compared to long-read trees for mixed non-overlapping regions and single region controls. 
  - `sim_analysis_nonov_nanostructure_random_pairs_[2, ..., 10].R`:  Calculate sister pairs, triplets, ..., sister decuplets between short-read trees compared to long-read trees for mixed non-overlapping regions and single region controls. 
  - `sim_analysis_nonov_mast.R`:  Calculate MAST distances between short-read trees compared to long-read trees for mixed non-overlapping regions and single region controls. 

4E. Hypothesis 5: Do guides allow for mixing a diverse mix of VRs? 

  - `sim_analysis_nd_jrf.R`: Calculate JRF distances between short-read trees compared to long-read trees for mixed non-overlapping regions and single region controls.
  - `scripts/sim_analysis_nd_nye.R`Calculate Nye distances between short-read trees compared to long-read trees for mixed non-overlapping regions and single region controls.
  - `sim_analysis_nd_msd.R`: Calculate MSD distances between short-read trees compared to long-read trees for mixed non-overlapping regions and single region controls.
  - `sim_analysis_nd_pid.R`: Calculate PID distances between short-read trees compared to long-read trees for mixed non-overlapping regions and single region controls.
  - `sim_analysis_nd_mcd.R`: Calculate MCD distances between short-read trees compared to long-read trees for mixed non-overlapping regions and single region controls.
  - `sim_analysis_nd_path.R`: Calculate path distances between short-read trees compared to long-read trees for mixed non-overlapping regions and single region controls.
  - `sim_analysis_nd_nanostructure_random_pairs[2, .., 10].R`: Calculate percent shared sister pairs, triplets,..., decuplets between short-read trees compared to long-read trees for mixed non-overlapping regions and single region controls.
  - `sim_analysis_nd_mast.R`: Calculate MAST distances between short-read trees compared to long-read trees for mixed non-overlapping regions and single region controls.
  
5. Demo with experimental data

  - `subanalysis_HA_SC_EL.Rmd`: Analysis pathway of quality filtering experimental data raw reads. 
  - `experimental_demo_jrf_random_pairs.R`: Calculate JRF distances between short-read trees compared to long-read trees using experimental data. 
  - `experimental_demo_nye_random_paris.R`: Calculate Nye distances between short-read trees compared to long-read trees using experimental data. 
  - `experimental_demo_msd_random_pairs.R`: Calculate MSD distances between short-read trees compared to long-read trees using experimental data. 
  - `experimental_demo_pid_random_pairs.R`: Calculate PID distances between short-read trees compared to long-read trees using experimental data. 
  - `experimental_demo_mcd_random_pairs.R`: Calculate MCD distances between short-read trees compared to long-read trees using experimental data. 
  - `experimental_demo_path_random_pairings.R`: Calculate Path distances between short-read trees compared to long-read trees using experimental data. 
  - `experimental_demo_nanostructure_random_pairs.R`: Calculate sister pair similarity between short-read trees compared to long-read trees using experimental data. 
  - `experimental_demo_mast_random_pairs.R`: Calculate Mast distances between short-read trees compared to long-read trees using experimental data. 
  
6. Data visualization
  - `sim_analysis_visualization_V2.R`: Visualize results
  - `sim_analysis_visualization_functions.R`: Functions ot visualize reuslts.


## Data
Data used to generate figures is stored in the data/ folder. Here are some brief descriptions of each datastructure.

- Introductory figures:
    - `hvr_by_length.rds`: Length of each HVR for each guide sequence (Figure 2A VR Lengths)
    
- Hypothesis 1: Do guides improve accuracy of short-read trees at each VR? 
  - `hvr_rf_random_pairs.rds`: RF distances of short-read trees compared to long-read trees for VRs. 
  - `hvr_jrf_random_pairs.rds`: JRF distances of short-read trees compared to long-read trees for VRs.
  - `hvr_nye_random_pairs.rds`: Nye distances of short-read trees compared to long-read trees for VRs.
  - `hvr_msd_random_pairs.rds`: MSD distances of short-read trees compared to long-read trees for VRs.
  - `hvr_pid_random_pairs.rds`: PID distances of short-read trees compared to long-read trees for VRs. 
  - `hvr_mcd_random_pairs.rds`: MCD distances of short-read trees compared to long-read trees for VRs.
  - `hvr_path_random_pairings.rds`: Path distance between short-read and long-read trees for VRs.
  - `hvr_quartet_random_pairs.rds`: Quartet distances of short-read trees compared ot long-read trees for VRs. 
  - `hvr_quartet_normalized_random_pairs.rds`: Normalized quartet distances of short-read trees compared to long-read trees for VRs.
  - `hvr_s[2-10]_random_pairs.rds`: Cladal similarity of size 2, 3, ..., 10 of short-read and long-read trees for VRs. 
  - `hvr_mast_random_pairs.rds`: MAST of short-read and long-read trees for VRs. 

- Hypothesis 2: Do guides improve accuracy across short-read trees across different read lengths? 
  - `bp_rf_random_pairs.rds`: RF distances of short-read trees compared to long-read trees for VRs across different lengths.  
  - `bp_jrf_random_pairs.rds`: JRF distances of short-read trees compared to long-read trees for VRs across different lengths.  
  - `bp_nye_random_pairs.rds`: Nye distances of short-read trees compared to long-read trees for VRs across different lengths. 
  - `bp_msd_random_pairs.rds`: MSD distances of short-read trees compared to long-read trees for VRs across different lengths. 
  - `bp_pid_random_pairs.rds`: PID distances of short-read trees compared to long-read trees for VRs across different lengths. 
  - `bp_mcd_random_pairs.rds`: MCD distances of short-read trees compared to long-read trees for VRs across different lengths. 
  - `bp_path_random_pairs.rds`: Path distances of short-read trees compared to long-read trees for VRs across different lengths. 
  - `bp_quartet_normalized_random_pairs.rds`: Normalized quartet distances of short-read trees compared to long-read trees for VRs across different lengths. 
  - `bp_quartet_normalized.rds`: Normalized quartet distances of short-read trees compared to long-read trees for VRs across different lengths. 
  - `nanostructure_bp_[2, ..., 10]_random_pairs.rds`: Calculate similarity of sister pairs, triplets, ..., decuplets between long-read and short-read trees across VR and short-read lengths. 
  - `bp_mast.rds`: Calculate MAST distance between short-read and full-length control trees across different short-read lengths.
  
- Hypothesis 3: Do guides allow for integration of overlapping VRs? 
  - `ov_rf_random_pairs.rds`: RF distances between short-read trees compared to long-read trees for mixed overlapping regions and single region controls.
  - `ov_jrf_random_pairs.rds`: JRF distances between short-read trees compared to long-read trees for mixed overlapping regions and single region controls.
  - `ov_nye_random_pairs.rds`: Nye distances between short-read trees compared to long-read trees for mixed overlapping regions and single region controls.
  - `ov_msd_random_pairs.rds`: MSD distances between short-read trees compared to long-read trees for mixed overlapping regions and single region controls.
  - `ov_pid_random_pairs.rds`: PID distances between short-read trees compared to long-read trees for mixed overlapping regions and single region controls.
  - `ov_mcd_random_pairs.rds`:MCD distances between short-read trees compared to long-read trees for mixed overlapping regions and single region controls.
  - `ov_path_random_pairs.rds`: Path distances between short-read trees compared to long-read trees for mixed overlapping regions and single region controls.
  - `ov_quartet_normalized_random_pairs.rds`: Normalized quartet distances between short-read trees compared to long-read trees for mixed overlapping regions and single region controls.
  - `ov_quartet_random_pairs.rds`: Quartet distances between short-read trees compared to long-read trees for mixed overlapping regions and single region controls.
  - `ov_nanoclusters_random_pairs_[2, .., 10].rds`: Percent shared sister pairs, triplets, ..., decuplets between short-read trees compared to long-read trees for mixed overlapping regions and single region controls. 
  - `ov_mast_random_pairs.rds`: MAST distances between short-read trees compared to long-read trees for mixed overlapping regions and single region controls.

- Hypothesis 4: Do guides allow for integration of non-overlapping VRs? 
  - `nonov_rf_random_pairs.rds`: RF distances between short-read trees compared to long-read trees for mixed non-overlapping regions and single region controls.
  - `nonov_jrf_random_pairs.rds`:  JRF distances between short-read trees compared to long-read trees for mixed non-overlapping regions and single region controls.
  - `nonov_nye_random_pairs.rds`:  Nye distances between short-read trees compared to long-read trees for mixed non-overlapping regions and single region controls.
  - `nonov_msd_random_pairs.rds`:  MSD distances between short-read trees compared to long-read trees for mixed non-overlapping regions and single region controls.
  - `nonov_pid_random_pairs.rds` PID distances between short-read trees compared to long-read trees for mixed non-overlapping regions and single region controls.
  - `nonov_mcd_random_pairs.rds`:  MCD distances between short-read trees compared to long-read trees for mixed non-overlapping regions and single region controls.
  - `nonov_path_random_pairs.rds`:  PATH distances between short-read trees compared to long-read trees for mixed non-overlapping regions and single region controls.
  - `nonov_quartet_random_pairs.rds`:  Quartet distances between short-read trees compared to long-read trees for mixed non-overlapping regions and single region controls.
  - `nonov_quartet_normalized_random_pairs.rds`:  Normalized quartet distances between short-read trees compared to long-read trees for mixed non-overlapping regions and single region controls.
  - `nonov_nanostructure_[2, ..., 10]_random_pairs`: Shared sister pairs, triplets, ..., sister decuplets between short-read trees compared to long-read trees for mixed non-overlapping regions and single region controls. 
  - `nonov_mast_random_pairs.rds`: Normalized MAST distances between short-read trees compared to long-read trees for mixed non-overlapping regions and single region controls.

- Hypothesis 5: Do guides allow for mixing a diverse mix of VRs? 
  - `nd_rf_random_pairs.rds`: RF distances between short-read trees compared to long-read trees for mixed regions and single region controls.
  - `sim_analysis_nd_jrf.rds`: JRF distances between short-read trees compared to long-read trees for mixed regions and single region controls.
  - `sim_analysis_nd_nye.rds`: Nye distances between short-read trees compared to long-read trees for mixed regions and single region controls.
  - `sim_analysis_nd_msd.rds`: MSD distances between short-read trees compared to long-read trees for mixed regions and single region controls.
  - `sim_analysis_nd_pid.rds`: PID distances between short-read trees compared to long-read trees for mixed regions and single region controls.
  - `sim_analysis_nd_mcd.rds`: MCD distances between short-read trees compared to long-read trees for mixed regions and single region controls.
  - `nsim_analysis_nd_path.rds`: Path distances between short-read trees compared to long-read trees for mixed regions and single region controls.
  - `nd_quartet_random_pairs.rds`: Quartet distances between short-read trees compared to long-read trees for mixed regions and single region controls.
  - `nd_quartet_normalized_random_pairs.rds`: Normalized quartet distances between short-read trees compared to long-read trees for mixed regions and single region controls.
  - `nanostructure_nd_[2,...,10]_random_pairs.rds`: Percent shared sister pairs, triplets, ..., decuplets between short-read trees compared to long-read trees for mixed overlapping regions and single region controls. 
  - `sim_analysis_nd_mast.R`: MAST distances between short-read trees compared to long-read trees for mixed regions and single region controls.

- Experimental Data Demo
  - `experimental_demo_ps_long_reads.rds`: Long read experimental data. 
  - `experimental_demo_ps_shortreads_12fr.rds`: Short read experimental data. 
   - `experimental_demo_rf_random_pairs.rds`: RF distances between short-read trees compared to long-read trees using experimental data. 
   - `experimental_demo_jrf_random_pairs.rds`: JRF distances between short-read trees compared to long-read trees using experimental data. 
   - `experimental_demo_nye_random_pairs.rds`: Nye distances between short-read trees compared to long-read trees using experimental data. 
   - `experimental_demo_msd_random_pairs.rds`: MSD distances between short-read trees compared to long-read trees using experimental data. 
   - `experimental_demo_pid_random_pairs.rds`: PID distances between short-read trees compared to long-read trees using experimental data. 
   - `experimental_demo_mcd_random_pairs.rds`: MCD distances between short-read trees compared to long-read trees using experimental data. 
   - `experimental_demo_path_random_pairings.rds`: Path distances between short-read trees compared to long-read trees using experimental data. 
   - `experimental_demo_quartet.rds`: Quartet distances between short-read trees compared to long-read trees using experimental data. 
   - `experimental_demo_quartet_normalized.rds`: Normalized quartet distances between short-read trees compared to long-read trees using experimental data. 
   - `experimental_demo_subtrees_random_pairs.rds`: Similarity of sister pairs between short-read trees compared to long-read trees using experimental data. 
   - `experimental_demo_mast_random_pairs.rds`: MAST distances between short-read trees compared to long-read trees using experimental data. 

