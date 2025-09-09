################################################################################
#AUTHOR: ARNOLD
#DAY: Jan 24th, 2025
#SCRIPT: sim_analysis_visalization.R
#DESCRIPTION: Companion script to sim_analysis.R
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
library(tidyr)
library(tibble)
library(purrr)
library(ggpubr)
library(reshape2)
library(ggplot2)
library(dplyr)
library(patchwork)
library(stringr)
library(rstatix)
TREE_SIZE = 6846
################################################################################
# SET DIRECTORIES
################################################################################
# Set Home
#HOME = "/nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration"
#HOME_TMP = "/nfs3/Sharpton_Lab/tmp/projects/arnoldh/2025_hvr_guide_phylogenetic_integration/"
HOME = "/Users/arnoldhk/Desktop/Research/2025_HVR_Guide_Phylogenetic_Integration/"

# Make directories
setwd(HOME)
dir.scripts = file.path(HOME, "2025_hvr_guide_phylogenetic_integration/scripts/")
dir.ref = file.path(HOME, "reference_data/")
dir.out = file.path(HOME, "2025_hvr_guide_phylogenetic_integration/data_visulization/")
dir.sim2 = file.path(HOME, "2025_hvr_guide_phylogenetic_integration/out/sim_bp/")
dir.figures = file.path(HOME, "2025_hvr_guide_phylogenetic_integration/figuresNoSync/")
dir.sim3 = file.path(HOME, "2025_hvr_guide_phylogenetic_integration/out/sim_overlapping/")
dir.sim4 = file.path(HOME, "2025_hvr_guide_phylogenetic_integration/out/sim_nonoverlapping/")
dir.sim5 = file.path(HOME, "2025_hvr_guide_phylogenetic_integration/out/sim_nd_mix/")
  
# Get functions
setwd(dir.scripts)
source("sim_analysis_visualization_functions.R")

################################################################################
# COLORS
################################################################################
# Define a function to convert RGB to HEX
rgb_to_hex <- function(r, g, b) {
  rgb(r, g, b, maxColorValue = 255)
}

# Example usage
col_R1 = rgb_to_hex(145, 28, 67) 
col_R2 = rgb_to_hex(227, 117, 79) 
col_R3 = rgb_to_hex(249, 225, 150) 
col_R4 = rgb_to_hex(220, 48, 126) 
col_R5 = rgb_to_hex(180, 220, 168) 
col_R6 = rgb_to_hex(125, 192, 166)
col_R7 =rgb_to_hex(141, 213, 251)
col_R8 = rgb_to_hex(75, 134, 184)
col_R9 = rgb_to_hex(92, 80, 157)

regions <- data.frame(
  start = c(0, 77, 314, 488, 760, 881, 1020, 1149, 1363),           # Start of regions (x-axis)
  end = c(76, 313, 487, 759, 880, 1019, 1148, 1362, 1464),            # End of regions (x-axis)
  label = c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9"), # Labels
  color = c(col_R1, col_R2, col_R3, col_R4, col_R5, col_R6, col_R7, col_R8, col_R9) # Colors for regions
)

# Custom color palette
custom_palette <- setNames(regions$color, regions$label)
pal_no_guides = "#003F5C"
pal_guides= "#FFA600"

################################################################################
# PLOT HVR vs Length
################################################################################
setwd(dir.out)
hvrs = readRDS("hvr_by_length.rds")
hvrs

hvr_lengths_boxplot = ggplot(hvrs, aes(x = HVR, y = Length, color = HVR)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "",
       x = "Variable Region",
       y = "Base Pairs") +
  scale_color_manual(values = custom_palette) + # Apply custom colors
  theme(legend.position = "none", axis.text.x = element_text(size = 12)) + 
  theme(axis.text.x = element_text(size = 22, face = "bold"),
        axis.title.x = element_text(size = 22, face = "bold"),
        axis.title.y = element_text(size = 22, face = "bold"))  # Increase X-axis label size

hvr_lengths_boxplot

setwd(dir.figures)
ggsave(filename = "hvr_length.pdf", hvr_lengths_boxplot, width = 1000, height = 1000, units = "px")


################################################################################
# SIM 1: 
# Goal plot global measures of diversity for each HVR vs the other.
################################################################################

################################################################################
# Global Distances: 100 Random pairwise comparisons 
# Determine if guide sequence improve topological accuracy at each HVR
################################################################################

# Read in the distance metrics
setwd(dir.out)

## 1 RF
hvr_rf = as_tibble(readRDS("hvr_rf_random_pairs.rds")) %>% print(Inf)
labels = colnames(hvr_rf)

## 2 JRF
hvr_jrf = as_tibble(readRDS("hvr_jrf_random_pairs.rds")) %>% print(Inf)
if(setequal(colnames(hvr_jrf), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 3 Nye Distance
hvr_nye = as_tibble(readRDS("hvr_nye_random_pairs.rds")) %>% print(Inf)
if(setequal(colnames(hvr_nye), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 4 MSD
hvr_msd = as_tibble(readRDS("hvr_msd_random_pairs.rds")) %>% print(Inf)
if(setequal(colnames(hvr_msd), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 5 PID
hvr_pid = as_tibble(readRDS("hvr_pid_random_pairs.rds")) %>% print(Inf)
if(setequal(colnames(hvr_pid), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 6 MCD
hvr_mcd = as_tibble(readRDS("hvr_mcd_random_pairs.rds")) %>% print(Inf)
if(setequal(colnames(hvr_mcd), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 7 Path
hvr_path = as_tibble(readRDS("hvr_path_random_pairings.rds")) %>% print(Inf)
if(setequal(colnames(hvr_path), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 8 Quartet Distance
hvr_quartet = as_tibble(readRDS("hvr_quartet_normalized_random_pairs.rds")) %>% print(Inf)
if(setequal(colnames(hvr_quartet), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}


# test if there was a significant difference between guides and no guides for all. 
wilcox_test_hvr(hvr_rf)# All are significantly different
wilcox_test_hvr(hvr_jrf)# All are significantly different
wilcox_test_hvr(hvr_nye)# All are significantly different
wilcox_test_hvr(hvr_msd)# All are significantly different
wilcox_test_hvr(hvr_pid)# All are significantly different
wilcox_test_hvr(hvr_mcd)# All are significantly different
wilcox_test_hvr(hvr_path)# All are significantly different
wilcox_test_hvr(hvr_quartet) # All are significantly different aside from V6

## Viz Distances into six different plots
## 1 RF
# Reshape data to long format
hvr_rf_long <- hvr_rf %>%
  select(!VFull_guides) %>%
  rename("C" = VFull) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
  mutate(Group = ifelse(stringr::str_detect(Variable, "_guides"), "With_Guides", "Without_Guides")) %>%
  mutate(Variable = stringr::str_replace(Variable, "_guides", "")) %>% 
  mutate(Group = factor(case_when(Group == "Without_Guides" ~ "No Guides",
                                  Group == "With_Guides" ~ "Guides",
                                  TRUE ~ Group),
                        levels = c("No Guides", "Guides")))
hvr_rf_long

p_rf = ggplot(hvr_rf_long, aes(x = Variable, y = Value, color = Group)) +
  geom_boxplot() +
  labs(x = "Variable", y = "Value", fill = "Condition") +
  theme_minimal() +
  stat_compare_means(aes(group = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE) + 
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "RF", color = "Guides") +  
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 22)) + 
  coord_cartesian(clip = "off")
p_rf

## 2 JRF
hvr_jrf_long <- hvr_jrf %>%
  select(!VFull_guides) %>%
  rename("C" = VFull) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
  mutate(Group = ifelse(stringr::str_detect(Variable, "_guides"), "With_Guides", "Without_Guides"),
         Variable = stringr::str_replace(Variable, "_guides", "")) %>%
  mutate(Group = factor(case_when(Group == "Without_Guides" ~ "No Guides",
                                  Group == "With_Guides" ~ "Guides",
                                  TRUE ~ Group),
                        levels = c("No Guides", "Guides")))


p_jrf = ggplot(hvr_jrf_long, aes(x = Variable, y = Value, color = Group)) +
  geom_boxplot() +
  labs(x = "Variable", y = "Value", fill = "Condition") +
  theme_minimal() +
  stat_compare_means(aes(group = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE)+ 
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "JRF", color = "Guides") + 
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 22)) + 
  coord_cartesian(clip = "off")
p_jrf


## 3 NYE
hvr_nye_long <- hvr_nye %>%
  select(!VFull_guides) %>%
  rename("C" = VFull) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
  mutate(Group = ifelse(stringr::str_detect(Variable, "_guides"), "With_Guides", "Without_Guides"),
         Variable = stringr::str_replace(Variable, "_guides", "")) %>%
  mutate(Group = factor(case_when(Group == "Without_Guides" ~ "No Guides",
                                  Group == "With_Guides" ~ "Guides",
                                  TRUE ~ Group),
                        levels = c("No Guides", "Guides")))


p_nye = ggplot(hvr_nye_long, aes(x = Variable, y = Value, color = Group)) +
  geom_boxplot() +
  labs(x = "Variable", y = "Value", fill = "Condition") +
  theme_minimal() +
  stat_compare_means(aes(group = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE)+ 
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "NYE", color = "Guides") +  # Change the Y-axis label here
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 22)) + # Increase X-axis label size
  coord_cartesian(clip = "off")
p_nye

## 4 MSD
hvr_msd_long <- hvr_msd %>%
  select(!VFull_guides) %>%
  rename("C" = VFull) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
  mutate(Group = ifelse(stringr::str_detect(Variable, "_guides"), "With_Guides", "Without_Guides"),
         Variable = stringr::str_replace(Variable, "_guides", "")) %>%
  mutate(Group = factor(case_when(Group == "Without_Guides" ~ "No Guides",
                                  Group == "With_Guides" ~ "Guides",
                                  TRUE ~ Group),
                        levels = c("No Guides", "Guides")))


p_msd = ggplot(hvr_msd_long, aes(x = Variable, y = Value, color = Group)) +
  geom_boxplot() +
  labs(x = "Variable", y = "Value", fill = "Condition") +
  theme_minimal() +
  stat_compare_means(aes(group = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE)+ 
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "MSD", color = "Guides") +  # Change the Y-axis label here
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 22))+  # Increase X-axis label size
  coord_cartesian(clip = "off")
p_msd


## 5 PID
hvr_pid_long <- hvr_pid %>%
  select(!VFull_guides) %>%
  rename("C" = VFull) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
  mutate(Group = ifelse(stringr::str_detect(Variable, "_guides"), "With_Guides", "Without_Guides"),
         Variable = stringr::str_replace(Variable, "_guides", "")) %>%
  mutate(Group = factor(case_when(Group == "Without_Guides" ~ "No Guides",
                                  Group == "With_Guides" ~ "Guides",
                                  TRUE ~ Group),
                        levels = c("No Guides", "Guides")))


p_pid = ggplot(hvr_pid_long, aes(x = Variable, y = Value, color = Group)) +
  geom_boxplot() +
  labs(x = "Variable", y = "Value", fill = "Condition") +
  theme_minimal() +
  stat_compare_means(aes(group = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE)+ 
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "PID", color = "Guides") +  # Change the Y-axis label here
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 22)) + # Increase X-axis label size
  coord_cartesian(clip = "off")
p_pid

## 6 MCD
hvr_mcd_long <- hvr_mcd %>%
  select(!VFull_guides) %>%
  rename("C" = VFull) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
  mutate(Group = ifelse(stringr::str_detect(Variable, "_guides"), "With_Guides", "Without_Guides"),
         Variable = stringr::str_replace(Variable, "_guides", "")) %>%
  mutate(Group = factor(case_when(Group == "Without_Guides" ~ "No Guides",
                                  Group == "With_Guides" ~ "Guides",
                                  TRUE ~ Group),
                        levels = c("No Guides", "Guides")))


p_mcd = ggplot(hvr_mcd_long, aes(x = Variable, y = Value, color = Group)) +
  geom_boxplot() +
  labs(x = "Variable", y = "Value", fill = "Condition") +
  theme_minimal() +
  stat_compare_means(aes(group = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE)+ 
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "MCD", color = "Guides") +  # Change the Y-axis label here
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 22)) + # Increase X-axis label size
  coord_cartesian(clip = "off")
p_mcd


## 7 PATH
hvr_path_long <- hvr_path %>%
  select(!VFull_guides) %>%
  rename("C" = VFull) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
  mutate(Group = ifelse(stringr::str_detect(Variable, "_guides"), "With_Guides", "Without_Guides"),
         Variable = stringr::str_replace(Variable, "_guides", "")) %>%
  mutate(Group = factor(case_when(Group == "Without_Guides" ~ "No Guides",
                                  Group == "With_Guides" ~ "Guides",
                                  TRUE ~ Group),
                        levels = c("No Guides", "Guides")))


p_path = ggplot(hvr_path_long, aes(x = Variable, y = Value, color = Group)) +
  geom_boxplot() +
  labs(x = "Variable", y = "Value", fill = "Condition") +
  theme_minimal() +
  stat_compare_means(aes(group = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE)+ 
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "PATH", color = "Guides") +  # Change the Y-axis label here
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 22)) + # Increase X-axis label size
  coord_cartesian(clip = "off")
p_path


## 8 QUAR
hvr_quar_long <- hvr_quartet %>%
  select(!VFull_guides) %>%
  rename("C" = VFull) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
  mutate(Group = ifelse(stringr::str_detect(Variable, "_guides"), "With_Guides", "Without_Guides"),
         Variable = stringr::str_replace(Variable, "_guides", "")) %>%
  mutate(Group = factor(case_when(Group == "Without_Guides" ~ "No Guides",
                                  Group == "With_Guides" ~ "Guides",
                                  TRUE ~ Group),
                        levels = c("No Guides", "Guides")))


p_quar = ggplot(hvr_quar_long, aes(x = Variable, y = Value, color = Group)) +
  geom_boxplot() +
  labs(x = "Variable", y = "Value", fill = "Condition") +
  theme_minimal() +
  stat_compare_means(aes(group = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE)+ 
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "QUAR", color = "Guides") +  # Change the Y-axis label here
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 22)) +
  coord_cartesian(clip = "off")
p_quar

layout  = 
  p_rf / p_jrf / p_nye / p_msd / p_pid / p_mcd / p_path / p_quar +
  plot_annotation(title = "Topological Distances",
                  theme = theme(plot.title = element_text(size = 16, hjust = 0.5)))
layout

setwd(dir.figures)
ggsave("hvr_global_dists.pdf", plot = layout, height = 400, width = 300, units = "mm")
ggsave("hvr_global_dists_RF_only.pdf", plot = p_rf, height = 150, width = 150, units = "mm")

################################################################################
# Global Distances: 100 Random pairwise comparisons 
# Determine which HVR performed the best for each topological comparison compared 
# to control
################################################################################
# 1 RF
box_rf = get_boxplot_global_dist_ordered(df = hvr_rf, ylabel = "RF")
heat_rf = get_heatmap_pairwise_global_dist(df = hvr_rf, title = "RF")

# 2 JRF
box_jrf = get_boxplot_global_dist_ordered(df = hvr_jrf, ylabel = "JRF")
heat_jrf = get_heatmap_pairwise_global_dist(df = hvr_rf, title = "JRF")

# 3 NYE
box_nye = get_boxplot_global_dist_ordered(df = hvr_nye, ylabel = "NYE")
heat_nye = get_heatmap_pairwise_global_dist(df = hvr_nye, title = "NYE")

# 4 MSD
box_msd = get_boxplot_global_dist_ordered(df = hvr_msd, ylabel = "MSD")
heat_msd = get_heatmap_pairwise_global_dist(df = hvr_msd, title = "MSD")

# 5 PID
box_pid = get_boxplot_global_dist_ordered(df = hvr_pid, ylabel = "PID")
heat_pid = get_heatmap_pairwise_global_dist(df = hvr_pid, title = "PID")

# 6 MCD
box_mcd = get_boxplot_global_dist_ordered(df = hvr_mcd, ylabel = "MCD")
heat_mcd = get_heatmap_pairwise_global_dist(df = hvr_mcd, title = "MCD")

# 7 PD
box_pd = get_boxplot_global_dist_ordered(df = hvr_path, ylabel = "PATH")
heat_pd = get_heatmap_pairwise_global_dist(df = hvr_path, title = "PATH")

# 8 QUAR
box_quar = get_boxplot_global_dist_ordered(df = hvr_quartet, ylabel = "QUAR")
heat_quar = get_heatmap_pairwise_global_dist(df = hvr_quartet, title = "QUAR")

# Viz
# Create titles for the two columns
title_accuracy <- ggplot() + 
  annotate("text", x = 1, y = 1, label = "Ordered Global Accuracy", size = 6, fontface = "bold") +
  theme_void()

title_heatmaps <- ggplot() + 
  annotate("text", x = 1, y = 1, label = "Pairwise Accuracy (Q)", size = 6, fontface = "bold") +
  theme_void()

# Arrange plots in a 2x8 layout
final_plot_A <- (title_accuracy | title_heatmaps) /  # Column titles
  (box_rf  | heat_rf)  /
  (box_jrf | heat_jrf) /
  (box_nye | heat_nye) /
  (box_msd | heat_msd) 
final_plot_B <- (title_accuracy | title_heatmaps) /  # Column titles
  (box_pid | heat_pid) /
  (box_mcd | heat_mcd) /
  (box_pd  | heat_pd)  /
  (box_quar | heat_quar)
final_plot_A
final_plot_B

# Display the final plot
final_plot <- (title_accuracy | title_heatmaps) /
  (box_rf  | heat_rf)  /
  (box_jrf | heat_jrf) /
  (box_nye | heat_nye) /
  (box_msd | heat_msd) / 
  (box_pid | heat_pid) /
  (box_mcd | heat_mcd) /
  (box_pd  | heat_pd)  /
  (box_quar | heat_quar)
setwd(dir.figures)
ggsave(filename = "hvr_order_and_significance_comparisons.pdf", final_plot, height = 600, width  = 225, units = "mm")
#ggsave(filename = "hvr_order_and_significance_comparisons_A.pdf", final_plot_A, height = 400, width  = 225, units = "mm")
#ggsave(filename = "hvr_order_and_significance_comparisons_B.pdf", final_plot_B, height = 400, width  = 225, units = "mm")
ggsave("hvr_global_dists_ordered_RF_only.pdf", plot = box_rf, height = 150, width = 150, units = "mm")

################################################################################
# Nano Distances: 100 Random pairwise comparisons 
################################################################################
setwd(dir.out)
s2 = readRDS("hvr_s2_random_pairs.rds")
s3 = readRDS("hvr_s3_random_pairs.rds")
s4 = readRDS("hvr_s4_random_pairs.rds")
s5 = readRDS("hvr_s5_random_pairs.rds")
s6 = readRDS("hvr_s6_random_pairs.rds")
s7 = readRDS("hvr_s7_random_pairs.rds")
s8 = readRDS("hvr_s8_random_pairs.rds")
s9 = readRDS("hvr_s9_random_pairs.rds")
s10 = readRDS("hvr_s10_random_pairs.rds")

## S2
s2_long <- s2 %>%
  select(!C_guides) %>%
  rename("C" = C) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
  mutate(Group = ifelse(stringr::str_detect(Variable, "_guides"), "With_Guides", "Without_Guides"),
         Variable = stringr::str_replace(Variable, "_guides", "")) %>%
  mutate(Group = factor(case_when(Group == "Without_Guides" ~ "No Guides",
                                  Group == "With_Guides" ~ "Guides",
                                  TRUE ~ Group),
                        levels = c("No Guides", "Guides")))


p_s2 = ggplot(s2_long, aes(x = Variable, y = Value, color = Group)) +
  geom_boxplot() +
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.7, dodge.width = 0.8), alpha = 0.001, size = .005) +  # Jittered points
  labs(x = "Variable", y = "Value", fill = "Condition") +
  theme_minimal() +
  stat_compare_means(aes(group = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE)+ 
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "Nano 2", color = "Guides") +  # Change the Y-axis label here
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 22)) + # Increase X-axis label size
  coord_cartesian(clip = "off")
p_s2

## S3
s3_long <- s3 %>%
  select(!C_guides) %>%
  rename("C" = C) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
  mutate(Group = ifelse(stringr::str_detect(Variable, "_guides"), "With_Guides", "Without_Guides"),
         Variable = stringr::str_replace(Variable, "_guides", "")) %>%
  mutate(Group = factor(case_when(Group == "Without_Guides" ~ "No Guides",
                                  Group == "With_Guides" ~ "Guides",
                                  TRUE ~ Group),
                        levels = c("No Guides", "Guides")))


p_s3 = ggplot(s3_long, aes(x = Variable, y = Value, color = Group)) +
  geom_boxplot() +
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.7, dodge.width = 0.8), alpha = 0.001, size = .005) +  # Jittered points
  labs(x = "Variable", y = "Value", fill = "Condition") +
  theme_minimal() +
  stat_compare_means(aes(group = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE)+ 
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "Nano 3", color = "Guides") +  # Change the Y-axis label here
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 22)) + # Increase X-axis label size
  coord_cartesian(clip = "off")
p_s3

## S4
s4_long <- s4 %>%
  select(!C_guides) %>%
  rename("C" = C) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
  mutate(Group = ifelse(stringr::str_detect(Variable, "_guides"), "With_Guides", "Without_Guides"),
         Variable = stringr::str_replace(Variable, "_guides", "")) %>%
  mutate(Group = factor(case_when(Group == "Without_Guides" ~ "No Guides",
                                  Group == "With_Guides" ~ "Guides",
                                  TRUE ~ Group),
                        levels = c("No Guides", "Guides")))


p_s4 = ggplot(s4_long, aes(x = Variable, y = Value, color = Group)) +
  geom_boxplot() +
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.7, dodge.width = 0.8), alpha = 0.001, size = .005) +  # Jittered points
  labs(x = "Variable", y = "Value", fill = "Condition") +
  theme_minimal() +
  stat_compare_means(aes(group = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE)+ 
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "Nano 4", color = "Guides") +  # Change the Y-axis label here
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 22)) + # Increase X-axis label size
  coord_cartesian(clip = "off")
p_s4

## S5
s5_long <- s5 %>%
  select(!C_guides) %>%
  rename("C" = C) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
  mutate(Group = ifelse(stringr::str_detect(Variable, "_guides"), "With_Guides", "Without_Guides"),
         Variable = stringr::str_replace(Variable, "_guides", "")) %>%
  mutate(Group = factor(case_when(Group == "Without_Guides" ~ "No Guides",
                                  Group == "With_Guides" ~ "Guides",
                                  TRUE ~ Group),
                        levels = c("No Guides", "Guides")))


p_s5 = ggplot(s5_long, aes(x = Variable, y = Value, color = Group)) +
  geom_boxplot() +
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.7, dodge.width = 0.8), alpha = 0.001, size = .005) +  # Jittered points
  labs(x = "Variable", y = "Value", fill = "Condition") +
  theme_minimal() +
  stat_compare_means(aes(group = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE)+ 
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "Nano 5", color = "Guides") +  # Change the Y-axis label here
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 22)) + # Increase X-axis label size
  coord_cartesian(clip = "off")
p_s5

## S6
s6_long <- s6 %>%
  select(!C_guides) %>%
  rename("C" = C) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
  mutate(Group = ifelse(stringr::str_detect(Variable, "_guides"), "With_Guides", "Without_Guides"),
         Variable = stringr::str_replace(Variable, "_guides", "")) %>%
  mutate(Group = factor(case_when(Group == "Without_Guides" ~ "No Guides",
                                  Group == "With_Guides" ~ "Guides",
                                  TRUE ~ Group),
                        levels = c("No Guides", "Guides")))


p_s6 = ggplot(s6_long, aes(x = Variable, y = Value, color = Group)) +
  geom_boxplot() +
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.7, dodge.width = 0.8), alpha = 0.001, size = .005) +  # Jittered points
  labs(x = "Variable", y = "Value", fill = "Condition") +
  theme_minimal() +
  stat_compare_means(aes(group = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE)+ 
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "Nano 6", color = "Guides") +  # Change the Y-axis label here
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 22)) + # Increase X-axis label size
  coord_cartesian(clip = "off")
p_s6

## S7
s7_long <- s7 %>%
  select(!C_guides) %>%
  rename("C" = C) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
  mutate(Group = ifelse(stringr::str_detect(Variable, "_guides"), "With_Guides", "Without_Guides"),
         Variable = stringr::str_replace(Variable, "_guides", "")) %>%
  mutate(Group = factor(case_when(Group == "Without_Guides" ~ "No Guides",
                                  Group == "With_Guides" ~ "Guides",
                                  TRUE ~ Group),
                        levels = c("No Guides", "Guides")))


p_s7 = ggplot(s7_long, aes(x = Variable, y = Value, color = Group)) +
  geom_boxplot() +
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.7, dodge.width = 0.8), alpha = 0.001, size = .005) +  # Jittered points
  labs(x = "Variable", y = "Value", fill = "Condition") +
  theme_minimal() +
  stat_compare_means(aes(group = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE)+ 
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "Nano 7", color = "Guides") +  # Change the Y-axis label here
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 22)) + # Increase X-axis label size
  coord_cartesian(clip = "off")
p_s7

## S8
s8_long <- s8 %>%
  select(!C_guides) %>%
  rename("C" = C) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
  mutate(Group = ifelse(stringr::str_detect(Variable, "_guides"), "With_Guides", "Without_Guides"),
         Variable = stringr::str_replace(Variable, "_guides", "")) %>%
  mutate(Group = factor(case_when(Group == "Without_Guides" ~ "No Guides",
                                  Group == "With_Guides" ~ "Guides",
                                  TRUE ~ Group),
                        levels = c("No Guides", "Guides")))


p_s8 = ggplot(s8_long, aes(x = Variable, y = Value, color = Group)) +
  geom_boxplot() +
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.7, dodge.width = 0.8), alpha = 0.001, size = .005) +  # Jittered points
  labs(x = "Variable", y = "Value", fill = "Condition") +
  theme_minimal() +
  stat_compare_means(aes(group = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE)+ 
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "Nano 8", color = "Guides") +  # Change the Y-axis label here
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 22)) + # Increase X-axis label size
  coord_cartesian(clip = "off")
p_s8

## S9
s9_long <- s9 %>%
  select(!C_guides) %>%
  rename("C" = C) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
  mutate(Group = ifelse(stringr::str_detect(Variable, "_guides"), "With_Guides", "Without_Guides"),
         Variable = stringr::str_replace(Variable, "_guides", "")) %>%
  mutate(Group = factor(case_when(Group == "Without_Guides" ~ "No Guides",
                                  Group == "With_Guides" ~ "Guides",
                                  TRUE ~ Group),
                        levels = c("No Guides", "Guides")))


p_s9 = ggplot(s9_long, aes(x = Variable, y = Value, color = Group)) +
  geom_boxplot() +
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.7, dodge.width = 0.8), alpha = 0.001, size = .005) +  # Jittered points
  labs(x = "Variable", y = "Value", fill = "Condition") +
  theme_minimal() +
  stat_compare_means(aes(group = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE)+ 
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "Nano 9", color = "Guides") +  # Change the Y-axis label here
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 22)) + # Increase X-axis label size
  coord_cartesian(clip = "off")
p_s9

## S10
s10_long <- s10 %>%
  select(!C_guides) %>%
  rename("C" = C) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
  mutate(Group = ifelse(stringr::str_detect(Variable, "_guides"), "With_Guides", "Without_Guides"),
         Variable = stringr::str_replace(Variable, "_guides", "")) %>%
  mutate(Group = factor(case_when(Group == "Without_Guides" ~ "No Guides",
                                  Group == "With_Guides" ~ "Guides",
                                  TRUE ~ Group),
                        levels = c("No Guides", "Guides")))


p_s10 = ggplot(s10_long, aes(x = Variable, y = Value, color = Group)) +
  geom_boxplot() +
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.7, dodge.width = 0.8), alpha = 0.001, size = .005) +  # Jittered points
  labs(x = "Variable", y = "Value", fill = "Condition") +
  theme_minimal() +
  stat_compare_means(aes(group = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE)+ 
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "Nano 10", color = "Guides") +  # Change the Y-axis label here
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 22)) + # Increase X-axis label size
  coord_cartesian(clip = "off")
p_s10

hvr_nano = 
  p_s2 / p_s3 / p_s4 / p_s5 / p_s6 / p_s7 / p_s8 / p_s9 / p_s10
hvr_nano

mean(s2$C)
sd(s2$C)
mean(s2$V2_guides)
sd(s2$V2_guides)
mean(s2$V2)
sd(s2$V2)
mean(s2$V4_guides)
sd(s2$V4_guides)
mean(s2$V4)
sd(s2$V4)

mean(s2$V4_guides) - mean(s2$V4)

setwd(dir.figures)
#ggsave("hvr_nano_structure_s2_s10.pdf", hvr_nano, height = 500, width = 200, units = "mm")
ggsave("hvr_nano_structure_s2_only.pdf", p_s2, height = 150, width = 150, units = "mm")

setwd(dir.out)

################################################################################
# MicroStructure: MAST
################################################################################

## MAST
setwd(dir.out)
mast = as_tibble(readRDS("hvr_mast_random_pairs.rds"))
mast %>% print(Inf)

## MAST
mast_long <- mast %>%
  select(!C_guides) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
  mutate(Group = ifelse(stringr::str_detect(Variable, "_guides"), "With_Guides", "Without_Guides"),
         Variable = stringr::str_replace(Variable, "_guides", "")) %>%
  mutate(Group = factor(case_when(Group == "Without_Guides" ~ "No Guides",
                                  Group == "With_Guides" ~ "Guides",
                                  TRUE ~ Group),
                        levels = c("No Guides", "Guides")))


p_mast = ggplot(mast_long, aes(x = Variable, y = Value, color = Group)) +
  geom_boxplot() +
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.7, dodge.width = 0.8), alpha = 0.1, size = .05) +  # Jittered points
  labs(x = "Variable", y = "Value", fill = "Condition") +
  theme_minimal() +
  stat_compare_means(aes(group = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE)+ 
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "MAST", color = "Guides") +  # Change the Y-axis label here
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 10)) + # Increase X-axis label size
  coord_cartesian(clip = "off")
p_mast

mean(mast$C)
sd(mast$C)

mean(mast$V4_guides)
sd(mast$V4_guides)

mean(mast$V4)
sd(mast$V4)

mean(mast$V4_guides) - mean(mast$V4)


################################################################################
# MicroStructure: MAST Iterative
################################################################################
# TODO Add vizualization for iterative mast here.
setwd(dir.out)
micro_totals = as_tibble(data.frame("C" =  micro_total(readRDS("hvr_micro_random_pairs_vfull.rds")), 
                                    "C_guides" =  micro_total(readRDS("hvr_micro_random_pairs_vfull_guides.rds")),
                                    "V1" =  micro_total(readRDS("hvr_micro_random_pairs_v1.rds")), 
                                    "V2" =  micro_total(readRDS("hvr_micro_random_pairs_v2.rds")), 
                                    "V3" =  micro_total(readRDS("hvr_micro_random_pairs_v3.rds")), 
                                    "V4" =  micro_total(readRDS("hvr_micro_random_pairs_v4.rds")), 
                                    "V5" =  micro_total(readRDS("hvr_micro_random_pairs_v5.rds")), 
                                    "V6" =  micro_total(readRDS("hvr_micro_random_pairs_v6.rds")), 
                                    "V7" =  micro_total(readRDS("hvr_micro_random_pairs_v7.rds")), 
                                    "V8" =  micro_total(readRDS("hvr_micro_random_pairs_v8.rds")), 
                                    "V9" =  micro_total(readRDS("hvr_micro_random_pairs_v9.rds")), 
                                    "V1_guides" =  micro_total(readRDS("hvr_micro_random_pairs_v1_guides.rds")), 
                                    "V2_guides" =  micro_total(readRDS("hvr_micro_random_pairs_v2_guides.rds")), 
                                    "V3_guides" =  micro_total(readRDS("hvr_micro_random_pairs_v3_guides.rds")), 
                                    "V4_guides" =  micro_total(readRDS("hvr_micro_random_pairs_v4_guides.rds")), 
                                    "V5_guides" =  micro_total(readRDS("hvr_micro_random_pairs_v5_guides.rds")), 
                                    "V6_guides" =  micro_total(readRDS("hvr_micro_random_pairs_v6_guides.rds")), 
                                    "V7_guides" =  micro_total(readRDS("hvr_micro_random_pairs_v7_guides.rds")), 
                                    "V8_guides" =  micro_total(readRDS("hvr_micro_random_pairs_v8_guides.rds")), 
                                    "V9_guides" =  micro_total(readRDS("hvr_micro_random_pairs_v9_guides.rds"))))
micro_totals = micro_totals/TREE_SIZE

micro_mean = as_tibble(data.frame("C" =  micro_mean_size(readRDS("hvr_micro_random_pairs_vfull.rds")), 
                                  "C_guides" =  micro_mean_size(readRDS("hvr_micro_random_pairs_vfull_guides.rds")),
                                  "V1" =  micro_mean_size(readRDS("hvr_micro_random_pairs_v1.rds")), 
                                  "V2" =  micro_mean_size(readRDS("hvr_micro_random_pairs_v2.rds")), 
                                  "V3" =  micro_mean_size(readRDS("hvr_micro_random_pairs_v3.rds")), 
                                  "V4" =  micro_mean_size(readRDS("hvr_micro_random_pairs_v4.rds")), 
                                  "V5" =  micro_mean_size(readRDS("hvr_micro_random_pairs_v5.rds")), 
                                  "V6" =  micro_mean_size(readRDS("hvr_micro_random_pairs_v6.rds")), 
                                  "V7" =  micro_mean_size(readRDS("hvr_micro_random_pairs_v7.rds")), 
                                  "V8" =  micro_mean_size(readRDS("hvr_micro_random_pairs_v8.rds")), 
                                  "V9" =  micro_mean_size(readRDS("hvr_micro_random_pairs_v9.rds")), 
                                  "V1_guides" =  micro_mean_size(readRDS("hvr_micro_random_pairs_v1_guides.rds")), 
                                  "V2_guides" =  micro_mean_size(readRDS("hvr_micro_random_pairs_v2_guides.rds")), 
                                  "V3_guides" =  micro_mean_size(readRDS("hvr_micro_random_pairs_v3_guides.rds")), 
                                  "V4_guides" =  micro_mean_size(readRDS("hvr_micro_random_pairs_v4_guides.rds")), 
                                  "V5_guides" =  micro_mean_size(readRDS("hvr_micro_random_pairs_v5_guides.rds")), 
                                  "V6_guides" =  micro_mean_size(readRDS("hvr_micro_random_pairs_v6_guides.rds")), 
                                  "V7_guides" =  micro_mean_size(readRDS("hvr_micro_random_pairs_v7_guides.rds")), 
                                  "V8_guides" =  micro_mean_size(readRDS("hvr_micro_random_pairs_v8_guides.rds")), 
                                  "V9_guides" =  micro_mean_size(readRDS("hvr_micro_random_pairs_v9_guides.rds"))))

## ITERATIVE MAST TOTALS
micro_total_long <- micro_totals %>%
  select(!C_guides) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
  mutate(Group = ifelse(stringr::str_detect(Variable, "_guides"), "With_Guides", "Without_Guides"),
         Variable = stringr::str_replace(Variable, "_guides", "")) %>%
  mutate(Group = factor(case_when(Group == "Without_Guides" ~ "No Guides",
                                  Group == "With_Guides" ~ "Guides",
                                  TRUE ~ Group),
                        levels = c("No Guides", "Guides")))


p_total = ggplot(micro_total_long, aes(x = Variable, y = Value, color = Group)) +
  geom_boxplot() +
  labs(x = "Variable", y = "Value", fill = "Condition") +
  theme_minimal() +
  stat_compare_means(aes(group = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE)+ 
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "% Tips in Shared Microstructures", color = "Guides") +  # Change the Y-axis label here
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 10)) + # Increase X-axis label size
  coord_cartesian(clip = "off")
p_total

## ITERATIVE MAST MEAN CLUSTER SIZE
micro_mean_long <- micro_mean %>%
  select(!C_guides) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
  mutate(Group = ifelse(stringr::str_detect(Variable, "_guides"), "With_Guides", "Without_Guides"),
         Variable = stringr::str_replace(Variable, "_guides", "")) %>%
  mutate(Group = factor(case_when(Group == "Without_Guides" ~ "No Guides",
                                  Group == "With_Guides" ~ "Guides",
                                  TRUE ~ Group),
                        levels = c("No Guides", "Guides")))


p_mean = ggplot(micro_mean_long, aes(x = Variable, y = Value, color = Group)) +
  geom_boxplot() +
  labs(x = "Variable", y = "Value", fill = "Condition") +
  theme_minimal() +
  stat_compare_means(aes(group = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE)+ 
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "Log(Mean Microstructure Size)", color = "Guides") +  # Change the Y-axis label here
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 10)) + # Increase X-axis label size
  coord_cartesian(clip = "off")
p_mean
hvr_micro = 
  p_mast / p_total / p_mean
setwd(dir.figures)
ggsave(file = "hvr_microstructure.pdf", hvr_micro, width = 150, height = 250, units = "mm")
ggsave(file = "hvr_microstructure_mast.pdf", p_mast, width = 150, height = 150, units = "mm")

################################################################################
# Clean up
################################################################################
rm(list = c(ls(pattern = "^box_"), ls(pattern = "_long$"), ls(pattern = "^p_"), "final_plot", ls(pattern = "^heat_")))
rm(list = c(ls(pattern = "^col_")))
rm(list = c(ls(pattern = "s[0-9]"), "mast", "pearson_cor_hvr", "title_accuracy", "title_heatmaps", "final_plot_B", "final_plot_A", ls(pattern = "hvr_")))
rm(list = c("layout", "hvrs"))
rm(list = c("micro_mean", "micro_totals"))
rm(list = "labels")

################################################################################
# SIM 2: 
# Goal plot global measures of diversity for each length vs the other.
################################################################################

################################################################################
# Global Distances: 100 Random pairwise comparisons 
# Determine if guide sequence improve topological accuracy at each HVR
################################################################################

# Read in the distance metrics
setwd(dir.sim2)
MATCHES = c("VFull_control", "VFull_guides")

## 1 RF
bp_rf = as_tibble(readRDS("bp_rf_random_pairs.rds")) %>% rename_with(~ str_remove(., "_rf")) %>%
  select(!all_of(MATCHES))
labels = colnames(bp_rf)
labels

## 2 JRF
bp_jrf = as_tibble(readRDS("bp_jrf_random_pairs.rds")) %>% rename_with(~ str_remove(., "_jrf")) %>% 
  select(!all_of(MATCHES))
if(setequal(colnames(bp_jrf), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 3 Nye Distance
bp_nye = as_tibble(readRDS("bp_nye_random_pairs.rds")) %>% rename_with(~ str_remove(., "_nye")) %>% 
  select(!all_of(MATCHES))
if(setequal(colnames(bp_nye), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 4 MSD
bp_msd = as_tibble(readRDS("bp_msd_random_pairs.rds")) %>% rename_with(~ str_remove(., "_msd")) %>% 
  select(!all_of(MATCHES))
if(setequal(colnames(bp_msd), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 5 PID
bp_pid = as_tibble(readRDS("bp_pid_random_pairs.rds"))  %>% rename_with(~ str_remove(., "_pid"))%>% 
  select(!all_of(MATCHES))
if(setequal(colnames(bp_pid), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 6 MCD
bp_mcd = as_tibble(readRDS("bp_mcd_random_pairs.rds"))  %>% rename_with(~ str_remove(., "_mcd")) %>%
  select(!all_of(MATCHES))
if(setequal(colnames(bp_mcd), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 7 Path
bp_path = as_tibble(readRDS("bp_path_random_pairs.rds"))%>% rename_with(~ str_remove(., "_path")) %>% 
  select(!all_of(MATCHES))
if(setequal(colnames(bp_path), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 8 Quartet Distance
bp_quartet = as_tibble(readRDS("bp_quartet_normalized_random_pairs.rds")) %>% rename_with(~ str_remove(., "_quartet")) %>% 
  select(!all_of(MATCHES))
if(setequal(colnames(bp_quartet), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

# Check if each one of these is normally distributed
# Test for normality for each response variable
test_normality = function(dist){
  p = apply(dist, 2, function(x){
  return(shapiro.test(x)$p.value)
  })
}

# Test for normality of each of the 
normality = data.frame("rf" = test_normality(bp_rf[labels]),
                       "jrf" = test_normality(bp_jrf[labels]),
                       "nye" = test_normality(bp_nye[labels]),
                       "msd"= test_normality(bp_msd[labels]),
                       "pid" = test_normality((bp_pid[labels])),
                       "mcd" = test_normality((bp_mcd[labels])),
                       "path" = test_normality((bp_path[labels])),
                       "quar" = test_normality((bp_quartet[labels])))
normality_q = data.frame(apply(normality, 2, function(p) p.adjust(p, method = "bonferroni")))
normality_q$names = rownames(normality_q)
Q = 0.05
normality_q %>% filter(pid < Q) %>% pull(names)

# Mostly normal per level
normality_q %>% filter(rf < Q) %>% nrow()
normality_q %>% filter(jrf < Q) %>% nrow()
normality_q %>% filter(nye < Q) %>% nrow()

# Lots of non-normal distributions per level
normality_q %>% filter(pid < Q) %>% nrow() #23
normality_q %>% filter(msd < Q) %>% nrow() # 30
normality_q %>% filter(mcd < Q) %>% nrow() # 25
normality_q %>% filter(path < Q) %>% nrow() #26
normality_q %>% filter(quar < Q) %>% nrow() #71

## 1 RF
# Reshape data to long format
bp_rf_long <- 
  bp_rf %>%
  pivot_longer(cols = everything(), names_to = "Simulation", values_to = "Dist") %>% 
  mutate(Group = ifelse(stringr::str_detect(Simulation, "_guides"), "With_Guides", "Without_Guides"),
         VR = case_when(stringr::str_detect(Simulation, "V1") ~ "V1",
                        stringr::str_detect(Simulation, "V2") ~ "V2",
                        stringr::str_detect(Simulation, "V3") ~ "V3",
                        stringr::str_detect(Simulation, "V4") ~ "V4",
                        stringr::str_detect(Simulation, "V5") ~ "V5",
                        stringr::str_detect(Simulation, "V6") ~ "V6",
                        stringr::str_detect(Simulation, "V7") ~ "V7",
                        stringr::str_detect(Simulation, "V8") ~ "V8",
                        TRUE ~ NA_character_,),
         Length = stringr::str_replace(Simulation, "_guides", ""),
         Length = stringr::str_replace(Length, "^V[1-9]_", ""),
         Length = stringr::str_replace(Length, "^L", ""),
         Length = as.numeric(Length)) %>%
  mutate(Group = factor(case_when(Group == "Without_Guides" ~ "No Guides",
                                  Group == "With_Guides" ~ "Guides",
                                  TRUE ~ Group),
                        levels = c("No Guides", "Guides")))
bp_rf_long$VR <- factor(bp_rf_long$VR, levels = c("V4", "V1", "V2", "V3", "V5", "V6", "V7", "V8"))
bp_rf_long

# See if LM violates assumptions
rf_lm = lm(Dist ~ Length*VR*Group, data = bp_rf_long )
summary(rf_lm)


# RF violates several assumptions used for linear model
# Residual errors have a mean value of zero
plot(rf_lm, 1) # some structure present in residuals
car::durbinWatsonTest(rf_lm) #p = 0 accept the null hypothesis that there is significant autocorrelation in residuals
plot(rf_lm, 1) # the fitted line should be at zero, which it is not
car::ncvTest(rf_lm) # Heteroscedasticity is present.

# Will check that the beta logit function, which is an easy model to intepret is reasonable for AIC / BIC of others
rf_beta_logit = betareg::betareg(Dist ~ Length + VR + Group + VR*Group + VR*Length + Group*Length, data = bp_rf_long, link = "logit")
summary(rf_beta_logit)
rf_beta_probit = betareg::betareg(Dist ~ Length + VR + Group + VR*Group + VR*Length + Group*Length, data = bp_rf_long, link = "probit")
rf_beta_cloglog = betareg::betareg(Dist ~ Length + VR + Group + VR*Group + VR*Length + Group*Length, data = bp_rf_long, link = "cloglog")
rf_beta_cauchit = betareg::betareg(Dist ~ Length + VR + Group + VR*Group + VR*Length + Group*Length, data = bp_rf_long, link = "cauchit")
rf_beta_log = betareg::betareg(Dist ~ Length + VR + Group + VR*Group + VR*Length + Group*Length, data = bp_rf_long, link = "log")
rf_beta_loglog = betareg::betareg(Dist ~ Length + VR + Group + VR*Group + VR*Length + Group*Length, data = bp_rf_long, link = "loglog")
rf_gamma_log <- glm(Dist ~ Length + VR + Group + VR*Group + VR*Length + Group*Length, family = Gamma(link = "log"), data = bp_rf_long)

aics = AIC(rf_beta_logit, rf_beta_probit, rf_beta_cloglog, rf_beta_cauchit, rf_beta_loglog, rf_lm, rf_gamma_log, k = 2)
aics$model = rownames(aics)

bics = BIC(rf_beta_logit, rf_beta_probit, rf_beta_cloglog, rf_beta_cauchit, rf_beta_loglog, rf_lm, rf_gamma_log)
bics$model = rownames(bics)
as_tibble(aics) %>% arrange(AIC)
as_tibble(bics) %>% arrange(BIC)

# Final Model selected 
summary(rf_beta_logit)

# No plot out the results to look at 
bp_rf_long = 
  bp_rf_long %>%
  mutate(Length = as.character(Length)) %>%
  mutate(Length = stringr::str_replace(string = Length, pattern = "^", replacement = "L"))
bp_rf_long$Length = factor(bp_rf_long$Length, levels = c("L50", "L100", "L200", "L300", "L400", "L500"))


V1 = ggplot(bp_rf_long  %>% filter(VR == "V1")) +
  geom_boxplot(aes(x = Length, y = Dist,  color = Group), alpha = 0.5, outlier.shape = NA) + # Box plot without outlier dots
  labs(title = "VR1", 
       x = "Length", 
       y = "RF") +
  ylim(c(0, 1)) + 
  stat_compare_means(aes(x = Length, y = Dist,  color = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE) + 
  theme_minimal() +
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                   labels = c("-", "+")) +
  labs(x = "", y = "RF", color = "Guides") +  # Change the Y-axis label here
  
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        title = element_text(size = 20, face = "bold")) + # Increase X-axis label size
  coord_cartesian(clip = "off")
V1

V2 = ggplot(bp_rf_long %>% filter(VR == "V2")) +
  geom_boxplot(aes(x = Length, y = Dist,  color = Group), alpha = 0.5, outlier.shape = NA) + # Box plot without outlier dots
  labs(title = "VR2", 
       x = "Length", 
       y = "RF") +
  ylim(c(0, 1)) + 
  stat_compare_means(aes(x = Length, y = Dist,  color = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE) + 
  theme_minimal() +
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "RF", color = "Guides") +  # Change the Y-axis label here
  
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        title = element_text(size = 20, face = "bold")) + # Increase X-axis label size
  coord_cartesian(clip = "off")
V2

V3 = ggplot(bp_rf_long %>% filter(VR == "V3")) +
  geom_boxplot(aes(x = Length, y = Dist,  color = Group), alpha = 0.5, outlier.shape = NA) + # Box plot without outlier dots
  labs(title = "VR3", 
       x = "Length", 
       y = "RF") +
  ylim(c(0, 1)) + 
  stat_compare_means(aes(x = Length, y = Dist,  color = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE) + 
  theme_minimal() +
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "RF", color = "Guides") +  # Change the Y-axis label here
  
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        title = element_text(size = 20, face = "bold")) + # Increase X-axis label size
  coord_cartesian(clip = "off")
V3

V4 = ggplot(bp_rf_long %>% filter(VR == "V4")) +
  geom_boxplot(aes(x = Length, y = Dist,  color = Group), alpha = 0.5, outlier.shape = NA) + # Box plot without outlier dots
  labs(title = "VR4", 
       x = "Length", 
       y = "RF") +
  ylim(c(0, 1)) + 
  stat_compare_means(aes(x = Length, y = Dist,  color = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE) + 
  theme_minimal() +
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "RF", color = "Guides") +  # Change the Y-axis label here
  
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        title = element_text(size = 20, face = "bold")) + # Increase X-axis label size
  coord_cartesian(clip = "off")
V4


V5 = ggplot(bp_rf_long %>% filter(VR == "V5")) +
  geom_boxplot(aes(x = Length, y = Dist,  color = Group), alpha = 0.5, outlier.shape = NA) + # Box plot without outlier dots
  labs(title = "VR5", 
       x = "Length", 
       y = "RF") +
  ylim(c(0, 1)) + 
  stat_compare_means(aes(x = Length, y = Dist,  color = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE) + 
  theme_minimal() +
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "RF", color = "Guides") +  # Change the Y-axis label here
  
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        title = element_text(size = 20, face = "bold")) + # Increase X-axis label size
  coord_cartesian(clip = "off")
V5

V6 = ggplot(bp_rf_long %>% filter(VR == "V6")) +
  geom_boxplot(aes(x = Length, y = Dist,  color = Group), alpha = 0.5, outlier.shape = NA) + # Box plot without outlier dots
  labs(title = "VR6", 
       x = "Length", 
       y = "RF") +
  ylim(c(0, 1)) + 
  scale_x_discrete(limits = c("L50", "L100", "L200", "L300", "L400", "L500")) +
  stat_compare_means(aes(x = Length, y = Dist,  color = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE) + 
  theme_minimal() +
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "RF", color = "Guides") +  # Change the Y-axis label here
  
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        title = element_text(size = 20, face = "bold")) + # Increase X-axis label size
  coord_cartesian(clip = "off")
V6

V7 = ggplot(bp_rf_long %>% filter(VR == "V7")) +
  geom_boxplot(aes(x = Length, y = Dist,  color = Group), alpha = 0.5, outlier.shape = NA) + # Box plot without outlier dots
  labs(title = "VR7", 
       x = "Length", 
       y = "RF") +
  ylim(c(0, 1)) + 
  scale_x_discrete(limits = c("L50", "L100", "L200", "L300", "L400", "L500")) +
  stat_compare_means(aes(x = Length, y = Dist,  color = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE) + 
  theme_minimal() +
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "RF", color = "Guides") +  # Change the Y-axis label here
  
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        title = element_text(size = 20, face = "bold")) + # Increase X-axis label size
  coord_cartesian(clip = "off")
V7

V8 = ggplot(bp_rf_long %>% filter(VR == "V8")) +
  geom_boxplot(aes(x = Length, y = Dist,  color = Group), alpha = 0.5, outlier.shape = NA) + # Box plot without outlier dots
  labs(title = "VR8", 
       x = "Length", 
       y = "RF") +
  ylim(c(0, 1)) + 
  scale_x_discrete(limits = c("L50", "L100", "L200", "L300", "L400", "L500")) +
  stat_compare_means(aes(x = Length, y = Dist,  color = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE) + 
  theme_minimal() +
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "RF", color = "Guides") +  # Change the Y-axis label here
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        title = element_text(size = 20, face = "bold")) + # Increase X-axis label size 
  coord_cartesian(clip = "off")
V8

# Extract the legend to use
legend <- get_legend(V1 + theme(legend.position = "bottom"))

#Remove the legends from all the others
V1 = V1 + theme(legend.position = "none")
V2 = V2 + theme(legend.position = "none")
V3 = V3 + theme(legend.position = "none")
V4 = V4 + theme(legend.position = "none")
V5 = V5 + theme(legend.position = "none")
V6 = V6 + theme(legend.position = "none")
V7 = V7 + theme(legend.position = "none")
V8 = V8 + theme(legend.position = "none")
top_row <- (V1 + V2 + V3 + V4)
bottom_row <- (V5 + V6 + V7 + V8) 

# Final layout with the legend placed at the bottom
final_plot <- (top_row / bottom_row) / legend + plot_layout(heights = c(1, 1, 0.1))
final_plot 

setwd(dir.figures)
ggsave(filename = "length_vr_rf_plots.pdf", height = 500, width = 300, units = "mm")

# Now test if each one fof the global values is correlated
bp_rf = bp_rf %>% mutate(simulation = paste0("sim", row_number()))
bp_rf_long = 
  as_tibble(melt(bp_rf)) %>%
  unite(col = "ID", simulation, variable, sep = "_") %>%
  rename("RF" = value) 
bp_jrf = bp_jrf %>% mutate(simulation = paste0("sim", row_number()))
bp_jrf_long = 
  as_tibble(melt(bp_jrf)) %>%
  unite(col = "ID", simulation, variable, sep = "_") %>%
  rename("JRF" = value) 
bp_nye = bp_nye %>% mutate(simulation = paste0("sim", row_number()))
bp_nye_long = 
  as_tibble(melt(bp_nye)) %>%
  unite(col = "ID", simulation, variable, sep = "_") %>%
  rename("NYE" = value) 
bp_msd = bp_msd %>% mutate(simulation = paste0("sim", row_number()))
bp_msd_long = 
  as_tibble(melt(bp_msd)) %>%
  unite(col = "ID", simulation, variable, sep = "_") %>%
  rename("MSD" = value) 
bp_pid = bp_pid %>% mutate(simulation = paste0("sim", row_number()))
bp_pid_long = 
  as_tibble(melt(bp_pid)) %>%
  unite(col = "ID", simulation, variable, sep = "_") %>%
  rename("PID" = value) 
bp_mcd = bp_mcd %>% mutate(simulation = paste0("sim", row_number()))
bp_mcd_long = 
  as_tibble(melt(bp_mcd)) %>%
  unite(col = "ID", simulation, variable, sep = "_") %>%
  rename("MCD" = value) 
bp_path = bp_path %>% mutate(simulation = paste0("sim", row_number()))
bp_path_long = 
  as_tibble(melt(bp_path)) %>%
  unite(col = "ID", simulation, variable, sep = "_") %>%
  rename("PATH" = value) 
bp_quartet = bp_quartet %>% mutate(simulation = paste0("sim", row_number()))
bp_quartet_long = 
  as_tibble(melt(bp_quartet)) %>%
  unite(col = "ID", simulation, variable, sep = "_") %>%
  rename("QUAR" = value) 

g_dists = 
  bp_rf_long %>%
  left_join(bp_jrf_long, join_by("ID")) %>%
  left_join(bp_nye_long, join_by("ID")) %>%
  left_join(bp_msd_long, join_by("ID")) %>%
  left_join(bp_pid_long, join_by("ID")) %>%
  left_join(bp_mcd_long, join_by("ID")) %>%
  left_join(bp_path_long, join_by("ID")) %>%
  left_join(bp_quartet_long, join_by("ID"))
g_dists  
cor.test(g_dists$RF, g_dists$JRF)
cor.test(g_dists$RF, g_dists$NYE)
cor.test(g_dists$RF, g_dists$MSD)
cor.test(g_dists$RF, g_dists$PID)
cor.test(g_dists$RF, g_dists$MCD)
cor.test(g_dists$RF, g_dists$PATH)
cor.test(g_dists$RF, g_dists$QUAR)

# Check if average global dists correlated
mean_dist = data.frame("RF" = apply(bp_rf, 2, mean)[labels],
                       "JRF" = apply(bp_jrf, 2, mean)[labels],
                       "NYE" = apply(bp_nye, 2, mean)[labels],
                       "MSD" = apply(bp_msd, 2, mean)[labels],
                       "PID" = apply(bp_pid, 2, mean)[labels],
                       "MCD" = apply(bp_mcd, 2, mean)[labels],
                       "PATH" = apply(bp_path, 2, mean)[labels],
                       "QUAR" = apply(bp_quartet, 2, mean)[labels])
mean_dist
cor.test(mean_dist$RF, mean_dist$JRF)
cor.test(mean_dist$RF, mean_dist$NYE)
cor.test(mean_dist$RF, mean_dist$MSD)
cor.test(mean_dist$RF, mean_dist$PID)
cor.test(mean_dist$RF, mean_dist$MCD)
cor.test(mean_dist$RF, mean_dist$PATH)
cor.test(mean_dist$RF, mean_dist$QUAR)


################################################################################
# Nano Distances: 100 Random pairwise comparisons 
################################################################################
setwd(dir.sim2)
MATCHES = c("VFull_control_subtree", "VFull_guides_subtree")
s2 = as_tibble(readRDS("nanostructure_bp_2_random_pairs.rds")) %>% select(!all_of(MATCHES))
s2
s3 = as_tibble(readRDS("nanostructure_bp_3_random_pairs.rds"))%>% select(!all_of(MATCHES))
s3
s4 = as_tibble(readRDS("nanostructure_bp_4_random_pairs.rds"))%>% select(!all_of(MATCHES))
s4
s5 = as_tibble(readRDS("nanostructure_bp_5_random_pairs.rds"))%>% select(!all_of(MATCHES))
s5
s6 = as_tibble(readRDS("nanostructure_bp_6_random_pairs.rds"))%>% select(!all_of(MATCHES))
s6
s7 = as_tibble(readRDS("nanostructure_bp_7_random_pairs.rds"))%>% select(!all_of(MATCHES))
s7
s8 = as_tibble(readRDS("nanostructure_bp_8_random_pairs.rds"))%>% select(!all_of(MATCHES))
s8
s9 = as_tibble(readRDS("nanostructure_bp_9_random_pairs.rds"))%>% select(!all_of(MATCHES))
s9
s10 = as_tibble(readRDS("nanostructure_bp_10_random_pairs.rds"))%>% select(!all_of(MATCHES))
s10

bp_s2_long <- 
  s2 %>%
  pivot_longer(cols = everything(), names_to = "Simulation", values_to = "Dist") %>% 
  mutate(Group = ifelse(stringr::str_detect(Simulation, "_guides"), "With_Guides", "Without_Guides"),
         VR = case_when(stringr::str_detect(Simulation, "V1") ~ "V1",
                        stringr::str_detect(Simulation, "V2") ~ "V2",
                        stringr::str_detect(Simulation, "V3") ~ "V3",
                        stringr::str_detect(Simulation, "V4") ~ "V4",
                        stringr::str_detect(Simulation, "V5") ~ "V5",
                        stringr::str_detect(Simulation, "V6") ~ "V6",
                        stringr::str_detect(Simulation, "V7") ~ "V7",
                        stringr::str_detect(Simulation, "V8") ~ "V8",
                        TRUE ~ NA_character_,),
         Length = stringr::str_replace(Simulation, "_subtree$", ""),
         Length = stringr::str_replace(Length, "_guides", ""),
         Length = stringr::str_replace(Length, "^V[1-9]_", ""),
         Length = stringr::str_replace(Length, "^L", ""),
         Length = as.numeric(Length)) %>%
  mutate(Group = factor(case_when(Group == "Without_Guides" ~ "No Guides",
                                  Group == "With_Guides" ~ "Guides",
                                  TRUE ~ Group),
                        levels = c("No Guides", "Guides")))
bp_s2_long

bp_s2_long$VR <- factor(bp_s2_long$VR, levels = c("V4", "V1", "V2", "V3", "V5", "V6", "V7", "V8"))
bp_s2_long

# See if LM violates assumptions
rf_lm = lm(Dist ~ Length*VR*Group, data = bp_s2_long )
summary(rf_lm)


# RF violates several assumptions used for linear model
# Residual errors have a mean value of zero
plot(rf_lm, 1) # some structure present in residuals
car::durbinWatsonTest(rf_lm) #p = 0 accept the null hypothesis that there is significant autocorrelation in residuals
plot(rf_lm, 1) # the fitted line should be at zero, which it is not
car::ncvTest(rf_lm) # Heteroscedasticity is present.
shapiro.test(as.vector(bp_s2_long[sample(x = seq(from = 1, to = 5000, by = 1), replace = FALSE, size = 5000), "Dist"]$Dist))

# Will check that the beta logit function, which is an easy model to intepret is reasonable for AIC / BIC of others
rf_beta_logit = betareg::betareg(Dist ~ Length + VR + Group + VR*Group + VR*Length + Group*Length, data = bp_s2_long, link = "logit")
rf_beta_probit = betareg::betareg(Dist ~ Length + VR + Group + VR*Group + VR*Length + Group*Length, data = bp_s2_long, link = "probit")

rf_beta_cloglog = betareg::betareg(Dist ~ Length + VR + Group + VR*Group + VR*Length + Group*Length, data = bp_s2_long, link = "cloglog")
rf_beta_cauchit = betareg::betareg(Dist ~ Length + VR + Group + VR*Group + VR*Length + Group*Length, data = bp_s2_long, link = "cauchit")
rf_beta_loglog = betareg::betareg(Dist ~ Length + VR + Group + VR*Group + VR*Length + Group*Length, data = bp_s2_long, link = "loglog")
rf_gamma_log <- glm(Dist ~ Length + VR + Group + VR*Group + VR*Length + Group*Length, family = Gamma(link = "log"), data = bp_s2_long)

aics = AIC(rf_beta_logit, rf_beta_probit, rf_beta_cloglog, rf_beta_cauchit, rf_beta_loglog, rf_lm, rf_gamma_log, k = 2)
aics$model = rownames(aics)

bics = BIC(rf_beta_logit, rf_beta_probit, rf_beta_cloglog, rf_beta_cauchit, rf_beta_loglog, rf_lm, rf_gamma_log)
bics$model = rownames(bics)
as_tibble(aics) %>% arrange(AIC)
as_tibble(bics) %>% arrange(BIC)

# Final Model selected 
summary(rf_beta_logit)

# Graph out differences for figure for s2
bp_s2_long = 
  bp_s2_long %>%
  mutate(Length = as.character(Length)) %>%
  mutate(Length = stringr::str_replace(string = Length, pattern = "^", replacement = "L"))
bp_s2_long$Length = factor(bp_s2_long$Length, levels = c("L50", "L100", "L200", "L300", "L400", "L500"))
bp_s2_long

V1 = ggplot(bp_s2_long  %>% filter(VR == "V1")) +
  geom_boxplot(aes(x = Length, y = Dist,  color = Group), alpha = 0.5, outlier.shape = NA) + # Box plot without outlier dots
  labs(title = "VR1", 
       x = "Length", 
       y = "NSM (S=2)") +
  ylim(c(0, 1)) + 
  stat_compare_means(aes(x = Length, y = Dist,  color = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE) + 
  theme_minimal() +
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "NSM (S=2)", color = "Guides") +  # Change the Y-axis label here
  
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        title = element_text(size = 20, face = "bold")) + # Increase X-axis label size
  coord_cartesian(clip = "off")
V1

V2 = ggplot(bp_s2_long %>% filter(VR == "V2")) +
  geom_boxplot(aes(x = Length, y = Dist,  color = Group), alpha = 0.5, outlier.shape = NA) + # Box plot without outlier dots
  labs(title = "VR2", 
       x = "Length", 
       y = "NSM (S=2)") +
  ylim(c(0, 1)) + 
  stat_compare_means(aes(x = Length, y = Dist,  color = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE) + 
  theme_minimal() +
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "NSM (S=2)", color = "Guides") +  # Change the Y-axis label here
  
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        title = element_text(size = 20, face = "bold")) + # Increase X-axis label size
  coord_cartesian(clip = "off")
V2

V3 = ggplot(bp_s2_long %>% filter(VR == "V3")) +
  geom_boxplot(aes(x = Length, y = Dist,  color = Group), alpha = 0.5, outlier.shape = NA) + # Box plot without outlier dots
  labs(title = "VR3", 
       x = "Length", 
       y = "NSM (S=2)") +
  ylim(c(0, 1)) + 
  stat_compare_means(aes(x = Length, y = Dist,  color = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE) + 
  theme_minimal() +
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "NSM (S=2)", color = "Guides") +  # Change the Y-axis label here
  
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        title = element_text(size = 20, face = "bold")) + # Increase X-axis label size
  coord_cartesian(clip = "off")
V3

V4 = ggplot(bp_s2_long %>% filter(VR == "V4")) +
  geom_boxplot(aes(x = Length, y = Dist,  color = Group), alpha = 0.5, outlier.shape = NA) + # Box plot without outlier dots
  labs(title = "VR4", 
       x = "Length", 
       y = "NSM (S=2)") +
  ylim(c(0, 1)) + 
  stat_compare_means(aes(x = Length, y = Dist,  color = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE) + 
  theme_minimal() +
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "NSM (S=2)", color = "Guides") +  # Change the Y-axis label here
  
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        title = element_text(size = 20, face = "bold")) + # Increase X-axis label size
  coord_cartesian(clip = "off")
V4


V5 = ggplot(bp_s2_long %>% filter(VR == "V5")) +
  geom_boxplot(aes(x = Length, y = Dist,  color = Group), alpha = 0.5, outlier.shape = NA) + # Box plot without outlier dots
  labs(title = "VR5", 
       x = "Length", 
       y = "NSM (S=2)") +
  ylim(c(0, 1)) + 
  stat_compare_means(aes(x = Length, y = Dist,  color = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE) + 
  theme_minimal() +
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "NSM (S=2)", color = "Guides") +  # Change the Y-axis label here
  
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        title = element_text(size = 20, face = "bold")) + # Increase X-axis label size
  coord_cartesian(clip = "off")
V5

V6 = ggplot(bp_s2_long %>% filter(VR == "V6")) +
  geom_boxplot(aes(x = Length, y = Dist,  color = Group), alpha = 0.5, outlier.shape = NA) + # Box plot without outlier dots
  labs(title = "VR6", 
       x = "Length", 
       y = "NSM (S=2)") +
  ylim(c(0, 1)) + 
  scale_x_discrete(limits = c("L50", "L100", "L200", "L300", "L400", "L500")) +
  stat_compare_means(aes(x = Length, y = Dist,  color = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE) + 
  theme_minimal() +
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "NSM (S=2)", color = "Guides") +  # Change the Y-axis label here
  
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        title = element_text(size = 20, face = "bold")) + # Increase X-axis label size
  coord_cartesian(clip = "off")
V6

V7 = ggplot(bp_s2_long %>% filter(VR == "V7")) +
  geom_boxplot(aes(x = Length, y = Dist,  color = Group), alpha = 0.5, outlier.shape = NA) + # Box plot without outlier dots
  labs(title = "VR7", 
       x = "Length", 
       y = "NSM (S=2)") +
  ylim(c(0, 1)) + 
  scale_x_discrete(limits = c("L50", "L100", "L200", "L300", "L400", "L500")) +
  stat_compare_means(aes(x = Length, y = Dist,  color = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE) + 
  theme_minimal() +
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "NSM (S=2)", color = "Guides") +  # Change the Y-axis label here
  
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        title = element_text(size = 20, face = "bold")) + # Increase X-axis label size
  coord_cartesian(clip = "off")
V7

V8 = ggplot(bp_s2_long %>% filter(VR == "V8")) +
  geom_boxplot(aes(x = Length, y = Dist,  color = Group), alpha = 0.5, outlier.shape = NA) + # Box plot without outlier dots
  labs(title = "VR8", 
       x = "Length", 
       y = "NSM (S=2)") +
  ylim(c(0, 1)) + 
  scale_x_discrete(limits = c("L50", "L100", "L200", "L300", "L400", "L500")) +
  stat_compare_means(aes(x = Length, y = Dist,  color = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE) + 
  theme_minimal() +
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "NSM (S=2)", color = "Guides") +  # Change the Y-axis label here
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        title = element_text(size = 20, face = "bold")) + # Increase X-axis label size 
  coord_cartesian(clip = "off")
V8

# Extract the legend to use
legend <- get_legend(V1 + theme(legend.position = "bottom"))

#Remove the legends from all the others
V1 = V1 + theme(legend.position = "none")
V2 = V2 + theme(legend.position = "none")
V3 = V3 + theme(legend.position = "none")
V4 = V4 + theme(legend.position = "none")
V5 = V5 + theme(legend.position = "none")
V6 = V6 + theme(legend.position = "none")
V7 = V7 + theme(legend.position = "none")
V8 = V8 + theme(legend.position = "none")
top_row <- (V1 + V2 + V3 + V4)
bottom_row <- (V5 + V6 + V7 + V8) 

# Final layout with the legend placed at the bottom
final_plot <- (top_row / bottom_row) / legend + plot_layout(heights = c(1, 1, 0.1))
final_plot 

setwd(dir.figures)
ggsave(filename = "length_vr_s2_plots.pdf", height = 500, width = 300, units = "mm")


# Now look if nanostructure of s2 correlates with nanastructure of size 3, 4, ..., 9, 10
s2 = s2 %>% mutate(simulation = paste0("sim", row_number()))
s2_long = 
  as_tibble(melt(s2)) %>%
  unite(col = "ID", simulation, variable, sep = "_") %>%
  rename("S2" = value) 
s2_long

s3 = s3 %>% mutate(simulation = paste0("sim", row_number()))
s3_long = 
  as_tibble(melt(s3)) %>%
  unite(col = "ID", simulation, variable, sep = "_") %>%
  rename("S3" = value) 
s3_long

s4 = s4 %>% mutate(simulation = paste0("sim", row_number()))
s4_long = 
  as_tibble(melt(s4)) %>%
  unite(col = "ID", simulation, variable, sep = "_") %>%
  rename("S4" = value) 
s4_long

s5 = s5 %>% mutate(simulation = paste0("sim", row_number()))
s5_long = 
  as_tibble(melt(s5)) %>%
  unite(col = "ID", simulation, variable, sep = "_") %>%
  rename("S5" = value) 
s5_long

s6 = s6 %>% mutate(simulation = paste0("sim", row_number()))
s6_long = 
  as_tibble(melt(s6)) %>%
  unite(col = "ID", simulation, variable, sep = "_") %>%
  rename("S6" = value) 
s6_long

s7 = s7 %>% mutate(simulation = paste0("sim", row_number()))
s7_long = 
  as_tibble(melt(s7)) %>%
  unite(col = "ID", simulation, variable, sep = "_") %>%
  rename("S7" = value) 
s7_long

s8 = s8 %>% mutate(simulation = paste0("sim", row_number()))
s8_long = 
  as_tibble(melt(s8)) %>%
  unite(col = "ID", simulation, variable, sep = "_") %>%
  rename("S8" = value) 
s8_long

s9 = s9 %>% mutate(simulation = paste0("sim", row_number()))
s9_long = 
  as_tibble(melt(s9)) %>%
  unite(col = "ID", simulation, variable, sep = "_") %>%
  rename("S9" = value) 
s9_long

s10 = s10 %>% mutate(simulation = paste0("sim", row_number()))
s10_long = 
  as_tibble(melt(s10)) %>%
  unite(col = "ID", simulation, variable, sep = "_") %>%
  rename("S10" = value) 
s10_long

s_dists = 
  s2_long %>%
  left_join(s3_long, join_by("ID")) %>%
  left_join(s4_long, join_by("ID")) %>%
  left_join(s5_long, join_by("ID")) %>%
  left_join(s6_long, join_by("ID")) %>%
  left_join(s7_long, join_by("ID")) %>%
  left_join(s8_long, join_by("ID")) %>%
  left_join(s9_long, join_by("ID")) %>%
  left_join(s10_long, join_by("ID"))

cor.test(s_dists$S2, s_dists$S3)
cor.test(s_dists$S2, s_dists$S4)
cor.test(s_dists$S2, s_dists$S5)
cor.test(s_dists$S2, s_dists$S6)
cor.test(s_dists$S2, s_dists$S7)
cor.test(s_dists$S2, s_dists$S8)
cor.test(s_dists$S2, s_dists$S9)
cor.test(s_dists$S2, s_dists$S10)
################################################################################
# MicroStructure: MAST
################################################################################
setwd(dir.sim2)
mast = as_tibble(readRDS("bp_mast.rds")) %>%
  select(!VFull_control) %>% 
  select(!VFull_guides)
mast

mast_long <- 
  mast %>%
  pivot_longer(cols = everything(), names_to = "Simulation", values_to = "Dist") %>% 
  mutate(Group = ifelse(stringr::str_detect(Simulation, "_guides"), "With_Guides", "Without_Guides"),
         VR = case_when(stringr::str_detect(Simulation, "V1") ~ "V1",
                        stringr::str_detect(Simulation, "V2") ~ "V2",
                        stringr::str_detect(Simulation, "V3") ~ "V3",
                        stringr::str_detect(Simulation, "V4") ~ "V4",
                        stringr::str_detect(Simulation, "V5") ~ "V5",
                        stringr::str_detect(Simulation, "V6") ~ "V6",
                        stringr::str_detect(Simulation, "V7") ~ "V7",
                        stringr::str_detect(Simulation, "V8") ~ "V8",
                        TRUE ~ NA_character_,),
         Length = stringr::str_replace(Simulation, "_subtree$", ""),
         Length = stringr::str_replace(Length, "_guides", ""),
         Length = stringr::str_replace(Length, "^V[1-9]_", ""),
         Length = stringr::str_replace(Length, "^L", ""),
         Length = as.numeric(Length)) %>%
  mutate(Group = factor(case_when(Group == "Without_Guides" ~ "No Guides",
                                  Group == "With_Guides" ~ "Guides",
                                  TRUE ~ Group),
                        levels = c("No Guides", "Guides")))
mast_long

# Graph out differences for figure for s2
mast_long = 
  mast_long %>%
  mutate(Length = as.character(Length)) %>%
  mutate(Length = stringr::str_replace(string = Length, pattern = "^", replacement = "L"))
mast_long$Length = factor(mast_long$Length, levels = c("L50", "L100", "L200", "L300", "L400", "L500"))
mast_long

V1 = ggplot(mast_long  %>% filter(VR == "V1")) +
  geom_boxplot(aes(x = Length, y = Dist,  color = Group), alpha = 0.5, outlier.shape = NA) + # Box plot without outlier dots
  labs(title = "VR1", 
       x = "Length", 
       y = "MAST") +
  ylim(c(min(mast_long$Dist), max(mast_long$Dist))) +
  stat_compare_means(aes(x = Length, y = Dist,  color = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE) + 
  theme_minimal() +
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "MAST", color = "Guides") +  # Change the Y-axis label here
  
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        title = element_text(size = 20, face = "bold")) + # Increase X-axis label size
  coord_cartesian(clip = "off")
V1

V2 = ggplot(mast_long %>% filter(VR == "V2")) +
  geom_boxplot(aes(x = Length, y = Dist,  color = Group), alpha = 0.5, outlier.shape = NA) + # Box plot without outlier dots
  labs(title = "VR2", 
       x = "Length", 
       y = "MAST") +
  ylim(c(min(mast_long$Dist), max(mast_long$Dist))) +
  stat_compare_means(aes(x = Length, y = Dist,  color = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE) + 
  theme_minimal() +
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "MAST", color = "Guides") +  # Change the Y-axis label here
  
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        title = element_text(size = 20, face = "bold")) + # Increase X-axis label size
  coord_cartesian(clip = "off")
V2

V3 = ggplot(mast_long %>% filter(VR == "V3")) +
  geom_boxplot(aes(x = Length, y = Dist,  color = Group), alpha = 0.5, outlier.shape = NA) + # Box plot without outlier dots
  labs(title = "VR3", 
       x = "Length", 
       y = "MAST") +
  ylim(c(min(mast_long$Dist), max(mast_long$Dist))) +
  stat_compare_means(aes(x = Length, y = Dist,  color = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE) + 
  theme_minimal() +
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "MAST", color = "Guides") +  # Change the Y-axis label here
  
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        title = element_text(size = 20, face = "bold")) + # Increase X-axis label size
  coord_cartesian(clip = "off")
V3

V4 = ggplot(mast_long %>% filter(VR == "V4")) +
  geom_boxplot(aes(x = Length, y = Dist,  color = Group), alpha = 0.5, outlier.shape = NA) + # Box plot without outlier dots
  labs(title = "VR4", 
       x = "Length", 
       y = "MAST") +
  ylim(c(min(mast_long$Dist), max(mast_long$Dist))) +
  stat_compare_means(aes(x = Length, y = Dist,  color = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE) + 
  theme_minimal() +
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "MAST", color = "Guides") +  # Change the Y-axis label here
  
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        title = element_text(size = 20, face = "bold")) + # Increase X-axis label size
  coord_cartesian(clip = "off")
V4


V5 = ggplot(mast_long %>% filter(VR == "V5")) +
  geom_boxplot(aes(x = Length, y = Dist,  color = Group), alpha = 0.5, outlier.shape = NA) + # Box plot without outlier dots
  labs(title = "VR5", 
       x = "Length", 
       y = "MAST") +
  ylim(c(min(mast_long$Dist), max(mast_long$Dist))) +
  stat_compare_means(aes(x = Length, y = Dist,  color = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE) + 
  theme_minimal() +
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "MAST", color = "Guides") +  # Change the Y-axis label here
  
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        title = element_text(size = 20, face = "bold")) + # Increase X-axis label size
  coord_cartesian(clip = "off")
V5

V6 = ggplot(mast_long %>% filter(VR == "V6")) +
  geom_boxplot(aes(x = Length, y = Dist,  color = Group), alpha = 0.5, outlier.shape = NA) + # Box plot without outlier dots
  labs(title = "VR6", 
       x = "Length", 
       y = "MAST") +
  ylim(c(min(mast_long$Dist), max(mast_long$Dist))) +
  scale_x_discrete(limits = c("L50", "L100", "L200", "L300", "L400", "L500")) +
  stat_compare_means(aes(x = Length, y = Dist,  color = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE) + 
  theme_minimal() +
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "MAST", color = "Guides") +  # Change the Y-axis label here
  
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        title = element_text(size = 20, face = "bold")) + # Increase X-axis label size
  coord_cartesian(clip = "off")
V6

V7 = ggplot(mast_long %>% filter(VR == "V7")) +
  geom_boxplot(aes(x = Length, y = Dist,  color = Group), alpha = 0.5, outlier.shape = NA) + # Box plot without outlier dots
  labs(title = "VR7", 
       x = "Length", 
       y = "MAST") +
  ylim(c(min(mast_long$Dist), max(mast_long$Dist))) +
  scale_x_discrete(limits = c("L50", "L100", "L200", "L300", "L400", "L500")) +
  stat_compare_means(aes(x = Length, y = Dist,  color = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE) + 
  theme_minimal() +
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "MAST", color = "Guides") +  # Change the Y-axis label here
  
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        title = element_text(size = 20, face = "bold")) + # Increase X-axis label size
  coord_cartesian(clip = "off")
V7

V8 = ggplot(mast_long %>% filter(VR == "V8")) +
  geom_boxplot(aes(x = Length, y = Dist,  color = Group), alpha = 0.5, outlier.shape = NA) + # Box plot without outlier dots
  labs(title = "VR8", 
       x = "Length", 
       y = "MAST") +
  ylim(c(min(mast_long$Dist), max(mast_long$Dist))) +
  scale_x_discrete(limits = c("L50", "L100", "L200", "L300", "L400", "L500")) +
  stat_compare_means(aes(x = Length, y = Dist,  color = Group), method = "wilcox.test", paired = FALSE, label = "p.signif", show.legend = FALSE) + 
  theme_minimal() +
  scale_color_manual(values = c("#003F5C", "#FFA600"),
                     labels = c("-", "+")) +
  labs(x = "", y = "MAST", color = "Guides") +  # Change the Y-axis label here
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        title = element_text(size = 20, face = "bold")) + # Increase X-axis label size 
  coord_cartesian(clip = "off")
V8

# Extract the legend to use
legend <- get_legend(V1 + theme(legend.position = "bottom"))

#Remove the legends from all the others
V1 = V1 + theme(legend.position = "none")
V2 = V2 + theme(legend.position = "none")
V3 = V3 + theme(legend.position = "none")
V4 = V4 + theme(legend.position = "none")
V5 = V5 + theme(legend.position = "none")
V6 = V6 + theme(legend.position = "none")
V7 = V7 + theme(legend.position = "none")
V8 = V8 + theme(legend.position = "none")
top_row <- (V1 + V2 + V3 + V4)
bottom_row <- (V5 + V6 + V7 + V8) 

# Final layout with the legend placed at the bottom
final_plot <- (top_row / bottom_row) / legend + plot_layout(heights = c(1, 1, 0.1))
final_plot 

setwd(dir.figures)
ggsave(filename = "length_vr_mast_plots.pdf", height = 500, width = 300, units = "mm")


################################################################################
# Clean Up
################################################################################
ls()
rm(list = c(ls(pattern = "^bp_"), ls(pattern = "s[0-9]"), "aics", "bics", "legend", "top_row", "bottom_row", ls(pattern = "V[0-9]"), ls(pattern = "rf_")))
rm(list = c("s_dists", "mast", "mast_long", "MATCHES", ls(pattern = "col")))
rm(list = c("final_plot", "micro_mean_size", "micro_total", "pearson_cor_hvr"))
ls()
################################################################################
# SIM 3: 
# 
################################################################################
setwd(dir.sim3)

################################################################################
# Global Distances: 100 Random pairwise comparisons 
################################################################################
# Read in the distance metrics
CONTROLS = c("VFull", "VFull_guides")

## 1 RF
ov_rf = as_tibble(readRDS("ov_rf_random_pairs.rds")) %>% rename_with(~ str_remove(., "_rf")) %>% select(!all_of(CONTROLS)) %>%  print(Inf)
labels = colnames(ov_rf)
labels

## 2 JRF
ov_jrf = as_tibble(readRDS("ov_jrf_random_pairs.rds")) %>% rename_with(~ str_remove(., "_jrf"))%>% select(!all_of(CONTROLS)) %>% print(Inf)
if(setequal(colnames(ov_jrf), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 3 Nye Distance
ov_nye = as_tibble(readRDS("ov_nye_random_pairs.rds")) %>% rename_with(~ str_remove(., "_nye")) %>% select(!all_of(CONTROLS))%>% print(Inf)
if(setequal(colnames(ov_nye), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 4 MSD
ov_msd = as_tibble(readRDS("ov_msd_random_pairs.rds")) %>% rename_with(~ str_remove(., "_msd")) %>% select(!all_of(CONTROLS))%>% print(Inf)
if(setequal(colnames(ov_msd), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 5 PID
ov_pid = as_tibble(readRDS("ov_pid_random_pairs.rds"))  %>% rename_with(~ str_remove(., "_pid"))%>% select(!all_of(CONTROLS))%>% print(Inf)
if(setequal(colnames(ov_pid), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 6 MCD
ov_mcd = as_tibble(readRDS("ov_mcd_random_pairs.rds"))  %>% rename_with(~ str_remove(., "_mcd"))%>% select(!all_of(CONTROLS))%>% print(Inf)
if(setequal(colnames(ov_mcd), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 7 Path
ov_path = as_tibble(readRDS("ov_path_random_pairs.rds"))%>% rename_with(~ str_remove(., "_path")) %>% select(!all_of(CONTROLS))%>% print(Inf)
if(setequal(colnames(ov_path), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 8 Quartet Distance
ov_quartet = as_tibble(readRDS("ov_quartet_normalized_random_pairs.rds")) %>% rename_with(~ str_remove(., "_quartet")) %>% rename_with(~ str_remove(., "_control")) %>% select(!all_of(CONTROLS))%>% print(Inf)
if(setequal(colnames(ov_quartet), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames"); setdiff(labels, colnames(ov_quartet))}

X = 12; Y = 18; TITLE_SIZE = 14
p_rf = plot_sim3_global(df = ov_rf, ylabel = "RF", x_size = X, y_size = Y, jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_jrf = plot_sim3_global(df = ov_jrf, ylabe = "JRF", x_size = X, y_size = Y, jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_nye = plot_sim3_global(df = ov_nye, ylabe = "NYE", x_size = X, y_size = Y, jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_mcd = plot_sim3_global(df = ov_mcd, ylabe = "MCD", x_size = X, y_size = Y, jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_msd = plot_sim3_global(df = ov_msd, ylabe = "MSD", x_size = X, y_size = Y, jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_path = plot_sim3_global(df = ov_path, ylabe = "PATH",x_size = X, y_size = Y,  y_val_sig_bars = c(69000), jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_pid = plot_sim3_global(df = ov_pid, ylabe = "PID", x_size = X, y_size = Y, y_val_sig_bars = c(0.32, 0.34, 0.30, 0.46, 0.48, 0.5, 0.40), jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_quar = plot_sim3_global(df = ov_quartet, ylabe = "QUAR", x_size = X, y_size = Y, y_val_sig_bars = c(0.25), jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)

# Extract the legend to use
legend <- get_legend(p_rf + theme(legend.position = "bottom"))

#Remove the legends from all the others
p_rf = p_rf + theme(legend.position = "none")
p_jrf = p_jrf + theme(legend.position = "none")
p_nye = p_nye + theme(legend.position = "none")
p_mcd = p_mcd + theme(legend.position = "none")
p_msd = p_msd + theme(legend.position = "none")
p_path = p_path + theme(legend.position = "none")
p_pid = p_pid + theme(legend.position = "none")
p_quar = p_quar + theme(legend.position = "none")
top_row <- (p_rf + p_jrf + p_nye + p_msd)
bottom_row <- (p_pid + p_mcd + p_path + p_quar) 

# Final layout with the legend placed at the bottom
final_plot <- (top_row / bottom_row) / legend + plot_layout(heights = c(1, 1, 0.1))
final_plot 

setwd(dir.figures)
ggsave(filename = "overlapping_vr_global_plots.pdf", height = 10, width = 7, units = "in")

X = 24; Y = 18; TITLE_SIZE = 14
p_rf = plot_sim3_global(df = ov_rf, ylabel = "RF", x_size = X, y_size = Y, jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
ggsave(filename = "overlapping_vr_global_RF_only.pdf", p_rf, height = 150, width = 150, units = "mm")

################################################################################
# MicroStructure: MAST 100 random pairwise comparisons
################################################################################
setwd(dir.sim3)
ov_mast = as_tibble(readRDS("ov_mast_random_pairs.rds")) %>% select(!all_of(CONTROLS)) %>% rename_with(~ str_remove(., "_rf"))%>%  print(Inf)
ov_mast

X = 24; Y = 18; TITLE_SIZE = 14
p_mast = plot_sim3_global_increasing(df = ov_mast, ylabe = "MAST", c(1500), x_size = X, y_size = Y, title_size = TITLE_SIZE)
p_mast
setwd(dir.figures)
ggsave(filename = "ov_mast_plots.pdf", height = 150, width = 150, units = "mm")

top_row <- p_rf + p_jrf + p_nye
middle_row <- p_mcd + p_msd + p_path
bottom_row <- p_pid + p_quar + p_mast

# Combine rows into a 3x3 grid
final_plot <- (top_row / middle_row / bottom_row) + 
  plot_layout(heights = c(1, 1, 1))  # Equal row heights
final_plot


setwd(dir.figures)
ggsave(filename = "ov_global_mast_plots.pdf", height = 300, width = 500, units = "mm")



################################################################################
# Nano Distances: 100 Random pairwise comparisons 
################################################################################

setwd(dir = dir.sim3)
s2 = as_tibble(readRDS("ov_nanoclusters_random_pairs_2.rds")) %>% rename_with(~ str_remove(., "_rf")) %>% select(!all_of(CONTROLS)) %>%  print(Inf)
s2
s3 = as_tibble(readRDS("ov_nanoclusters_random_pairs_3.rds")) %>% rename_with(~ str_remove(., "_rf"))%>% select(!all_of(CONTROLS)) %>%  print(Inf)
s3
s4 = as_tibble(readRDS("ov_nanoclusters_random_pairs_4.rds")) %>% rename_with(~ str_remove(., "_rf"))%>% select(!all_of(CONTROLS)) %>%  print(Inf)
s4
s5 = as_tibble(readRDS("ov_nanoclusters_random_pairs_5.rds"))%>% rename_with(~ str_remove(., "_rf"))%>% select(!all_of(CONTROLS)) %>%  print(Inf)
s5
s6 = as_tibble(readRDS("ov_nanoclusters_random_pairs_6.rds"))%>% rename_with(~ str_remove(., "_rf"))%>% select(!all_of(CONTROLS)) %>%  print(Inf)
s6
s7 = as_tibble(readRDS("ov_nanoclusters_random_pairs_7.rds"))%>% rename_with(~ str_remove(., "_rf"))%>% select(!all_of(CONTROLS)) %>%  print(Inf)
s7
s8 = as_tibble(readRDS("ov_nanoclusters_random_pairs_8.rds"))%>% rename_with(~ str_remove(., "_rf"))%>% select(!all_of(CONTROLS)) %>%  print(Inf)
s8
s9 = as_tibble(readRDS("ov_nanoclusters_random_pairs_9.rds"))%>% rename_with(~ str_remove(., "_rf"))%>% select(!all_of(CONTROLS)) %>%  print(Inf)
s9
s10 = as_tibble(readRDS("ov_nanoclusters_random_pairs_10.rds"))%>% rename_with(~ str_remove(., "_rf"))%>% select(!all_of(CONTROLS)) %>%  print(Inf)
s10
MAX = max(c(max(s2), max(s3), max(s4), max(s5), max(s6), max(s7), max(s8), max(s9), max(s10))) + 0.1
MAX
MIN = min(c(min(s2), min(s3), min(s4), min(s5), min(s6), min(s7), min(s8), min(s9), min(s10))) - 0.1
MIN

X = 16; Y = 18; TITLE_SIZE = 18
p_s2 = plot_sim3_local(df = s2, ylabel = "S2", y_min = MIN, y_max = MAX, x_size = X, title_size = TITLE_SIZE)
p_s2

p_s3 = plot_sim3_local(df = s3, ylabel = "S3", y_min = MIN, y_max = MAX, y_val_sig_bars = c(0.5), x_size = X, title_size = TITLE_SIZE)
p_s3

p_s4 = plot_sim3_local(df = s4, ylabel = "S4", y_min = MIN, y_max = MAX, x_size = X, title_size = TITLE_SIZE)
p_s4

p_s5 = plot_sim3_local(df = s5, ylabel = "S5",y_val_sig_bars = c(0.6), y_min = MIN, y_max = MAX, x_size = X, title_size = TITLE_SIZE)
p_s5

p_s6 = plot_sim3_local(df = s6, ylabel = "S6", y_min = MIN, y_max = MAX, x_size = X, title_size = TITLE_SIZE)
p_s6

p_s7 = plot_sim3_local(df = s7, ylabel = "S7",y_val_sig_bars = c(0.6), y_min = MIN, y_max = MAX, x_size = X, title_size = TITLE_SIZE)
p_s7

p_s8 = plot_sim3_local(df = s8, ylabel = "S8",y_val_sig_bars = c(0.55), y_min = MIN, y_max = MAX, x_size = X, title_size = TITLE_SIZE)
p_s8

p_s9 = plot_sim3_local(df = s9, ylabel = "S9",y_val_sig_bars = c(0.3), y_min = MIN, y_max = MAX, x_size = X, title_size = TITLE_SIZE)
p_s9

p_s10 = plot_sim3_local(df = s10, ylabel = "S10", y_min = MIN, y_max = MAX, x_size = X, title_size = TITLE_SIZE)
p_s10


# Extract the legend to use
legend <- get_legend(p_s2 + theme(legend.position = "bottom"))

#Remove the legends from all the others
p_s2 = p_s2 + theme(legend.position = "none")
p_s3 = p_s3 + theme(legend.position = "none")
p_s4 = p_s4 + theme(legend.position = "none")
p_s5 = p_s5 + theme(legend.position = "none")
p_s6 = p_s6 + theme(legend.position = "none")
p_s7 = p_s7 + theme(legend.position = "none")
p_s8 = p_s8 + theme(legend.position = "none")
p_s9 = p_s9 + theme(legend.position = "none")
p_s10 = p_s10 + theme(legend.position = "none")

top_row <- p_s2 + p_s3 + p_s4
middle_row <- p_s5 + p_s6 + p_s7
bottom_row <- p_s8 + p_s9 + p_s10

# Combine rows into a 3x3 grid
final_plot <- (top_row / middle_row / bottom_row) /legend + 
  plot_layout(heights = c(1, 1, 1, 0.1))  # Equal row heights
final_plot

setwd(dir.figures)
ggsave(filename = "ov_local_plots.pdf", height = 300, width = 500, units = "mm")

X = 24; Y = 18; TITLE_SIZE = 14

p_s2 = plot_sim3_local_increasing(df = s2, ylabel = "S2", y_min = MIN, y_max = MAX, x_size = X, title_size = TITLE_SIZE)
p_s2
ggsave(filename = "ov_s2_plots.pdf", p_s2, height = 150, width = 150, units = "mm")

################################################################################
# Clean Up
################################################################################
rm(list = c(ls(pattern = "^ov_"), ls(pattern = "s[0-9]")))


################################################################################
# SIM 4: Non-overlapping regions
################################################################################
setwd(dir.sim4)
################################################################################
# Global Distances: 100 Random pairwise comparisons 
################################################################################
# Read in the distance metrics
# Read in the distance metrics
## 1 RF
CONTROLS = c("VFull", "VFull_guides")
nonov_rf = as_tibble(readRDS("nonov_rf_random_pairs.rds"))  %>% select(!all_of(CONTROLS))%>% rename_with(~ str_remove(., "_rf")) %>% print(Inf)
labels = colnames(nonov_rf)
labels

## 2 JRF
nonov_jrf = as_tibble(readRDS("nonov_jrf_random_pairs.rds")) %>% select(!all_of(CONTROLS)) %>%  print(Inf)
if(setequal(colnames(nonov_jrf), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 3 Nye Distance
nonov_nye = as_tibble(readRDS("nonov_nye_random_pairs.rds")) %>% select(!all_of(CONTROLS)) %>%  print(Inf)
if(setequal(colnames(nonov_nye), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 4 MSD
nonov_msd = as_tibble(readRDS("nonov_msd_random_pairs.rds"))%>% select(!all_of(CONTROLS)) %>%  print(Inf)
if(setequal(colnames(nonov_msd), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 5 PID
nonov_pid = as_tibble(readRDS("nonov_pid_random_pairs.rds"))  %>% select(!all_of(CONTROLS)) %>%  print(Inf)
if(setequal(colnames(nonov_pid), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 6 MCD
nonov_mcd = as_tibble(readRDS("nonov_mcd_random_pairs.rds")) %>% select(!all_of(CONTROLS)) %>%  print(Inf)
if(setequal(colnames(nonov_mcd), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 7 Path
nonov_path = as_tibble(readRDS("nonov_path_random_pairs.rds")) %>% select(!all_of(CONTROLS)) %>%  print(Inf)
if(setequal(colnames(nonov_path), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 8 Quartet Distance
nonov_quartet = as_tibble(readRDS("nonov_quartet_normalized_random_pairs.rds"))%>% select(!VFull_control) %>% select(!VFull_guides) %>% rename_with(~ str_remove(., "_quartet"))%>%  print(Inf)
if(setequal(colnames(nonov_quartet), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}


X = 12; Y = 18; TITLE_SIZE = 14; P = 0.05
p_rf = plot_sim4_global(df = nonov_rf, ylabel = "RF", x_size = X, y_size = Y, thresh = P, jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_rf
p_jrf = plot_sim4_global(df = nonov_jrf, ylabe = "JRF", x_size = X, y_size = Y, thresh = P,  jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_jrf
p_nye = plot_sim4_global(df = nonov_nye, ylabel = "NYE", x_size = X, y_size = Y, thresh = P,  jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_nye
p_mcd = plot_sim4_global(df = nonov_mcd, ylabel = "MCD", x_size = X, y_size = Y, thresh = P,  jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_mcd
p_msd = plot_sim4_global(df = nonov_msd, ylabe = "MSD",  x_size = X, y_size = Y, thresh = P, c(70000, 90000, 80000),  jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_msd
p_path = plot_sim4_global(df = nonov_path, ylabe = "PATH",  x_size = X, y_size = Y, thresh = P, c(70000),  jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_path
p_pid = plot_sim4_global(df = nonov_pid, ylabe = "PID",  x_size = X, y_size = Y, thresh = P, c(0.4),  jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_pid
p_quar = plot_sim4_global(df = nonov_quartet, ylabe = "QUAR", x_size = X, y_size = Y, thresh = P, c(0.2, 0.2),  jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_quar

# Extract the legend to use
legend <- get_legend(p_rf + theme(legend.position = "bottom"))

#Remove the legends from all the others
p_rf = p_rf + theme(legend.position = "none")
p_jrf = p_jrf + theme(legend.position = "none")
p_nye = p_nye + theme(legend.position = "none")
p_mcd = p_mcd + theme(legend.position = "none")
p_msd = p_msd + theme(legend.position = "none")
p_path = p_path + theme(legend.position = "none")
p_pid = p_pid + theme(legend.position = "none")
p_quar = p_quar + theme(legend.position = "none")
top_row <- (p_rf + p_jrf + p_nye + p_msd)
bottom_row <- (p_pid + p_mcd + p_path + p_quar) 

# Final layout with the legend placed at the bottom
final_plot <- (top_row / bottom_row) / legend + plot_layout(heights = c(1, 1, 0.1))
final_plot 

setwd(dir.figures)
ggsave(filename = "nonoverlapping_vr_global_plots.pdf", height = 10, width = 7, units = "in")

X = 24; Y = 18; TITLE_SIZE = 14; P = 0.05
p_rf = plot_sim4_global(df = nonov_rf, ylabel = "RF", x_size = X, y_size = Y, thresh = P, jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
ggsave(filename = "nonov_rf_only.pdf", p_rf, height = 150, width = 150, units = "mm")


################################################################################
# Nano Distances: 100 Random pairwise comparisons 
################################################################################
setwd(dir.sim4)
X = 10; Y = 18; TITLE_SIZE = 14; P = 0.05
CONTROLS = c("VFull_guides", "VFull")
s2 = as_tibble(readRDS("nonov_nanostructure_2_random_pairs.rds")) %>% select(!all_of(CONTROLS))  
s2

s3 = as_tibble(readRDS("nonov_nanostructure_3_random_pairs.rds")) %>% select(!all_of(CONTROLS))  
s3

s4 = as_tibble(readRDS("nonov_nanostructure_4_random_pairs.rds")) %>% select(!all_of(CONTROLS))  
s4

s5 = as_tibble(readRDS("nonov_nanostructure_5_random_pairs.rds")) %>% select(!all_of(CONTROLS))  
s5

s6 = as_tibble(readRDS("nonov_nanostructure_6_random_pairs.rds")) %>% select(!all_of(CONTROLS))  
s6

s7 = as_tibble(readRDS("nonov_nanostructure_7_random_pairs.rds")) %>% select(!all_of(CONTROLS))  
s7

s8 = as_tibble(readRDS("nonov_nanostructure_8_random_pairs.rds")) %>% select(!all_of(CONTROLS))  
s8

s9 = as_tibble(readRDS("nonov_nanostructure_9_random_pairs.rds")) %>% select(!all_of(CONTROLS))  
s9

s10 = as_tibble(readRDS("nonov_nanostructure_10_random_pairs.rds")) %>% select(!all_of(CONTROLS))  
s10

MAX = max(c(max(s2), max(s3), max(s4), max(s5), max(s6), max(s7), max(s8), max(s9), max(s10))) + 0.1
MAX
MIN = min(c(min(s2), min(s3), min(s4), min(s5), min(s6), min(s7), min(s8), min(s9), min(s10))) - 0.1
MIN
p_s2 = plot_sim4_local(df = s2, ylabel = "S2", y_min = MIN, y_max = MAX, x_size = X, y_size = Y, thresh = P, jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_s2

p_s3 = plot_sim4_local(df = s3, ylabel = "S3", y_min = MIN, y_max = MAX, x_size = X, y_size = Y, thresh = P, jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_s3

p_s4 = plot_sim4_local(df = s4, ylabel = "S4", y_min = MIN, y_max = MAX, c(0.7), x_size = X, y_size = Y, thresh = P, jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_s4

p_s5 = plot_sim4_local(df = s5, ylabel = "S5",y_val_sig_bars = c(0.6), y_min = MIN, y_max = MAX, x_size = X, y_size = Y, thresh = P, jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_s5

p_s6 = plot_sim4_local(df = s6, ylabel = "S6", y_min = MIN, y_max = MAX, x_size = X, y_size = Y, thresh = P, jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_s6

p_s7 = plot_sim4_local(df = s7, ylabel = "S7", y_min = MIN, y_max = MAX, x_size = X, y_size = Y, thresh = P, jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_s7

p_s8 = plot_sim4_local(df = s8, ylabel = "S8",y_val_sig_bars = c(0.55), y_min = MIN, y_max = MAX, x_size = X, y_size = Y, thresh = P, jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_s8

p_s9 = plot_sim4_local(df = s9, ylabel = "S9",y_val_sig_bars = c(0.3), y_min = MIN, y_max = MAX, x_size = X, y_size = Y, thresh = P, jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_s9

p_s10 = plot_sim4_local(df = s10, ylabel = "S10", y_min = MIN, y_max = MAX, x_size = X, y_size = Y, thresh = P, jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_s10

# Extract the legend to use
legend <- get_legend(p_s2 + theme(legend.position = "bottom"))

#Remove the legends from all the others
p_s2 = p_s2 + theme(legend.position = "none")
p_s3 = p_s3 + theme(legend.position = "none")
p_s4 = p_s4 + theme(legend.position = "none")
p_s5 = p_s5 + theme(legend.position = "none")
p_s6 = p_s6 + theme(legend.position = "none")
p_s7 = p_s7 + theme(legend.position = "none")
p_s8 = p_s8 + theme(legend.position = "none")
p_s9 = p_s9 + theme(legend.position = "none")
p_s10 = p_s10 + theme(legend.position = "none")


top_row <- p_s2 + p_s3 + p_s4
middle_row <- p_s5 + p_s6 + p_s7
bottom_row <- p_s8 + p_s9 + p_s10

# Combine rows into a 3x3 grid
final_plot <- (top_row / middle_row / bottom_row) / legend + plot_layout(heights = c(1, 1, 1, 0.1))
final_plot


setwd(dir.figures)
ggsave(filename = "nonov_local_plots.pdf", height = 300, width = 500, units = "mm")

X = 24; Y = 18; TITLE_SIZE = 14; P = 0.05
p_s2 = plot_sim4_local_increasing(df = s2, ylabel = "S2", x_size = X, y_size = Y, thresh = P, jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_s2
ggsave(filename = "nonov_s2_only.pdf", p_s2, height = 150, width = 150, units = "mm")

################################################################################
# MicroStructure: MAST 100 random pairwise comparisons
################################################################################
X = 12; Y = 18; P = 0.05
CONTROLS = c("VFull_guides", "VFull")
setwd(dir.sim4)
non_mast = as_tibble(readRDS("nonov_mast_random_pairs.rds"))  %>% select(!all_of(CONTROLS))  
non_mast
P = 0.05
p_mast = plot_sim4_global(df = non_mast, ylabe = "MAST", x_size = X, y_size = Y, thresh = P, c(3000, 2600, 2100))
p_mast

setwd(dir.figures)
ggsave(filename = "nonov_mast_plots.pdf", height = 200, width = 100, units = "mm")

X = 24; Y = 18; TITLE_SIZE = 14; P = 0.05
p_mast = plot_sim4_global_increasing(df = non_mast, ylabel = "MAST", x_size = X, y_size = Y, thresh = P, jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE, c(3000, 2600, 2100))
p_mast
ggsave(filename = "nonov_mast_only.pdf", p_mast, height = 150, width = 150, units = "mm")


################################################################################
# SIM 5: 
################################################################################
setwd(dir.sim5)
################################################################################
# Global Distances: 100 Random pairwise comparisons 
################################################################################
CONTROLS = c("VFull_guides", "VFull")

# Read in the distance metrics
## 1 RF
nd_rf = as_tibble(readRDS("nd_rf_random_pairs.rds")) %>% rename_with(~ str_remove(., "_rf")) %>% select(!all_of(CONTROLS)) %>% print(Inf)
labels = colnames(nd_rf)
labels

## 2 JRF
nd_jrf = as_tibble(readRDS("sim_analysis_nd_jrf.rds"))  %>% rename_with(~ str_remove(., "_rf")) %>% select(!all_of(CONTROLS)) %>% print(Inf)
if(setequal(colnames(nd_jrf), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 3 Nye Distance
nd_nye = as_tibble(readRDS("sim_analysis_nd_nye.rds"))  %>% rename_with(~ str_remove(., "_rf")) %>% select(!all_of(CONTROLS)) %>% print(Inf)
if(setequal(colnames(nd_nye), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 4 MSD
nd_msd = as_tibble(readRDS("sim_analysis_nd_msd.rds"))  %>% rename_with(~ str_remove(., "_rf")) %>% select(!all_of(CONTROLS)) %>% print(Inf)
if(setequal(colnames(nd_msd), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 5 PID
nd_pid = as_tibble(readRDS("sim_analysis_nd_pid.rds"))  %>% rename_with(~ str_remove(., "_rf")) %>% select(!all_of(CONTROLS)) %>% print(Inf)
if(setequal(colnames(nd_pid), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 6 MCD
nd_mcd = as_tibble(readRDS("sim_analysis_nd_mcd.rds"))  %>% rename_with(~ str_remove(., "_rf")) %>% select(!all_of(CONTROLS)) %>% print(Inf)
if(setequal(colnames(nd_mcd), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 7 Path
nd_path = as_tibble(readRDS("sim_analysis_nd_path.rds"))  %>% rename_with(~ str_remove(., "_rf")) %>% select(!all_of(CONTROLS)) %>% print(Inf)
if(setequal(colnames(nd_path), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 8 Quartet Distance
nd_quartet = as_tibble(readRDS("nd_quartet_normalized_random_pairs.rds"))  %>% rename_with(~ str_remove(., "_quartet")) %>% select(!all_of(CONTROLS)) %>% print(Inf)
if(setequal(colnames(nd_quartet), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}


X = 12; Y = 18; TITLE_SIZE = 14; P = 0.05
p_rf = plot_sim5_global(df = nd_rf, ylabel = "RF", x_size = X, y_size = Y, thresh = P, jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE, c(0.4))
p_rf
p_jrf = plot_sim5_global(df = nd_jrf, ylabe = "JRF", x_size = X, y_size = Y, thresh = P,  jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_jrf
p_nye = plot_sim5_global(df = nd_nye, ylabel = "NYE", x_size = X, y_size = Y, thresh = P,  jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_nye
p_mcd = plot_sim5_global(df = nd_mcd, ylabel = "MCD", x_size = X, y_size = Y, thresh = P,  jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE, c(0.28, 0.28))
p_mcd
p_msd = plot_sim5_global(df = nd_msd, ylabe = "MSD",  x_size = X, y_size = Y, thresh = P, c(90000, 70000, 50000),  jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_msd
p_path = plot_sim5_global(df = nd_path, ylabe = "PATH",  x_size = X, y_size = Y, thresh = P, c(55000, 65000, 60000, 75000, 70000, 60000, 55000, 55000),  jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_path
p_pid = plot_sim5_global(df = nd_pid, ylabe = "PID",  x_size = X, y_size = Y, thresh = P, c(0.4),  jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_pid
p_quar = plot_sim5_global(df = nd_quartet, ylabe = "QUAR", x_size = X, y_size = Y, thresh = P, c(0.2, 0.15, 0.17),  jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_quar

# Extract the legend to use
legend <- get_legend(p_rf + theme(legend.position = "bottom"))

#Remove the legends from all the others
p_rf = p_rf + theme(legend.position = "none")
p_jrf = p_jrf + theme(legend.position = "none")
p_nye = p_nye + theme(legend.position = "none")
p_mcd = p_mcd + theme(legend.position = "none")
p_msd = p_msd + theme(legend.position = "none")
p_path = p_path + theme(legend.position = "none")
p_pid = p_pid + theme(legend.position = "none")
p_quar = p_quar + theme(legend.position = "none")
top_row <- (p_rf + p_jrf + p_nye + p_msd)
bottom_row <- (p_pid + p_mcd + p_path + p_quar) 

# Final layout with the legend placed at the bottom
final_plot <- (top_row / bottom_row) / legend + plot_layout(heights = c(1, 1, 0.1))
final_plot 

setwd(dir.figures)
ggsave(filename = "ndmix_global_plots.pdf", height = 11, width = 7, units = "in")


################################################################################
# Nano Distances: 100 Random pairwise comparisons 
################################################################################
setwd(dir.sim)
s2 = as_tibble(readRDS("nanostructure_nd_2_random_pairs.rds")) %>% select(!all_of(CONTROLS))
s2
s3 = as_tibble(readRDS("nanostructure_nd_3_random_pairs.rds"))%>% select(!all_of(CONTROLS))
s3
s4 = as_tibble(readRDS("nanostructure_nd_4_random_pairs.rds"))%>% select(!all_of(CONTROLS))
s4
s5 = as_tibble(readRDS("nanostructure_nd_5_random_pairs.rds"))%>% select(!all_of(CONTROLS))
s5
s6 = as_tibble(readRDS("nanostructure_nd_6_random_pairs.rds"))%>% select(!all_of(CONTROLS))
s6
s7 = as_tibble(readRDS("nanostructure_nd_7_random_pairs.rds"))%>% select(!all_of(CONTROLS))
s7
s8 = as_tibble(readRDS("nanostructure_nd_8_random_pairs.rds"))%>% select(!all_of(CONTROLS))
s8
s9 = as_tibble(readRDS("nanostructure_nd_9_random_pairs.rds"))%>% select(!all_of(CONTROLS))
s9
s10 = as_tibble(readRDS("nanostructure_nd_10_random_pairs.rds"))%>% select(!all_of(CONTROLS))
s10 %>% print(Inf)


MAX = max(c(max(s2), max(s3), max(s4), max(s5), max(s6), max(s7), max(s8), max(s9), max(s10))) + 0.1
MAX
MIN = min(c(min(s2), min(s3), min(s4), min(s5), min(s6), min(s7), min(s8), min(s9), min(s10))) - 0.1
MIN
p_s2 = plot_sim5_local(df = s2, ylabel = "S2", y_min = MIN, y_max = MAX, x_size = X, y_size = Y, thresh = P, jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_s2

p_s3 = plot_sim5_local(df = s3, ylabel = "S3", y_min = MIN, y_max = MAX, x_size = X, y_size = Y, thresh = P, jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_s3

p_s4 = plot_sim5_local(df = s4, ylabel = "S4", y_min = MIN, y_max = MAX, c(0.7), x_size = X, y_size = Y, thresh = P, jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_s4

p_s5 = plot_sim5_local(df = s5, ylabel = "S5",y_val_sig_bars = c(0.6), y_min = MIN, y_max = MAX, x_size = X, y_size = Y, thresh = P, jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_s5

p_s6 = plot_sim5_local(df = s6, ylabel = "S6", y_min = MIN, y_max = MAX, x_size = X, y_size = Y, thresh = P, jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE, c(0.6))
p_s6

p_s7 = plot_sim5_local(df = s7, ylabel = "S7", y_min = MIN, y_max = MAX, x_size = X, y_size = Y, thresh = P, jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_s7

p_s8 = plot_sim5_local(df = s8, ylabel = "S8",y_val_sig_bars = c(0.55), y_min = MIN, y_max = MAX, x_size = X, y_size = Y, thresh = P, jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_s8

p_s9 = plot_sim5_local(df = s9, ylabel = "S9",y_val_sig_bars = c(0.3), y_min = MIN, y_max = MAX, x_size = X, y_size = Y, thresh = P, jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_s9

p_s10 = plot_sim5_local(df = s10, ylabel = "S10", y_min = MIN, y_max = MAX, x_size = X, y_size = Y, thresh = P, jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_s10

# Extract the legend to use
legend <- get_legend(p_s2 + theme(legend.position = "bottom"))

#Remove the legends from all the others
p_s2 = p_s2 + theme(legend.position = "none")
p_s3 = p_s3 + theme(legend.position = "none")
p_s4 = p_s4 + theme(legend.position = "none")
p_s5 = p_s5 + theme(legend.position = "none")
p_s6 = p_s6 + theme(legend.position = "none")
p_s7 = p_s7 + theme(legend.position = "none")
p_s8 = p_s8 + theme(legend.position = "none")
p_s9 = p_s9 + theme(legend.position = "none")
p_s10 = p_s10 + theme(legend.position = "none")

# Prepare plot rows.
top_row <- p_s2 + p_s3 + p_s4
middle_row <- p_s5 + p_s6 + p_s7
bottom_row <- p_s8 + p_s9 + p_s10

# Combine rows into a 3x3 grid
final_plot <- (top_row / middle_row / bottom_row) / legend + plot_layout(heights = c(1, 1, 1, 0.1))
final_plot

#Save File
setwd(dir.figures)
ggsave(filename = "nd_local_plots.pdf", height = 300, width = 500, units = "mm")


################################################################################
# MicroStructure: MAST 100 random pairwise comparisons
################################################################################
setwd(dir.sim5)
CONTROLS = c("VFull_guides", "VFull")
mast = readRDS("nd_mast_random_pairs.rds") %>% select(!all_of(CONTROLS))  
mast

P = 0.05
p_mast = plot_sim5_global(df = mast, ylabe = "MAST", x_size = X, y_size = Y, thresh = P, c(2500))
p_mast

setwd(dir.figures)
ggsave(filename = "nd_mast_plots.pdf", height = 200, width = 100, units = "mm")

################################################################################
# EXP: August 20th, 2025
################################################################################
rm(list = ls())

library(phangorn)
library(tidyr)
library(tibble)
library(purrr)
library(ggpubr)
library(reshape2)
library(ggplot2)
library(dplyr)
library(patchwork)
library(stringr)
library(rstatix)

################################################################################
# SET DIRECTORIES
################################################################################
# Set Home
HOME = "/Users/arnoldhk/Desktop/Research/2025_HVR_Guide_Phylogenetic_Integration/"

# Make directories
setwd(HOME)
dir.scripts = file.path(HOME, "2025_hvr_guide_phylogenetic_integration/scripts/")
dir.ref = file.path(HOME, "reference_data/")
dir.out = file.path(HOME, "2025_hvr_guide_phylogenetic_integration/data_visulization/")
dir.sim2 = file.path(HOME, "2025_hvr_guide_phylogenetic_integration/out/sim_bp/")
dir.figures = file.path(HOME, "2025_hvr_guide_phylogenetic_integration/figuresNoSync/")
dir.sim3 = file.path(HOME, "2025_hvr_guide_phylogenetic_integration/out/sim_overlapping/")
dir.sim4 = file.path(HOME, "2025_hvr_guide_phylogenetic_integration/out/sim_nonoverlapping/")
dir.sim5 = file.path(HOME, "2025_hvr_guide_phylogenetic_integration/out/sim_nd_mix/")
dir.exp = file.path(HOME, "2025_hvr_guide_phylogenetic_integration/out/experimental_demo/")

# Get functions
setwd(dir.scripts)
source("sim_analysis_visualization_functions.R")

# Read in the distance metrics
## 1 RF
setwd(dir.out)
CONTROLS = c("control", "long_control")
exp_rf = as_tibble(readRDS("experimental_demo_rf_random_pairs.rds"))  %>% select(!all_of(CONTROLS))%>% rename ("short_guides" = short_guides_opt) %>% print(Inf)
labels = colnames(exp_rf)
labels

## 2 JRF
exp_jrf = as_tibble(readRDS("experimental_demo_jrf_random_pairs.rds")) %>% select(!any_of(CONTROLS)) %>%  print(Inf)
if(setequal(colnames(exp_jrf), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 3 Nye Distance
exp_nye = as_tibble(readRDS("experimental_demo_nye_random_pairs.rds")) %>% select(!any_of(CONTROLS)) %>%  print(Inf)
if(setequal(colnames(exp_nye), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 4 MSD
exp_msd = as_tibble(readRDS("experimental_demo_msd_random_pairs.rds"))%>% select(!any_of(CONTROLS)) %>%  print(Inf)
if(setequal(colnames(exp_msd), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 5 PID
exp_pid = as_tibble(readRDS("experimental_demo_pid_random_pairs.rds"))  %>% select(!any_of(CONTROLS)) %>%  print(Inf)
if(setequal(colnames(exp_pid), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 6 MCD
exp_mcd = as_tibble(readRDS("experimental_demo_mcd_random_pairs.rds")) %>% select(!any_of(CONTROLS)) %>%  print(Inf)
if(setequal(colnames(exp_mcd), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 7 Path
exp_path = as_tibble(readRDS("experimental_demo_path_random_pairings.rds")) %>% select(!any_of(CONTROLS)) %>%  print(Inf)
if(setequal(colnames(exp_path), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}

## 8 Quartet Distance
exp_quartet = as_tibble(readRDS("experimental_demo_quartet_normalized.rds"))%>% select(!any_of(CONTORLS)) %>% print(Inf)
if(setequal(colnames(exp_quartet), labels)){print("PASS: Column names are standardized")}else{print("FAIL: Standardize colnames")}


X = 12; Y = 18; TITLE_SIZE = 14; 
p_rf = plot_exp_global(df = exp_rf, ylabel = "RF", x_size = X, y_size = Y, jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_rf
p_jrf = plot_exp_global(df = exp_jrf, ylabe = "JRF", x_size = X, y_size = Y, jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_jrf
p_nye = plot_exp_global(df = exp_nye, ylabel = "NYE", x_size = X, y_size = Y,  jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_nye
p_mcd = plot_exp_global(df = exp_mcd, ylabel = "MCD", x_size = X, y_size = Y,  jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_mcd
p_msd = plot_exp_global(df = exp_msd, ylabe = "MSD",  x_size = X, y_size = Y,   jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_msd
p_path = plot_exp_global(df = exp_path, ylabe = "PATH",  x_size = X, y_size = Y,  jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_path
p_pid = plot_exp_global(df = exp_pid, ylabe = "PID",  x_size = X, y_size = Y, jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_pid
p_quar = plot_exp_global(df = exp_quartet, ylabe = "QUAR", x_size = X, y_size = Y,  jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_quar

setwd(dir.figures)
ggsave(filename = "exp_global_rf.pdf", p_rf, height = 150, width = 150, units = "mm")

# Extract the legend to use
legend <- get_legend(p_rf + theme(legend.position = "bottom"))

#Remove the legends from all the others
p_rf = p_rf + theme(legend.position = "none")
p_jrf = p_jrf + theme(legend.position = "none")
p_nye = p_nye + theme(legend.position = "none")
p_mcd = p_mcd + theme(legend.position = "none")
p_msd = p_msd + theme(legend.position = "none")
p_path = p_path + theme(legend.position = "none")
p_pid = p_pid + theme(legend.position = "none")
p_quar = p_quar + theme(legend.position = "none")

final_plot <-
  (p_rf + p_jrf + p_nye + p_msd +
     p_pid + p_mcd + p_path + p_quar) +
  plot_layout(ncol = 4, nrow = 2, guides = "collect") &
  theme(legend.position = "bottom")

final_plot

setwd(dir.figures)
ggsave(filename = "exp_global_plots.pdf", height = 10, width = 7, units = "in")



################################################################################
# Nano Distances: 100 Random pairwise comparisons 
################################################################################
setwd(dir.out)
X = 10; Y = 18; TITLE_SIZE = 14;
CONTROLS = c(CONTROLS, "s2_long_control")
s2 = as_tibble(readRDS("experimental_demo_subtrees_random_pairs.rds")) %>% select(!any_of(CONTROLS))  
s2

p_s2 = plot_exp_local(df = s2, ylabe = "S2", x_size = X, y_size = Y, jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_s2

setwd(dir.figures)
ggsave(filename = "exp_local_s2.pdf", p_s2,   height = 100, width = 100, units = "mm")
################################################################################
# MicroStructure: MAST 100 random pairwise comparisons
################################################################################
X = 12; Y = 18; 
setwd(dir.out)
CONTROLS = c(CONTROLS, "long_control_micro")
exp_mast = as_tibble(readRDS("experimental_demo_mast_random_pairs.rds"))  %>% select(!any_of(CONTROLS))  
exp_mast

p_mast = plot_exp_local(df = exp_mast, ylabe = "MAST", x_size = X, y_size = Y, jitter_size = 0.1, geom_alpha = 0.6, title_size = TITLE_SIZE)
p_mast

setwd(dir.figures)
ggsave(filename = "exp_mast_plots.pdf", height = 100, width = 100, units = "mm")



