#AUTHOR: ARNOLD
#DAY: Jan 13th, 2025
# calculate_entrpy.R
# The purpose of this script is to calculate the entropy along the 16S gene
# and to fit a loess function to the data by optimizing the span parameter.

# LIBRARIES
.libPaths(c(.libPaths(), "/home/micro/arnoldho/R/x86_64-pc-linux-gnu-library/4.1"))
library(phylotools)
library(dplyr)
library(purr)
library("PIMMA")
library(ggplot2)

# DIRECTORIES
HOME = "/nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration"
setwd(HOME)
HOME = "/Users/arnoldhk/Desktop/Research/2025_HVR_Guide_Phylogenetic_Integration/"


dir.scripts = file.path(HOME, "2025_hvr_guide_phylogenetic_integration/scripts/")
dir.ref = file.path(HOME, "reference_data/")
dir.figures = file.path(HOME, "2025_hvr_guide_phylogenetic_integration/figuresNoSync/")

# Read in files
setwd(dir.ref)
genes = read.fasta("ref.align")
head(genes)

# # Flatten the alignment into a sequence matrix
# N = stringr::str_length(genes[1,2])
# d = as.data.frame(matrix(ncol = N, nrow = nrow(genes)))
# colnames(d) <- paste0("col", seq_len(ncol(d)))
# 
# for(i in 1:nrow(d)){
#   print(i/N)
#   rownames(d)[i] = genes[i, "seq.name"]
#   cur_string = genes[i, "seq.text"]
#   for(j in 1:N){
#     d[i,j] = substring(cur_string, j, j)
#   }
# }

# Save the matrix
setwd(dir.ref)
#saveRDS(d, "flattened_alignment.rds", compress = TRUE)
d = readRDS("flattened_alignment.rds")
# Traspose so that columns are strains, rows are alignment columns
t_d = t(d)
colnames(t_d)

# See what possible bases are present.
bases = vector()
for(i in 1:nrow(t_d)){
  bases = unique(c(bases, names(table(t_d[i,]))))
}
bases # Possible characters are A, C, T, G, -, .

# Calculate entropy across each base position.
entropy = vector(length = nrow(t_d))
for(i in 1:length(entropy)){
  cur_bases = names(table(t_d[i,]))
  n_bases = 0
  if("A" %in% cur_bases){
    n_bases = sum(table(t_d[i,])[c("A")]) + n_bases
  }
  if("C" %in% cur_bases){
    n_bases = sum(table(t_d[i,])[c("C")])+ n_bases
  }
  
  if("G" %in% cur_bases){
    n_bases = sum(table(t_d[i,])[c("G")])+ n_bases
  }
  
  if("T" %in% cur_bases){
    n_bases = sum(table(t_d[i,])[c("T")])+ n_bases
  }
  
  (table(t_d[i,])[c("A", "T", "G", "C")])
  n_bases
  if(n_bases < nrow(t_d)*0){
    entropy[i] = NA
  }else{
    if("A" %in% cur_bases){
      p_A = sum(table(t_d[i,])[c("A")])/n_bases
      A = p_A*log2(p_A)
    }else{A = 0}
    if("C" %in% cur_bases){
      p_C = sum(table(t_d[i,])[c("C")])/n_bases
      C = p_C*log2(p_C)
      
    }else{C = 0}
    
    if("G" %in% cur_bases){
      p_G = sum(table(t_d[i,])[c("G")])/n_bases
      G = p_G*log2(p_G)
      
    }else{G = 0}
    
    if("T" %in% cur_bases){
      p_T = sum(table(t_d[i,])[c("T")])/n_bases
      T = p_T*log2(p_T)
      
    }else{T = 0}
    
  }
  entropy[i] = -1*(A + C + G + T)
}
entropy


entropy = data.frame("base" = seq(from = 1, to = length(entropy), by = 1), "entropy" = entropy)
head(entropy)


# Now remove those bases that are not in our e. coli model 
setwd(dir.ref)
ec = read.fasta("e_coli_ref.fasta")
ec
which(t_d[,"AE005174.Esch1125"] != "-")

entropy$ecoli = NA
table(t_d[,"AE005174.Esch1125"])
entropy[which(t_d[,"AE005174.Esch1125"] != "-"),"ecoli"] = "Base"
head(entropy)
entropy_filtered = entropy %>% filter(ecoli == "Base")
entropy_filtered$basepair = seq(from = 1, to = nrow(entropy_filtered), by = 1)
head(entropy_filtered)

# Define range of span values to test
span_values <- seq(from = 0.1, to = 1, by = 0.1)
mse_values <- numeric(length(span_values))

# Loop through different spans and compute MSE
for (i in seq_along(span_values)) {
  loess_fit <- loess(entropy ~ basepair, data = entropy_filtered, span = span_values[i])
  pred <- predict(loess_fit, entropy_filtered$basepair)
  mse_values[i] <- Metrics::mse(entropy_filtered$entropy, pred)
}

# Find the best span (minimizing MSE)
best_span <- span_values[which.min(mse_values)]
SPAN = best_span
print(paste("Best smoothing span:", best_span))

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
  label = c("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9"), # Labels
  color = c(col_R1, col_R2, col_R3, col_R4, col_R5, col_R6, col_R7, col_R8, col_R9) # Colors for regions
)



# Custom color palette
custom_palette <- setNames(regions$color, regions$label)

# Create the plot
p = ggplot(entropy_filtered, aes(x = basepair, y = entropy)) +
  geom_point(color = "black", size = 2) +  # Data points
  geom_smooth(method = "loess", span = SPAN, color = "black", fill = "lightgrey", se = TRUE) + # Smoothed line
  geom_rect(
    data = regions, aes(
      xmin = start, xmax = end, ymin = -2, ymax = 2.0, fill = label
    ),
    alpha = 0.2, inherit.aes = FALSE
  ) +
  scale_fill_manual(values = custom_palette) + # Apply custom colors
  labs(
    title = "Entropy Accross the 16S Gene",
    x = "Basepair Position",
    y = "Entropy",
    fill = "Regions"
  ) +
  scale_x_continuous(breaks = seq(min(entropy_filtered$basepair), max(entropy_filtered$basepair), by = 50)) +
  theme_minimal()
p

setwd(dir.figures)
ggsave(filename = "EntropyAlong16S.pdf", p, height = 20, width =40, units = "cm")

l = data.frame("Variable Region" = c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9"),
               "VR" = c(45.8, 224.0, 165.9, 272.0, 122.9, 139.4, 126.0, 208.6, 83.4),
               "CF" = c(96.1, 258.9, 175.1, 283.0, 149.9, 154.5, 132.5, 227.0, 113.6))

apply(l[,2:3], 2, sum)

# Fit a linear model
lm_fit <- lm(CF ~ VR, data = l)
lm_fit
lm_fit# Get R² and p-value
summary_fit <- summary(lm_fit)
r_squared <- summary_fit$r.squared
p_value <- summary_fit$coefficients[2,4]  # p-value for the slope

# Create the plot with regression line and R², p-value annotation
p = ggplot(l, aes(x = VR, y = CF)) +
  geom_point(size = 3, color = "black") +  # Scatter plot points
  geom_smooth(method = "lm", color = "red", se = FALSE) +  # Linear fit line
  annotate("text", x = 50, y = 250, 
           label = paste0("y = ", signif(as.vector(lm_fit$coefficients["VR"]), 3), 
                          "x + ",
                          signif(as.vector(lm_fit$coefficients["(Intercept)"]), 3), 
                          "\nR² = ", round(r_squared, 3), 
                          "\np = ", signif(p_value, 3)), 
           size = 5, hjust = 0) +
  ylim(c(50, 300)) + 
  geom_text(aes(label = Variable.Region), vjust = -1, size = 3.5) +  # Add labels above points
  labs(title = "Reported vs. Guide Sequence Lengths",
       x = "VR Guide Sequence Base Pairs", y = "VR Base Pairs (Vargas-Albores et al. 2017)") 
p
setwd(dir.figures)
ggsave(p, filename = "vr_vs_cr_length.png", width = 100, height = 100, units = "mm")
 
