plot_exp_global <- function(df, ylabel = "Dist", 
                            x_size = 20, y_size = 18, 
                            jitter_size = 0.2, geom_alpha = 0.6, 
                            title_size = 18) {
  
  # If df isn't already long with Dist + Group_Type, make it so
  if (!all(c("Dist", "Group_Type") %in% names(df))) {
    df <- df %>%
      tidyr::pivot_longer(cols = tidyselect::everything(),
                          names_to = "Simulation", values_to = "Dist") %>%
      dplyr::mutate(Group_Type = dplyr::if_else(stringr::str_detect(Simulation, "_guides"),
                                                "With Guides", "Without Guides")) %>%
      dplyr::select(Dist, Group_Type)
  }
  
  # lock factor order/mapping
  df <- df %>%
    dplyr::mutate(Group_Type = factor(Group_Type, levels = c("With Guides", "Without Guides")))
  
  color_map <- c("With Guides" = "#FFA600", "Without Guides" = "#003F5C")
  label_map <- c("With Guides" = "+", "Without Guides" = "-")
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Group_Type, y = Dist, color = Group_Type)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(size = jitter_size, width = 0.2, alpha = geom_alpha) +
    ggplot2::stat_summary(fun = mean, geom = "point", 
                          shape = 18, size = 4, color = "black", fill = "black") + # add group means
    ggpubr::stat_compare_means(method = "wilcox.test", 
                               comparisons = list(c("With Guides", "Without Guides")),
                               label = "p.format") +  # prints p-value
    ggplot2::scale_color_manual(values = color_map, labels = label_map, name = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = NULL, x = "", y = ylabel) +
    ggplot2::theme(
      axis.text.x  = element_blank(),
      axis.title.y = ggplot2::element_text(size = y_size),
      plot.title   = ggplot2::element_text(size = title_size, face = "bold")
    )
  
  return(p)
}

plot_exp_local <- function(df, ylabel = "Dist", 
                            x_size = 20, y_size = 18, 
                            jitter_size = 0.2, geom_alpha = 0.6, 
                            title_size = 18) {
  
  # If df isn't already long with Dist + Group_Type, make it so
  if (!all(c("Dist", "Group_Type") %in% names(df))) {
    df <- df %>%
      tidyr::pivot_longer(cols = tidyselect::everything(),
                          names_to = "Simulation", values_to = "Dist") %>%
      dplyr::mutate(Group_Type = dplyr::if_else(stringr::str_detect(Simulation, "_guides"),
                                                "With Guides", "Without Guides")) %>%
      dplyr::select(Dist, Group_Type)
  }
  
  # lock factor order/mapping
  df <- df %>%
    dplyr::mutate(Group_Type = factor(Group_Type, levels = c("With Guides", "Without Guides")))
  
  color_map <- c("With Guides" = "#FFA600", "Without Guides" = "#003F5C")
  label_map <- c("With Guides" = "+", "Without Guides" = "-")
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Group_Type, y = Dist, color = Group_Type)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(size = jitter_size, width = 0.2, alpha = geom_alpha) +
    ggplot2::stat_summary(fun = mean, geom = "point", 
                          shape = 18, size = 4, color = "black", fill = "black") + # add group means
    ggpubr::stat_compare_means(method = "wilcox.test", 
                               comparisons = list(c("With Guides", "Without Guides")),
                               label = "p.format") +  # prints p-value
    ggplot2::scale_color_manual(values = color_map, labels = label_map, name = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = NULL, x = "", y = ylabel) +
    ggplot2::theme(
      axis.text.x  =  ggplot2::element_text(size = x_size),
      axis.title.y = ggplot2::element_text(size = y_size),
      plot.title   = ggplot2::element_text(size = title_size, face = "bold")
    )
  
  return(p)
}

plot_sim5_local = function(df,  ylabel = "Dist", 
                           y_val_sig_bars = rep(1, times = 13), 
                           thresh = 0.05,
                           x_size = 20,
                           y_size = 18,
                           y_min = 0,
                           y_max = 1, 
                           jitter_size = 0.2, 
                           geom_alpha = 0.6, 
                           title_size = 18){
  
  df <- 
    df %>%
    pivot_longer(cols = everything(), names_to = "Simulation", values_to = "Dist") %>%
    mutate(Group_Type = ifelse(stringr::str_detect(Simulation, '_guides'), "With Guides", "Without Guides")) %>%
    mutate(Group = case_when(Simulation == "V1V2" ~ "V1V2-",
                             Simulation == "V1V2_guides" ~ "V1V2+",
                             Simulation == "V1V3" ~ "V1V3-",
                             Simulation == "V1V3_guides" ~ "V1V3+",
                             
                             Simulation == "V3V4" ~ "V3V4-",
                             Simulation == "V3V4_guides" ~ "V3V4+",
                             
                             Simulation == "V3V5" ~ "V3V5-",
                             Simulation == "V3V5_guides" ~ "V3V5+",
                             
                             Simulation == "V4" ~ "V4-",
                             Simulation == "V4_guides" ~ "V4+",
                             Simulation == "V4V5" ~ "V4V5-",
                             Simulation == "V4V5_guides" ~ "V4V5+",
                             
                             Simulation == "nd_mix" ~ "Mix-",
                             Simulation == "nd_mix_guides" ~ "Mix+",
                             TRUE ~ Simulation)) %>%
    select(!Simulation) %>%
    rename("Value" = Dist)
  df
  
  # Order by median value
  df <- df %>%
    group_by(Group) %>%
    mutate(median_value = median(Value)) %>%
    arrange(Group_Type, desc(median_value))
  df
  
  # Perform pairwise Wilcoxon Rank Sum tests within each Group_Type
  pairwise_tests = pairwise_wilcox_test(data = data.frame(df), 
                                        formula = Value ~ Group, 
                                        p.adjust.method = "fdr", 
                                        paired = FALSE)
  as_tibble(pairwise_tests) %>% print(n = 100)
  
  
  
  if(pairwise_tests %>% filter(p.adj > thresh) %>% nrow() > 0){
    pairwise_tests= pairwise_tests %>% filter(p.adj > thresh)
    pairwise_tests$y.position = y_val_sig_bars
    p = ggplot(df, aes(x = reorder(Group, -median_value), y = Value)) +
      geom_boxplot(aes(color = Group_Type)) +
      ylim(y_min, y_max) + 
      theme_minimal() +
      labs(title = ylabel,
           x = "",
           y = "") +
      theme(axis.text.x = element_text(size = x_size, angle = 45, hjust = 1),
            axis.title.y = element_text(size = y_size),
            plot.title = element_text(size = title_size, face = "bold")) +
      scale_color_manual(values = c("#FFA600", "#003F5C"),
                         labels = c("+", "-")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      geom_jitter(aes(color = Group_Type), size = jitter_size, width = 0.2, alpha = geom_alpha) +  # Add jitter points
      stat_pvalue_manual(pairwise_tests, label = "p.adj.signif", tip.length = 0.02) 
  }else{
    p = ggplot(df, aes(x = reorder(Group, -median_value), y = Value)) +
      geom_boxplot(aes(color = Group_Type)) +
      ylim(y_min, y_max) + 
      theme_minimal() +
      labs(title = ylabel,
           x = "",
           y = "") +
      theme(axis.text.x = element_text(size = x_size, angle = 45, hjust = 1),
            axis.title.y = element_text(size = y_size),
            plot.title = element_text(size = title_size, face = "bold")) +
      scale_color_manual(values = c("#FFA600", "#003F5C"),
                         labels = c("+", "-")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_jitter(aes(color = Group_Type), size = jitter_size, width = 0.2, alpha = geom_alpha)   # Add jitter points
    
  }
  
  return(p)
}

plot_sim5_global = function(df,  ylabel = "Dist", 
                           y_val_sig_bars = NA, 
                           x_size = 20,
                           y_size = 18,
                           thresh = 0.05,
                           jitter_size = 0.2, 
                           geom_alpha = 0.6, 
                           title_size = 18){
  
  df <- 
    df %>%
    pivot_longer(cols = everything(), names_to = "Simulation", values_to = "Dist") %>%
    mutate(Group_Type = ifelse(stringr::str_detect(Simulation, '_guides'), "With Guides", "Without Guides")) %>%
    mutate(Group = case_when(Simulation == "V1V2" ~ "V1V2-",
                             Simulation == "V1V2_guides" ~ "V1V2+",
                             Simulation == "V1V3" ~ "V1V3-",
                             Simulation == "V1V3_guides" ~ "V1V3+",
                             
                             Simulation == "V3V4" ~ "V3V4-",
                             Simulation == "V3V4_guides" ~ "V3V4+",
                             Simulation == "V3V5" ~ "V3V5-",
                             Simulation == "V3V5_guides" ~ "V3V5+",
                             Simulation == "V4" ~ "V4-",
                             Simulation == "V4_guides" ~ "V4+",
                             Simulation == "V4V5" ~ "V4V5-",
                             Simulation == "V4V5_guides" ~ "V4V5+",
                             
                             Simulation == "nd_mix" ~ "Mix-",
                             Simulation == "nd_mix_guides" ~ "Mix+",
                             TRUE ~ Simulation)) %>%
    select(!Simulation) %>%
    rename("Value" = Dist)
  df
  
  # Order by median value
  df <- df %>%
    group_by(Group) %>%
    mutate(median_value = median(Value)) %>%
    arrange(Group_Type, desc(median_value))
  df
  
  # Perform pairwise Wilcoxon Rank Sum tests within each Group_Type
  pairwise_tests = pairwise_wilcox_test(data = data.frame(df), 
                                        formula = Value ~ Group, 
                                        p.adjust.method = "fdr", 
                                        paired = FALSE)
  print(pairwise_tests)
  
  
  if(pairwise_tests %>% filter(p.adj > thresh) %>% nrow() > 0){
    pairwise_tests= pairwise_tests %>% filter(p.adj > thresh)
    pairwise_tests$y.position = y_val_sig_bars
    p = ggplot(df, aes(x = reorder(Group, -median_value), y = Value)) +
      geom_boxplot(aes(color = Group_Type)) +
      theme_minimal() +
      labs(title = ylabel,
           x = "",
           y = "") +
      theme(axis.text.x = element_text(size = x_size, angle = 45, hjust = 1),
            axis.title.y = element_text(size = y_size),
            plot.title = element_text(size = title_size, face = "bold")) +
      scale_color_manual(values = c("#FFA600", "#003F5C"),
                         labels = c("+", "-")) +
      geom_jitter(aes(color = Group_Type), size = jitter_size, width = 0.2, alpha = geom_alpha) +  # Add jitter points
      stat_pvalue_manual(pairwise_tests, label = "p.adj.signif", tip.length = 0.02) 
    
  }else{
    p = ggplot(df, aes(x = reorder(Group, -median_value), y = Value)) +
      geom_boxplot(aes(color = Group_Type)) +
      theme_minimal() +
      labs(title = ylabel,
           x = "",
           y = "") +
      theme(axis.text.x = element_text(size = x_size, angle = 45, hjust = 1),
            axis.title.y = element_text(size = y_size),
            plot.title = element_text(size = title_size, face = "bold")) +
      geom_jitter(aes(color = Group_Type), size = jitter_size, width = 0.2, alpha = geom_alpha) +  # Add jitter point
      scale_color_manual(values = c("#FFA600", "#003F5C"),
                         labels = c("+", "-")) 
  }
  
  return(p)
}

plot_sim4_local_increasing = function(df,  ylabel = "Dist", 
                           y_val_sig_bars = rep(1, times = 13), 
                           thresh = 0.05,
                           x_size = 20,
                           y_size = 18,
                           y_min = 0,
                           y_max = 1, 
                           jitter_size = 0.2, 
                           geom_alpha = 0.6, 
                           title_size = 18){
  
  df <- 
    df %>%
    pivot_longer(cols = everything(), names_to = "Simulation", values_to = "Dist") %>%
    mutate(Group_Type = ifelse(stringr::str_detect(Simulation, '_guides'), "With Guides", "Without Guides")) %>%
    mutate(Group = case_when(Simulation == "V1V2" ~ "V1V2-",
                             Simulation == "V1V2_guides" ~ "V1V2+",
                             Simulation == "V3V4" ~ "V3V4-",
                             Simulation == "V3V4_guides" ~ "V3V4+",
                             Simulation == "V5V7" ~ "V5V7-",
                             Simulation == "V5V7_guides" ~ "V5V7+",
                             Simulation == "V1V2_V3V4_V5V7_mix" ~ "Mix-",
                             Simulation == "V1V2_V3V4_V5V7_mix_guides" ~ "Mix+",
                             TRUE ~ Simulation)) %>%
    select(!Simulation) %>%
    rename("Value" = Dist)
  df
  
  # Order by median value
  df <- df %>%
    group_by(Group) %>%
    mutate(median_value = median(Value)) %>%
    arrange(Group_Type, desc(median_value))
  df
  
  # Perform pairwise Wilcoxon Rank Sum tests within each Group_Type
  pairwise_tests = pairwise_wilcox_test(data = data.frame(df), 
                                        formula = Value ~ Group, 
                                        p.adjust.method = "fdr", 
                                        paired = FALSE)
  as_tibble(pairwise_tests) %>% print(n = 100)
  
  
  
  if(pairwise_tests %>% filter(p.adj > thresh) %>% nrow() > 0){
    pairwise_tests= pairwise_tests %>% filter(p.adj > thresh)
    pairwise_tests$y.position = y_val_sig_bars
    p = ggplot(df, aes(x = reorder(Group, median_value), y = Value)) +
      geom_boxplot(aes(color = Group_Type)) +
      ylim(y_min, y_max) + 
      theme_minimal() +
      labs(title = ylabel,
           x = "",
           y = "") +
      theme(axis.text.x = element_text(size = x_size, angle = 45, hjust = 1),
            axis.title.y = element_text(size = y_size),
            plot.title = element_text(size = title_size, face = "bold")) +
      scale_color_manual(values = c("#FFA600", "#003F5C"),
                         labels = c("+", "-")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      geom_jitter(aes(color = Group_Type), size = jitter_size, width = 0.2, alpha = geom_alpha) +  # Add jitter points
      stat_pvalue_manual(pairwise_tests, label = "p.adj.signif", tip.length = 0.02) 
  }else{
    p = ggplot(df, aes(x = reorder(Group, median_value), y = Value)) +
      geom_boxplot(aes(color = Group_Type)) +
      ylim(y_min, y_max) + 
      theme_minimal() +
      labs(title = ylabel,
           x = "",
           y = "") +
      theme(axis.text.x = element_text(size = x_size, angle = 45, hjust = 1),
            axis.title.y = element_text(size = y_size),
            plot.title = element_text(size = title_size, face = "bold")) +
      scale_color_manual(values = c("#FFA600", "#003F5C"),
                         labels = c("+", "-")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_jitter(aes(color = Group_Type), size = jitter_size, width = 0.2, alpha = geom_alpha)   # Add jitter points
    
  }
  
  return(p)
}

plot_sim4_local = function(df,  ylabel = "Dist", 
                           y_val_sig_bars = rep(1, times = 13), 
                           thresh = 0.05,
                           x_size = 20,
                           y_size = 18,
                           y_min = 0,
                           y_max = 1, 
                           jitter_size = 0.2, 
                           geom_alpha = 0.6, 
                           title_size = 18){
  
  df <- 
    df %>%
    pivot_longer(cols = everything(), names_to = "Simulation", values_to = "Dist") %>%
    mutate(Group_Type = ifelse(stringr::str_detect(Simulation, '_guides'), "With Guides", "Without Guides")) %>%
    mutate(Group = case_when(Simulation == "V1V2" ~ "V1V2-",
                             Simulation == "V1V2_guides" ~ "V1V2+",
                             Simulation == "V3V4" ~ "V3V4-",
                             Simulation == "V3V4_guides" ~ "V3V4+",
                             Simulation == "V5V7" ~ "V5V7-",
                             Simulation == "V5V7_guides" ~ "V5V7+",
                             Simulation == "V1V2_V3V4_V5V7_mix" ~ "Mix-",
                             Simulation == "V1V2_V3V4_V5V7_mix_guides" ~ "Mix+",
                             TRUE ~ Simulation)) %>%
    select(!Simulation) %>%
    rename("Value" = Dist)
  df
  
  # Order by median value
  df <- df %>%
    group_by(Group) %>%
    mutate(median_value = median(Value)) %>%
    arrange(Group_Type, desc(median_value))
  df
  
  # Perform pairwise Wilcoxon Rank Sum tests within each Group_Type
  pairwise_tests = pairwise_wilcox_test(data = data.frame(df), 
                                        formula = Value ~ Group, 
                                        p.adjust.method = "fdr", 
                                        paired = FALSE)
  as_tibble(pairwise_tests) %>% print(n = 100)
  
  
  
  if(pairwise_tests %>% filter(p.adj > thresh) %>% nrow() > 0){
    pairwise_tests= pairwise_tests %>% filter(p.adj > thresh)
    pairwise_tests$y.position = y_val_sig_bars
    p = ggplot(df, aes(x = reorder(Group, -median_value), y = Value)) +
      geom_boxplot(aes(color = Group_Type)) +
      ylim(y_min, y_max) + 
      theme_minimal() +
      labs(title = ylabel,
           x = "",
           y = "") +
      theme(axis.text.x = element_text(size = x_size, angle = 45, hjust = 1),
            axis.title.y = element_text(size = y_size),
            plot.title = element_text(size = title_size, face = "bold")) +
      scale_color_manual(values = c("#FFA600", "#003F5C"),
                         labels = c("+", "-")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      geom_jitter(aes(color = Group_Type), size = jitter_size, width = 0.2, alpha = geom_alpha) +  # Add jitter points
      stat_pvalue_manual(pairwise_tests, label = "p.adj.signif", tip.length = 0.02) 
  }else{
    p = ggplot(df, aes(x = reorder(Group, -median_value), y = Value)) +
      geom_boxplot(aes(color = Group_Type)) +
      ylim(y_min, y_max) + 
      theme_minimal() +
      labs(title = ylabel,
           x = "",
           y = "") +
      theme(axis.text.x = element_text(size = x_size, angle = 45, hjust = 1),
            axis.title.y = element_text(size = y_size),
            plot.title = element_text(size = title_size, face = "bold")) +
      scale_color_manual(values = c("#FFA600", "#003F5C"),
                         labels = c("+", "-")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_jitter(aes(color = Group_Type), size = jitter_size, width = 0.2, alpha = geom_alpha)   # Add jitter points
      
  }
  
  return(p)
}


plot_sim4_global = function(df,  ylabel = "Dist", 
                            y_val_sig_bars = rep(1, times = 13), 
                            thresh = 0.05,
                            x_size = 20,
                            y_size = 18, 
                            jitter_size = 0.2, 
                            geom_alpha = 0.6, 
                            title_size = 18){
  
  df <- 
    df %>%
    pivot_longer(cols = everything(), names_to = "Simulation", values_to = "Dist") %>%
    mutate(Group_Type = ifelse(stringr::str_detect(Simulation, '_guides'), "With Guides", "Without Guides")) %>%
    mutate(Group = case_when(Simulation == "V1V2" ~ "V1V2-",
                             Simulation == "V1V2_guides" ~ "V1V2+",
                             Simulation == "V3V4" ~ "V3V4-",
                             Simulation == "V3V4_guides" ~ "V3V4+",
                             Simulation == "V5V7" ~ "V5V7-",
                             Simulation == "V5V7_guides" ~ "V5V7+",
                             Simulation == "V1V2_V3V4_V5V7_mix" ~ "Mix-",
                             Simulation == "V1V2_V3V4_V5V7_mix_guides" ~ "Mix+",
                             TRUE ~ Simulation)) %>%
    select(!Simulation) %>%
    rename("Value" = Dist)
  df
  
  # Order by median value
  df <- df %>%
    group_by(Group) %>%
    mutate(median_value = median(Value)) %>%
    arrange(Group_Type, desc(median_value))
  df
  
  # Perform pairwise Wilcoxon Rank Sum tests within each Group_Type
  pairwise_tests = pairwise_wilcox_test(data = data.frame(df), 
                                        formula = Value ~ Group, 
                                        p.adjust.method = "fdr", 
                                        paired = FALSE)
  as_tibble(pairwise_tests) %>% print(n = 100)
  
  
  
  if(pairwise_tests %>% filter(p.adj > thresh) %>% nrow() > 0){
    pairwise_tests= pairwise_tests %>% filter(p.adj > thresh)
    pairwise_tests$y.position = y_val_sig_bars
    p = ggplot(df, aes(x = reorder(Group, -median_value), y = Value)) +
      geom_boxplot(aes(color = Group_Type)) +
      
      theme_minimal() +
      labs(title = ylabel,
           x = "",
           y = "") +
      theme(axis.text.x = element_text(size = x_size, angle = 45, hjust = 1),
            axis.title.y = element_text(size = y_size),
            plot.title = element_text(size = title_size, face = "bold")) +
      scale_color_manual(values = c("#FFA600", "#003F5C"),
                         labels = c("+", "-")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      geom_jitter(aes(color = Group_Type), size = jitter_size, width = 0.2, alpha = geom_alpha) +  # Add jitter points
      stat_pvalue_manual(pairwise_tests, label = "p.adj.signif", tip.length = 0.02) 
  }else{
    p = ggplot(df, aes(x = reorder(Group, -median_value), y = Value)) +
      geom_boxplot(aes(color = Group_Type)) +
      
      theme_minimal() +
      labs(title = ylabel,
           x = "",
           y = "") +
      theme(axis.text.x = element_text(size = x_size, angle = 45, hjust = 1),
            axis.title.y = element_text(size = y_size),
            plot.title = element_text(size = title_size, face = "bold")) +
      scale_color_manual(values = c("#FFA600", "#003F5C"),
                         labels = c("+", "-")) +
      geom_jitter(aes(color = Group_Type), size = jitter_size, width = 0.2, alpha = geom_alpha) +  # Add jitter points
    
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  }
  
  return(p)
}
plot_sim4_global_increasing = function(df,  ylabel = "Dist", 
                                       y_val_sig_bars = rep(1, times = 13), 
                                       thresh = 0.05,
                                       x_size = 20,
                                       y_size = 18, 
                                       jitter_size = 0.2, 
                                       geom_alpha = 0.6, 
                                       title_size = 18){
  
  df <- 
    df %>%
    pivot_longer(cols = everything(), names_to = "Simulation", values_to = "Dist") %>%
    mutate(Group_Type = ifelse(stringr::str_detect(Simulation, '_guides'), "With Guides", "Without Guides")) %>%
    mutate(Group = case_when(Simulation == "V1V2" ~ "V1V2-",
                             Simulation == "V1V2_guides" ~ "V1V2+",
                             Simulation == "V3V4" ~ "V3V4-",
                             Simulation == "V3V4_guides" ~ "V3V4+",
                             Simulation == "V5V7" ~ "V5V7-",
                             Simulation == "V5V7_guides" ~ "V5V7+",
                             Simulation == "V1V2_V3V4_V5V7_mix" ~ "Mix-",
                             Simulation == "V1V2_V3V4_V5V7_mix_guides" ~ "Mix+",
                             TRUE ~ Simulation)) %>%
    select(!Simulation) %>%
    rename("Value" = Dist)
  df
  
  # Order by median value
  df <- df %>%
    group_by(Group) %>%
    mutate(median_value = median(Value)) %>%
    arrange(Group_Type, desc(median_value))
  df
  
  # Perform pairwise Wilcoxon Rank Sum tests within each Group_Type
  pairwise_tests = pairwise_wilcox_test(data = data.frame(df), 
                                        formula = Value ~ Group, 
                                        p.adjust.method = "fdr", 
                                        paired = FALSE)
  as_tibble(pairwise_tests) %>% print(n = 100)
  
  
  
  if(pairwise_tests %>% filter(p.adj > thresh) %>% nrow() > 0){
    pairwise_tests= pairwise_tests %>% filter(p.adj > thresh)
    pairwise_tests$y.position = y_val_sig_bars
    p = ggplot(df, aes(x = reorder(Group, median_value), y = Value)) +
      geom_boxplot(aes(color = Group_Type)) +
      
      theme_minimal() +
      labs(title = ylabel,
           x = "",
           y = "") +
      theme(axis.text.x = element_text(size = x_size, angle = 45, hjust = 1),
            axis.title.y = element_text(size = y_size),
            plot.title = element_text(size = title_size, face = "bold")) +
      scale_color_manual(values = c("#FFA600", "#003F5C"),
                         labels = c("+", "-")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      geom_jitter(aes(color = Group_Type), size = jitter_size, width = 0.2, alpha = geom_alpha) +  # Add jitter points
      stat_pvalue_manual(pairwise_tests, label = "p.adj.signif", tip.length = 0.02) 
  }else{
    p = ggplot(df, aes(x = reorder(Group, median_value), y = Value)) +
      geom_boxplot(aes(color = Group_Type)) +
      
      theme_minimal() +
      labs(title = ylabel,
           x = "",
           y = "") +
      theme(axis.text.x = element_text(size = x_size, angle = 45, hjust = 1),
            axis.title.y = element_text(size = y_size),
            plot.title = element_text(size = title_size, face = "bold")) +
      scale_color_manual(values = c("#FFA600", "#003F5C"),
                         labels = c("+", "-")) +
      geom_jitter(aes(color = Group_Type), size = jitter_size, width = 0.2, alpha = geom_alpha) +  # Add jitter points
      
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  }
  
  return(p)
}



plot_sim3_local = function(df,  ylabel = "Dist", 
                           y_val_sig_bars = NA, 
                           x_size = 20,
                           y_size = 18,
                           thresh = 0.05,
                           y_min = 0,
                           y_max = 1, 
                           jitter_size = 0.2, 
                           geom_alpha = 0.6, 
                           title_size = 18){
  
  df <- 
    df %>%
    pivot_longer(cols = everything(), names_to = "Simulation", values_to = "Dist") %>%
    mutate(Group_Type = ifelse(stringr::str_detect(Simulation, '_guides'), "With Guides", "Without Guides")) %>%
    mutate(Group = case_when(Simulation == "V4" ~ "V4-",
                             Simulation == "V4_guides" ~ "V4+",
                             Simulation == "V3V4" ~ "V3V4-",
                             Simulation == "V3V4_guides" ~ "V3V4+",
                             Simulation == "V4_V3V4_mix" ~ "Mix-",
                             Simulation == "V4_V3V4_mix_guides" ~ "Mix+",
                             TRUE ~ Simulation)) %>%
    select(!Simulation) %>%
    rename("Value" = Dist)
  df
  
  # Order by median value
  df <- df %>%
    group_by(Group) %>%
    mutate(median_value = median(Value)) %>%
    arrange(Group_Type, desc(median_value))
  df
  
  # Perform pairwise Wilcoxon Rank Sum tests within each Group_Type
  pairwise_tests = pairwise_wilcox_test(data = data.frame(df), 
                                        formula = Value ~ Group, 
                                        p.adjust.method = "fdr", 
                                        paired = FALSE)
  print(pairwise_tests)
  
  
  if(pairwise_tests %>% filter(p.adj > thresh) %>% nrow() > 0){
    pairwise_tests= pairwise_tests %>% filter(p.adj > thresh)
    pairwise_tests$y.position = y_val_sig_bars
    p = ggplot(df, aes(x = reorder(Group, -median_value), y = Value)) +
      geom_boxplot(aes(color = Group_Type)) +
      theme_minimal() +
      ylim(y_min, y_max) + 
      labs(title = ylabel,
           x = "",
           y = "") +
      theme(axis.text.x = element_text(size = x_size, angle = 45, hjust = 1),
            axis.title.y = element_text(size = y_size),
            plot.title = element_text(size = title_size, face = "bold")) +
      scale_color_manual(values = c("#FFA600", "#003F5C"),
                         labels = c("+", "-")) +
      geom_jitter(aes(color = Group_Type), size = jitter_size, width = 0.2, alpha = geom_alpha) +  # Add jitter points
      stat_pvalue_manual(pairwise_tests, label = "p.adj.signif", tip.length = 0.02) 
    
  }else{
    p = ggplot(df, aes(x = reorder(Group, -median_value), y = Value)) +
      geom_boxplot(aes(color = Group_Type)) +
      theme_minimal() +
      ylim(y_min, y_max) + 
      labs(title = ylabel,
           x = "",
           y = "") +
      theme(axis.text.x = element_text(size = x_size, angle = 45, hjust = 1),
            axis.title.y = element_text(size = y_size),
            plot.title = element_text(size = title_size, face = "bold")) +
      geom_jitter(aes(color = Group_Type), size = jitter_size, width = 0.2, alpha = geom_alpha) +  # Add jitter point
      scale_color_manual(values = c("#FFA600", "#003F5C"),
                         labels = c("+", "-")) 
  }
  
  return(p)
}
plot_sim3_local_increasing = function(df,  ylabel = "Dist", 
                           y_val_sig_bars = NA, 
                           x_size = 20,
                           y_size = 18,
                           thresh = 0.05,
                           y_min = 0,
                           y_max = 1, 
                           jitter_size = 0.2, 
                           geom_alpha = 0.6, 
                           title_size = 18){
  
  df <- 
    df %>%
    pivot_longer(cols = everything(), names_to = "Simulation", values_to = "Dist") %>%
    mutate(Group_Type = ifelse(stringr::str_detect(Simulation, '_guides'), "With Guides", "Without Guides")) %>%
    mutate(Group = case_when(Simulation == "V4" ~ "V4-",
                             Simulation == "V4_guides" ~ "V4+",
                             Simulation == "V3V4" ~ "V3V4-",
                             Simulation == "V3V4_guides" ~ "V3V4+",
                             Simulation == "V4_V3V4_mix" ~ "Mix-",
                             Simulation == "V4_V3V4_mix_guides" ~ "Mix+",
                             TRUE ~ Simulation)) %>%
    select(!Simulation) %>%
    rename("Value" = Dist)
  df
  
  # Order by median value
  df <- df %>%
    group_by(Group) %>%
    mutate(median_value = median(Value)) %>%
    arrange(Group_Type, desc(median_value))
  df
  
  # Perform pairwise Wilcoxon Rank Sum tests within each Group_Type
  pairwise_tests = pairwise_wilcox_test(data = data.frame(df), 
                                        formula = Value ~ Group, 
                                        p.adjust.method = "fdr", 
                                        paired = FALSE)
  print(pairwise_tests)
  
  
  if(pairwise_tests %>% filter(p.adj > thresh) %>% nrow() > 0){
    pairwise_tests= pairwise_tests %>% filter(p.adj > thresh)
    pairwise_tests$y.position = y_val_sig_bars
    p = ggplot(df, aes(x = reorder(Group, median_value), y = Value)) +
      geom_boxplot(aes(color = Group_Type)) +
      theme_minimal() +
      ylim(y_min, y_max) + 
      labs(title = ylabel,
           x = "",
           y = "") +
      theme(axis.text.x = element_text(size = x_size, angle = 45, hjust = 1),
            axis.title.y = element_text(size = y_size),
            plot.title = element_text(size = title_size, face = "bold")) +
      scale_color_manual(values = c("#FFA600", "#003F5C"),
                         labels = c("+", "-")) +
      geom_jitter(aes(color = Group_Type), size = jitter_size, width = 0.2, alpha = geom_alpha) +  # Add jitter points
      stat_pvalue_manual(pairwise_tests, label = "p.adj.signif", tip.length = 0.02) 
    
  }else{
    p = ggplot(df, aes(x = reorder(Group, median_value), y = Value)) +
      geom_boxplot(aes(color = Group_Type)) +
      theme_minimal() +
      ylim(y_min, y_max) + 
      labs(title = ylabel,
           x = "",
           y = "") +
      theme(axis.text.x = element_text(size = x_size, angle = 45, hjust = 1),
            axis.title.y = element_text(size = y_size),
            plot.title = element_text(size = title_size, face = "bold")) +
      geom_jitter(aes(color = Group_Type), size = jitter_size, width = 0.2, alpha = geom_alpha) +  # Add jitter point
      scale_color_manual(values = c("#FFA600", "#003F5C"),
                         labels = c("+", "-")) 
  }
  
  return(p)
}

plot_sim3_global = function(df,  ylabel = "Dist", 
                            y_val_sig_bars = NA, 
                            x_size = 20,
                            y_size = 18,
                            thresh = 0.05, 
                            jitter_size = 0.2, 
                            geom_alpha = 0.6, 
                            title_size = 18){
  
  df <- 
    df %>%
    pivot_longer(cols = everything(), names_to = "Simulation", values_to = "Dist") %>%
    mutate(Group_Type = ifelse(stringr::str_detect(Simulation, '_guides'), "With Guides", "Without Guides")) %>%
    mutate(Group = case_when(Simulation == "V4" ~ "V4-",
                             Simulation == "V4_guides" ~ "V4+",
                             Simulation == "V3V4" ~ "V3V4-",
                             Simulation == "V3V4_guides" ~ "V3V4+",
                             Simulation == "V4_V3V4_mix" ~ "Mix-",
                             Simulation == "V4_V3V4_mix_guides" ~ "Mix+",
                             TRUE ~ Simulation)) %>%
    select(!Simulation) %>%
    rename("Value" = Dist)
  df
  
  # Order by median value
  df <- df %>%
    group_by(Group) %>%
    mutate(median_value = median(Value)) %>%
    arrange(Group_Type, desc(median_value))
  df
  
  # Perform pairwise Wilcoxon Rank Sum tests within each Group_Type
  pairwise_tests = pairwise_wilcox_test(data = data.frame(df), 
                                        formula = Value ~ Group, 
                                        p.adjust.method = "fdr", 
                                        paired = FALSE)
  print(pairwise_tests)
  
  
  if(pairwise_tests %>% filter(p.adj > thresh) %>% nrow() > 0){
    pairwise_tests= pairwise_tests %>% filter(p.adj > thresh)
    pairwise_tests$y.position = y_val_sig_bars
    
    p = ggplot(df, aes(x = reorder(Group, -median_value), y = Value)) +
      geom_boxplot(aes(color = Group_Type)) +
      theme_minimal() +
      labs(title = ylabel,
           x = "",
           y = "") +
      theme(axis.text.x = element_text(size = x_size, angle = 45, hjust = 1),
            axis.title.y = element_text(size = y_size),
            plot.title = element_text(size = title_size, face = "bold")) +
      scale_color_manual(values = c("#FFA600", "#003F5C"),
                         labels = c("+", "-")) +
      geom_jitter(aes(color = Group_Type), size = jitter_size, width = 0.2, alpha = geom_alpha) +  # Add jitter points
      stat_pvalue_manual(pairwise_tests, label = "p.adj.signif", tip.length = 0.02) 
    
  }else{
    p = ggplot(df, aes(x = reorder(Group, -median_value), y = Value)) +
      geom_boxplot(aes(color = Group_Type)) +
      theme_minimal() +
      labs(title = ylabel,
           x = "",
           y = "") +
      theme(axis.text.x = element_text(size = x_size, angle = 45, hjust = 1),
            axis.title.y = element_text(size = y_size),
            plot.title = element_text(size = title_size, face = "bold")) +
      scale_color_manual(values = c("#FFA600", "#003F5C"),
                         labels = c("+", "-")) + 
      geom_jitter(aes(color = Group_Type), size = jitter_size, width = 0.2, alpha = geom_alpha)   # Add jitter points
  }
  
  return(p)
}

plot_sim3_global_increasing = function(df,  ylabel = "Dist", 
                            y_val_sig_bars = NA, 
                            x_size = 20,
                            y_size = 18,
                            thresh = 0.05, 
                            jitter_size = 0.2, 
                            geom_alpha = 0.6, 
                            title_size = 18){
  
  df <- 
    df %>%
    pivot_longer(cols = everything(), names_to = "Simulation", values_to = "Dist") %>%
    mutate(Group_Type = ifelse(stringr::str_detect(Simulation, '_guides'), "With Guides", "Without Guides")) %>%
    mutate(Group = case_when(Simulation == "V4" ~ "V4-",
                             Simulation == "V4_guides" ~ "V4+",
                             Simulation == "V3V4" ~ "V3V4-",
                             Simulation == "V3V4_guides" ~ "V3V4+",
                             Simulation == "V4_V3V4_mix" ~ "Mix-",
                             Simulation == "V4_V3V4_mix_guides" ~ "Mix+",
                             TRUE ~ Simulation)) %>%
    select(!Simulation) %>%
    rename("Value" = Dist)
  df
  
  # Order by median value
  df <- df %>%
    group_by(Group) %>%
    mutate(median_value = median(Value)) %>%
    arrange(Group_Type, desc(median_value))
  df
  
  # Perform pairwise Wilcoxon Rank Sum tests within each Group_Type
  pairwise_tests = pairwise_wilcox_test(data = data.frame(df), 
                                        formula = Value ~ Group, 
                                        p.adjust.method = "fdr", 
                                        paired = FALSE)
  print(pairwise_tests)
  
  
  if(pairwise_tests %>% filter(p.adj > thresh) %>% nrow() > 0){
    pairwise_tests= pairwise_tests %>% filter(p.adj > thresh)
    pairwise_tests$y.position = y_val_sig_bars
    
    p = ggplot(df, aes(x = reorder(Group, median_value), y = Value)) +
      geom_boxplot(aes(color = Group_Type)) +
      theme_minimal() +
      labs(title = ylabel,
           x = "",
           y = "") +
      theme(axis.text.x = element_text(size = x_size, angle = 45, hjust = 1),
            axis.title.y = element_text(size = y_size),
            plot.title = element_text(size = title_size, face = "bold")) +
      scale_color_manual(values = c("#FFA600", "#003F5C"),
                         labels = c("+", "-")) +
      geom_jitter(aes(color = Group_Type), size = jitter_size, width = 0.2, alpha = geom_alpha) +  # Add jitter points
      stat_pvalue_manual(pairwise_tests, label = "p.adj.signif", tip.length = 0.02) 
    
  }else{
    p = ggplot(df, aes(x = reorder(Group, median_value), y = Value)) +
      geom_boxplot(aes(color = Group_Type)) +
      theme_minimal() +
      labs(title = ylabel,
           x = "",
           y = "") +
      theme(axis.text.x = element_text(size = x_size, angle = 45, hjust = 1),
            axis.title.y = element_text(size = y_size),
            plot.title = element_text(size = title_size, face = "bold")) +
      scale_color_manual(values = c("#FFA600", "#003F5C"),
                         labels = c("+", "-")) + 
      geom_jitter(aes(color = Group_Type), size = jitter_size, width = 0.2, alpha = geom_alpha)   # Add jitter points
  }
  
  return(p)
}



# For an interative mast list, sum all of the numbers in each list within the list
micro_total = function(list){
  totals = vector()
  for(i in 1:length(list)){
    totals = c(totals, sum(list[[i]]))
  }
  return(totals)
}
# For an interative mast list, calculate mean for each list within the list

micro_mean_size = function(list){
  totals = vector()
  for(i in 1:length(list)){
    totals = c(totals, mean(log(list[[i]])))
  }
  return(totals)
}
# Code to make figures hvr_order_and_significance_comparisons_A.pdf and hvr_order_and_significance_comparisons_B.pdf
get_boxplot_global_dist_ordered = function(df, ylabel){
  # Rename columns in hvr_rf: Change "_guides" to "_+"
  hvr_rf_renamed <- df %>%
    select(!"VFull_guides") %>%
    rename("C" = VFull) %>%
    rename_with(~ stringr::str_replace(., "_guides$", "+"))
  
  # Convert to long format
  df_long <- hvr_rf_renamed %>%
    pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
    mutate(
      Guide_Status = ifelse(grepl("\\+$", Variable), "Guides", "Non-Guides")  # Match "_+" for guides
    )
  
  # Reorder variables by median value
  df_long <- df_long %>%
    group_by(Variable) %>%
    mutate(Median_Value = median(Value, na.rm = TRUE)) %>%
    ungroup()
  
  # Create the boxplot
  p_rf_ordered = ggplot(df_long, aes(x = reorder(Variable, Median_Value), y = Value, color = Guide_Status)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # Boxplot with semi-transparent fill
    scale_color_manual(values = c("Non-Guides" = "#003F5C", "Guides" = "#FFA600")) + # Custom colors
    theme_minimal() +
    labs(x = "", y = ylabel, color = "Category") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  return(p_rf_ordered)
}

# Code to make figures hvr_order_and_significance_comparisons_A.pdf and hvr_order_and_significance_comparisons_B.pdf
get_heatmap_pairwise_global_dist = function(df, q_threshold = 0.01, q_threshold_label = "q<=0.01", correction_method = "fdr", title){
  
  hvr_rf_renamed <- df %>%
    select(!"VFull_guides") %>%
    rename("C" = VFull) %>%
    rename_with(~ stringr::str_replace(., "_guides$", "+"))
  
  # Generate all unique column pairs
  column_names <- colnames(hvr_rf_renamed)
  
  # Perform pairwise tests
  pairwise_tests <- expand.grid(column_names, column_names) %>%
    filter(Var1 != Var2) %>%
    rowwise() %>%
    mutate(p_value = wilcox.test(hvr_rf_renamed[[Var1]], hvr_rf_renamed[[Var2]], paired = FALSE)$p.value) %>%
    ungroup()
  
  # Adjust p-values for multiple testing (optional)
  pairwise_tests <- pairwise_tests %>%
    mutate(p_adj = p.adjust(p_value, method = correction_method))  # False Discovery Rate correction
  
  # Convert to matrix format for heatmap
  heatmap_matrix <- dcast(pairwise_tests, Var1 ~ Var2, value.var = "p_adj")
  
  # Convert to long format for ggplot2
  heatmap_long <- melt(heatmap_matrix, id.vars = "Var1") %>% 
    mutate(q = case_when(value <= q_threshold ~ "*",
                         value > q_threshold ~ "ns",
                         TRUE ~ NA_character_))
  
  p_rf_ordered_p = ggplot(heatmap_long, aes(x = variable, y = Var1, fill = q)) +
    geom_tile() +
    scale_fill_manual(values = c("*" = "#B30000", "ns" = "#253494"), na.value = "#F0F0F0") +  
    theme_minimal() +
    labs(x = "", y = "", fill = q_threshold_label, title = title) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  return(p_rf_ordered_p)
}

pearson_cor_hvr = function(dist1, dist2){
  d = data.frame(matrix(nrow = 20, ncol = 1))
  colnames(d) = c("Pearson")
  rownames(d) = c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "VFull", 
                  "V1_guides", "V2_guides", "V3_guides", "V4_guides", "V5_guides", "V6_guides", 
                  "V7_guides", "V8_guides", "V9_guides", "VFull_guides")
  
  d["V1", "Pearson"] = cor(x = dist1$V1, y = dist2$V1, method = "pearson")
  d["V2", "Pearson"] = cor(x = dist1$V2, y = dist2$V2, method = "pearson")
  d["V3", "Pearson"] = cor(x = dist1$V3, y = dist2$V3, method = "pearson")
  d["V4", "Pearson"] = cor(x = dist1$V4, y = dist2$V4, method = "pearson")
  d["V5", "Pearson"] = cor(x = dist1$V5, y = dist2$V5, method = "pearson")
  d["V6", "Pearson"] = cor(x = dist1$V6, y = dist2$V6, method = "pearson")
  d["V7", "Pearson"] = cor(x = dist1$V7, y = dist2$V7, method = "pearson")
  d["V8", "Pearson"] = cor(x = dist1$V8, y = dist2$V8, method = "pearson")
  d["V9", "Pearson"] = cor(x = dist1$V9, y = dist2$V9, method = "pearson")
  d["VFull", "Pearson"] = cor(x = dist1$VFull, y = dist2$VFull, method = "pearson")
  
  d["V1_guides", "Pearson"] = cor(x = dist1$V1_guides, y = dist2$V1_guides, method = "pearson")
  d["V2_guides", "Pearson"] = cor(x = dist1$V2_guides, y = dist2$V2_guides, method = "pearson")
  d["V3_guides", "Pearson"] = cor(x = dist1$V3_guides, y = dist2$V3_guides, method = "pearson")
  d["V4_guides", "Pearson"] = cor(x = dist1$V4_guides, y = dist2$V4_guides, method = "pearson")
  d["V5_guides", "Pearson"] = cor(x = dist1$V5_guides, y = dist2$V5_guides, method = "pearson")
  d["V6_guides", "Pearson"] = cor(x = dist1$V6_guides, y = dist2$V6_guides, method = "pearson")
  d["V7_guides", "Pearson"] = cor(x = dist1$V7_guides, y = dist2$V7_guides, method = "pearson")
  d["V8_guides", "Pearson"] = cor(x = dist1$V8_guides, y = dist2$V8_guides, method = "pearson")
  d["V9_guides", "Pearson"] = cor(x = dist1$V9_guides, y = dist2$V9_guides, method = "pearson")
  d["VFull_guides", "Pearson"] = cor(x = dist1$VFull_guides, y = dist2$VFull_guides, method = "pearson")
  
  return(d)
}

# Code to test every V1 - 9 with V1_guides - V9_guides with wilcox.test.
wilcox_test_hvr = function(dist){
  d = data.frame(matrix(nrow = 9, ncol = 2))
  colnames(d) = c("W", "p")
  rownames(d) = c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9")
  
  v1 = wilcox.test(dist$V1, dist$V1_guides, paired = FALSE)
  d["V1", "W"] = v1$statistic
  d["V1", "p"] = v1$p.value
  
  v2 = wilcox.test(dist$V2, dist$V2_guides, paired = FALSE)
  d["V2", "W"] = v2$statistic
  d["V2", "p"] = v2$p.value
  
  v3 = wilcox.test(dist$V3, dist$V3_guides, paired = FALSE)
  d["V3", "W"] = v3$statistic
  d["V3", "p"] = v3$p.value
  
  v4 = wilcox.test(dist$V4, dist$V4_guides, paired = FALSE)
  d["V4", "W"] = v4$statistic
  d["V4", "p"] = v4$p.value
  
  v5 = wilcox.test(dist$V5, dist$V5_guides, paired = FALSE)
  d["V5", "W"] = v5$statistic
  d["V5", "p"] = v5$p.value
  
  v6 = wilcox.test(dist$V6, dist$V6_guides, paired = FALSE)
  d["V6", "W"] = v6$statistic
  d["V6", "p"] = v6$p.value
  
  v7 = wilcox.test(dist$V7, dist$V7_guides, paired = FALSE)
  d["V7", "W"] = v7$statistic
  d["V7", "p"] = v7$p.value
  
  v8 = wilcox.test(dist$V8, dist$V8_guides, paired = FALSE)
  d["V8", "W"] = v8$statistic
  d["V8", "p"] = v8$p.value
  
  v9 = wilcox.test(dist$V9, dist$V9_guides, paired = FALSE)
  d["V9", "W"] = v9$statistic
  d["V9", "p"] = v9$p.value
  
  return(d)
}