# Plot statistically significant time points
  # runs on all ERP components, over all participants

  # from uncorrected univariate significance
  # from cluster-sum permutation testing
  # from maximum statistics correction

  # input: t2 and HT2 results from all participants in folder 'perm_results' 
  # (created in Permutations_and_analysis.m)
  # input: onset estimates, each EPR component one file 'r_onsets_all_methods.xlsx'
  # (created in onset_estimation.R)
  # input: results from statistical results from each participants in folder 'r_averages'
  # (created in onset_estimation.R)

  # no further external dependencies

library(R.matlab)
library(tidyverse)
library(readxl)
library(ggplot2)

SUB <- as.character(1:40) #missing 1-9
alpha <- 0.05

excludeMap <- list(
  N170  = c("1", "5", "16"),
  MMN   = c("7"),
  N2pc  = c("7", "9", "10", "12", "28"),
  P3    = c("6", "9", "10", "30", "35", "40"),
  N400  = c("40"),
  LRP   = c("6", "30", "40"),
  ERN   = c("5", "6", "30", "40"))

component_channels <- list(
  "N170" = 26,
  "MMN" = 20,
  "N2pc" = 9, 
  "P3" = 13,
  "N400" = 14,
  "LRP" = 5, 
  "ERN" = 20)

component_electrodes <- list(
  "N170" = 'PO8',
  "MMN" = 'FCz',
  "N2pc" = 'PO7', 
  "P3" = 'Pz',
  "N400" = 'CPz',
  "LRP" = 'C3', 
  "ERN" = 'FCz')

for (comp in components) { 
  
# Specify component parameters 
  exclude <- excludeMap[[comp]]
  chanIdx <- component_channels[[comp]]
 spec_elec <- component_electrodes[[comp]]


#Define baseline period
if (comp == "LRP") {
  end_of_baseline = -600
} else if (comp == "ERN") {
  end_of_baseline = -200
} else {
  end_of_baseline = 0}

# Specify main and component directories
main_dir <- "E:/Preprocessed_ERP_Core/All_components_files"
comp_dir <- file.path(main_dir, comp)
results_dir <- file.path(comp_dir, "perm_results")
averages_dir <- file.path(comp_dir, "r_averages")
fig_dir <- file.path(comp_dir, 'r_figures_univ_significance', comp) 
dir.create(fig_dir, showWarnings = FALSE)

# Load onset estimation results
onsets_path <- file.path(comp_dir, "r_onsets_all_methods.xlsx")
df_onsets <- read_excel(onsets_path)

# Loop over subjects
for (sub in SUB) {
  
  if (sub %in% exclude) {
    cat(sprintf("Skipping subject %s for component %s\n", sub, comp))
    next
  }
  
  # sub <- '1' # for running single participants
  
  filepath <- file.path(results_dir, paste0(sub, "_results.mat"))
  
  res <- readMat(filepath)$res
  field_names <- unlist(attr(res, "dimnames")[[1]])
  names(res) <- field_names
  
  Xf <- as.numeric(res$Xf)
  postStimIdx <- Xf >= end_of_baseline 
  Xf_post <- Xf[postStimIdx]
  
  max_tval <- as.numeric(res$max.tval)
  max_tcrit <- as.numeric(res$max.tcrit)
  max_pval <- as.numeric(res$max.pval)
  max_abs_tval <- abs(max_tval)
  
  max_t2_bn <- as.numeric(res$maxt2.bn)
  max_t2_perm_bn <- as.data.frame(res$maxt2.perm.bn)
  
  spec_tval <- as.numeric(res$spec.tval)
  spec_tcrit <- as.numeric(res$spec.tcrit)
  spec_pval <- as.numeric(res$spec.pval)
  spec_abs_tval <- abs(spec_tval)
  
  spec_t2_bn <- as.numeric(res$spec.t2.bn)
  spec_t2_perm_bn <- as.data.frame(res$spec.t2.perm.bn)
  
  # Load data from statistical analysis
  averages_path <- file.path(averages_dir, paste0("average_", sub, ".xlsx"))
  average_res <- read_excel(averages_path)
  
  # Maximum statistic correction
  max_t2_max_thresh <- average_res$max_t2_max_thresh
  spec_t2_max_thresh <- average_res$spec_t2_max_thresh
  
  # Cluster-sum permutation testing
  max_t2_cs_thresh <- average_res$max_t2_cs_thresh
  max_t2_cs_test <- average_res$max_t2_cs_test
  
  
  # Plot 1: t-values from virtual maximum t²electrode
  df_p1 <- data.frame(
    Xf = Xf,
    tval = max_abs_tval,
    pval = max_pval,
    alpha = alpha,
    sig = max_pval < alpha,  # TRUE/FALSE for significance
    onset_val = df_onsets$maxt2_cpd[df_onsets$participant == sub])
  
  p1 <- ggplot(df_p1, aes(x = Xf)) +
    geom_line(aes(y = tval, color = "Absolute t-values")) +
    geom_line(aes(y = pval, color = "Two-tailed p-values")) +
    geom_hline(aes(yintercept = max_tcrit, color = "Two-tailed critical threshold")) +
    geom_hline(aes(yintercept = alpha, color = "alpha = 0.05")) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
    geom_point(
      data = subset(df_p1, sig),
      aes(y = 0, color = "Significant timepoints"),
      size = 0.5) +
    geom_vline(aes(xintercept = onset_val, color = "CPD on max t²"), linetype = "dashed") +
    scale_color_manual(
      values = c(
        "Absolute t-values" = "black",
        "Two-tailed p-values" = "red",
        "Two-tailed critical threshold" = "blue",
        "alpha = 0.05" = "black",
        "Significant timepoints" = "green",
        "CPD on max t²" = "black" ),
      breaks = c(
        "Absolute t-values",
        "Two-tailed p-values",
        "Two-tailed critical threshold",
        "alpha = 0.05",
        "Significant timepoints",
        "CPD on max t²")) +
    labs(
      title = paste0(comp, " | P", sub, " | Univariate Significance at max t² electrode"),
      x = "Time (ms)",
      y = "Absolute t-value",
      color = "",
      linetype = "") +
    annotate("text",
             x = Inf, y = Inf,
             label = paste("Onset:", df_p1$onset_val[1], "ms"),
             hjust = 1, vjust = 1.5) +
    theme_light()
  
  
  ggsave(
    filename = file.path(fig_dir, paste0(sub, "_fig_ttest_max.png")),
    plot = p1, width = 6, height = 4, dpi = 300
  )
  
  # Plot 2: t-values from specified electrode of interest
  df_p2 <- data.frame(
    Xf = Xf,
    tval = spec_abs_tval,
    pval = spec_pval,
    alpha = alpha,
    sig = spec_pval < alpha,
    onset_val = df_onsets$spec_t2_cpd[df_onsets$participant == sub])
  
  p2 <- ggplot(df_p2, aes(x = Xf)) +
    geom_line(aes(y = tval, color = "Absolute t-values")) +
    geom_line(aes(y = pval, color = "Two-tailed p-values")) +
    geom_hline(aes(yintercept = max_tcrit, color = "Two-tailed critical threshold")) +
    geom_hline(aes(yintercept = alpha, color = "alpha = 0.05"), linetype = "dotted") +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
    geom_point(
      data = subset(df_p2, sig),
      aes(y = 0, color = "Significant timepoints"),
      size = 0.5) +
    geom_vline(aes(xintercept = onset_val, color = "CPD on t²"), linetype = "dashed") +
    scale_color_manual(
      values = c(
        "Absolute t-values" = "black",
        "Two-tailed p-values" = "red",
        "Two-tailed critical threshold" = "blue",
        "alpha = 0.05" = "black",
        "Significant timepoints" = "green",
        "CPD on t²" = "black" ),
      breaks = c(
        "Absolute t-values",
        "Two-tailed critical threshold",
        "Two-tailed p-values",
        "alpha = 0.05",
        "Significant timepoints",
        "CPD on t²")) +
    labs(
      title = paste0(comp, " | P", sub, " | Univariate Significance at Electrode ", spec_elec),
      x = "Time (ms)",
      y = "Absolute t-value",
      color = "") +
    annotate("text",
             x = Inf, y = Inf,
             label = paste("Onset:", df_p2$onset_val[1], "ms"),
             hjust = 1, vjust = 1.5) +
    theme_light()
  
  ggsave(
    filename = file.path(fig_dir, paste0(sub, "_fig_ttest_spec.png")),
    plot = p2, width = 6, height = 4, dpi = 300)
  
  
  # Plot 3: t²-values and permutation threshold for MAX Correction (max t2)
  
  df_perm_long <- max_t2_perm_bn %>%
    mutate(Xf = Xf) %>%
    pivot_longer(
      cols = -Xf,
      names_to = "perm_id",
      values_to = "value")
  df_p3 <- data.frame(
    Xf = Xf,
    t2 = max_t2_bn,
    sig = max_t2_bn > max_t2_max_thresh,
    max_thresh =  max_t2_max_thresh,
    onset_val = df_onsets$maxt2_MAX[df_onsets$participant == sub])
  
  p3 <- ggplot() +
    geom_line(data = df_perm_long, aes(x = Xf, y = value, group = perm_id, color = "Permutation t² values"), alpha = 0.2) +
    geom_line(data = df_p3, aes(x = Xf, y = t2, color = "Squared t² values")) +
    geom_hline(aes(yintercept = max_t2_max_thresh, color = "MAX threshold (95th pct)")) +
    geom_point(
      data = subset(df_p3, sig),
      aes(x = Xf, y = 0, color = "Significant timepoints"),
      size = 0.5) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
    geom_vline(data = df_p3, aes(xintercept = onset_val, color = "First significant point"), linetype = "dashed") +
    scale_color_manual(
      values = c(
        "Squared t² values" = "red",
        "Permutation t² values" = "grey",
        "MAX threshold (95th pct)" = "blue",
        "Significant timepoints" = "green",
        "First significant point" = "black"),
      breaks = c(
        "Squared t² values",
        "Permutation t² values",
        "MAX threshold (95th pct)",
        "Significant timepoints",
        "First significant point")) +
    labs(
      title = paste0(comp, " | P", sub, " | Max Stats Correction at max t² electrode"),
      x = "Time (ms)",
      y = "t²-value",
      color = "") +
    annotate("text",
             x = Inf, y = Inf,
             label = paste("Onset:", df_p3$onset_val[1], "ms"),
             hjust = 1, vjust = 1.5) +
    theme_light()
  
  ggsave(
    filename = file.path(fig_dir, paste0(sub, "_fig_max_t2_MAXthresh.png")),
    plot = p3, width = 6, height = 4, dpi = 300)
  
  # Plot 4: t²-values and permutation threshold for MAX Correction (spec t2)
  
  df_perm_long_spec <- spec_t2_perm_bn %>%
    mutate(Xf = Xf) %>%
    pivot_longer(
      cols = -Xf,
      names_to = "perm_id",
      values_to = "value")
  
  df_p4 <- data.frame(
    Xf = Xf,
    t2 = spec_t2_bn,
    sig = spec_t2_bn > spec_t2_max_thresh,
    max_thresh = spec_t2_max_thresh,
    onset_val = df_onsets$spec_t2_MAX[df_onsets$participant == sub])
  
  p4 <- ggplot() +
    geom_line(data = df_perm_long_spec, aes(x = Xf, y = value, group = perm_id, color = "Permutation t² values"), alpha = 0.2) +
    geom_line(data = df_p4, aes(x = Xf, y = t2, color = "Squared t² values")) +
    geom_hline(aes(yintercept = spec_t2_max_thresh, color = "MAX threshold (95th pct)")) +
    geom_point(
      data = subset(df_p4, sig),
      aes(x = Xf, y = 0, color = "Significant timepoints"),
      size = 0.5) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
    #geom_vline(data = df_p4, aes(xintercept = onset_val, color = "First significant point"), linetype = "dashed") +
    (if (!all(is.na(df_p4$onset_val))) geom_vline(data = df_p4, aes(xintercept = onset_val, color = "First significant point"), linetype = "dashed") else NULL) +
    scale_color_manual(
      values = c(
        "Squared t² values" = "red",
        "Permutation t² values" = "grey",
        "MAX threshold (95th pct)" = "blue",
        "Significant timepoints" = "green",
        "First significant point" = "black"),
      breaks = c(
        "Squared t² values",
        "Permutation t² values",
        "MAX threshold (95th pct)",
        "Significant timepoints",
        "First significant point")) +
    labs(
      title = paste0(comp, " | P", sub, " | Max Stats Correction at electrode ", spec_elec),
      x = "Time (ms)",
      y = "t²-value",
      color = "") +
    annotate("text",
             x = Inf, y = Inf,
             label = paste("Onset:", df_p4$onset_val[1], "ms"),
             hjust = 1, vjust = 1.5) +
    theme_light()
  
  ggsave(
    filename = file.path(fig_dir, paste0(sub, "_fig_spec_t2_MAXthresh.png")),
    plot = p4, width = 6, height = 4, dpi = 300)
  
  # Plot 5: Cluster-based permutation testing
  df_p5 <- data.frame(
    Xf = Xf,
    t2 = max_t2_bn,
    sig = as.logical(max_t2_cs_test),
    onset_val = df_onsets$cluster_sum[df_onsets$participant == sub])
  
  p5 <- ggplot() +
    geom_line(data = df_p5, aes(x = Xf, y = t2, color = "Squared t² values")) +
    #geom_hline(aes(yintercept = average_res$max_t2_cs_thresh, color = "Cluster-sum threshold")) +
    geom_point(
      data = subset(df_p5, sig),
      aes(x = Xf, y = 0, color = "Significant timepoints"),
      size = 1.5) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
    geom_vline(data = df_p5, aes(xintercept = onset_val, color = "Cluster onset"), linetype = "dashed") +
    scale_color_manual(
      values = c(
        "Squared t² values" = "red",
        #"Cluster-sum threshold" = "blue",
        "Significant timepoints" = "green",
        "Cluster onset" = "black"), 
      breaks = c(
        "Squared t² values",
        #"Cluster-sum threshold",
        "Significant timepoints",
        "Cluster onset")) +
    labs(
      title = paste0(comp, " | P", sub, " | Cluster-sum Permutation Testing at max t² electrode"),
      x = "Time (ms)",
      y = "t²-value",
      color = "") +
    annotate("text",
             x = Inf, y = Inf,
             label = paste("Onset:", df_p5$onset_val[1], "ms"),
             hjust = 1, vjust = 1.5) +
    theme_light()
  
  ggsave(
    filename = file.path(fig_dir, paste0(sub, "_fig_cluster_sum.png")),
    plot = p5, width = 6, height = 4, dpi = 300)
  
  cat("Saved all plots for subject", sub, "\n")
  
} # end SUB loop
} # end components loop