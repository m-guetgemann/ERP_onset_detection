# Obtain Onset Estimates & Plots of Onset
  # runs on all ERP components, over all participants

  # from Change Point Detection
  # from cluster-sum permutation testing (CSPT)
  # from maximum statistics (MAX) correction

  # input: t2 and HT2 results from all participants in folder 'perm_results' 
  # (created in Permutations_and_analysis.m)
  # input: SVM results from all particpants in folder 'SVM_results'
  # (created in SVM_decoding.m)

  # no further external dependencies; functions are defined at the beginning of the script

library(R.matlab)
library(changepoint)
library(tidyverse) # includes rlang and dplyr
library(writexl)     
library(ggplot2)    
library(matrixStats)
library(bootES)
library(changepoint)

DIR <- "E:/Preprocessed_ERP_Core/All_components_files"
COMP <- c('N170','MMN', 'N2pc', 'N400', 'P3', 'LRP', 'ERN')

SUB <- as.character(1:40)
alpha <- 0.05 # specify alpha level of significance
pth <-  1 - alpha # percentile threshold for MAX Correction
set.seed(1)

# ERP Core defined electrodes of interest for each component
component_channels <- list(
  "N170" = 26,
  "MMN" = 20,
  "N2pc" = 9, 
  "P3" = 13,
  "N400" = 14,
  "LRP" = 5, 
  "ERN" = 20)

# Define list of participants to exclude per component
excludeMap <- list(
  N170  = c("1", "5", "16"),
  MMN   = c("7"),
  N2pc  = c("7", "9", "10", "12", "28"),
  P3    = c("6", "9", "10", "30", "35", "40"),
  N400  = c("40"),
  LRP   = c("6", "30", "40"),
  ERN   = c("5", "6", "30", "40"))

# Create function for CPD on mean and variance
get_cpts <- function(data, Xf) {
  cp <- suppressWarnings(cpt.meanvar(data, method = "BinSeg", Q = 2, penalty = "Asymptotic", pen.value = 0.05))
  cps <- cpts(cp) # extract indices of CPs from object "cp" returned by cpt.meanvar()
  return(Xf[cps[1]]) # return only the first CP
}

# Create function to detect significant time points for MAX Correction
find_onset <- function(sigmask, Xf) {
  onset <- NA
  above <- which(sigmask > 0)
  if (length(above) > 0) {
    onset <- Xf[above[1]]
  }
  return(onset)
}

## Functions from https://github.com/GRousselet/onsetsim/blob/main/code/functions.R
# Form clusters using binary vector: pvals < alpha
# clusters must be at least 2 time frames
cluster.make <- function(x){
  y <- rle(x)
  cmap <- vector(mode = "numeric", length = 0)
  nC <- length(y$values) # number of clusters
  indx <- 0 # cluster counter
  for(CL in 1:nC){
    if(y$values[CL] == 0 || y$lengths[CL] == 1){
      val <- 0
    } else {
      indx <- indx + 1
      val <- indx
    }
    cmap <- c(cmap, rep(val, y$lengths[CL]))
  }
  cmap
}

# Save sum for each cluster
# values = statistics (t values)
# cmap = output from cluster.make
cluster.sum <- function(values, cmap){
  csum <- vector(mode = "numeric", length = max(cmap))
  if(max(cmap)>0){
    for(CL in 1:max(cmap)){
      csum[CL] <- sum(values[cmap==CL])
    }
  } else {
    csum <- 0
  }
  csum
}

# Cluster test
# values = statistics (t values)
# cmap = output from cluster.make
# boot.th = bootstrap quantile of distribution of max cluster sums
cluster.test <- function(values, cmap, boot.th){
  csig <- vector(mode = "logical", length = length(cmap))
  if(max(cmap)>0){
    for(CL in 1:max(cmap)){
      csig[cmap==CL] <- sum(values[cmap==CL]) > boot.th
    }
  }
  csig
}

for (comp in COMP) {

  chanIdx <- component_channels[[comp]]
  currentDIR <- file.path(DIR, comp)
  resultsDIR <- file.path(currentDIR, "perm_results")
  
  figureDIR <- file.path(currentDIR, "r_figures")
  svmDIR <- file.path(currentDIR, "SVM_results")
  if (!dir.exists(figureDIR)) dir.create(figureDIR, showWarnings = FALSE, recursive = TRUE)
  if (!dir.exists(svmDIR)) dir.create(svmDIR, showWarnings = FALSE, recursive = TRUE)
  averagesDIR <- file.path(currentDIR, 'r_averages')
  if (!dir.exists(averagesDIR)) dir.create(averagesDIR, showWarnings = FALSE, recursive = TRUE)
  
  exclude <- excludeMap[[comp]]
  
  # Initialize result tables
  if (comp == "N2pc" || comp == 'LRP' ) {
    onset_res <- tibble(
      participant = SUB,
      spec_t2_cpd_prestim = vector("list", length(SUB)),
      spec_HT2_cpd_prestim = vector("list", length(SUB)),
      spec_t2_cpd = vector("list", length(SUB)), 
      spec_HT2_cpd = vector("list", length(SUB)), 
      spec_t2_MAX = vector("list", length(SUB)),
      spec_HT2_MAX = vector("list", length(SUB)),
      cluster_sum = vector("list", length(SUB)),
      svm_cpd_prestim = vector("list",length(SUB)),
      svm_cpd = vector("list",length(SUB)))
  } else {
    onset_res <- tibble(
      participant = SUB,
      maxt2_cpd_prestim = vector("list", length(SUB)),
      maxHT2_cpd_prestim = vector("list", length(SUB)),
      spec_t2_cpd_prestim = vector("list", length(SUB)),
      spec_HT2_cpd_prestim = vector("list", length(SUB)),
      maxt2_cpd = vector("list", length(SUB)),
      spec_t2_cpd = vector("list", length(SUB)), 
      maxHT2_cpd = vector("list", length(SUB)),
      spec_HT2_cpd = vector("list", length(SUB)), 
      maxt2_MAX = vector("list", length(SUB)),
      maxHT2_MAX = vector("list", length(SUB)),
      spec_t2_MAX = vector("list", length(SUB)),
      spec_HT2_MAX = vector("list", length(SUB)),
      cluster_sum = vector("list", length(SUB)),
      svm_cpd_prestim = vector("list",length(SUB)),
      svm_cpd = vector("list",length(SUB)) 
    )
  }
  
  average_res <- tibble(
    participant = SUB,
    cond1_av = vector("list", length(SUB)),
    cond1_uppCI = vector("list", length(SUB)),
    cond1_lowCI = vector("list", length(SUB)),
    cond2_av = vector("list", length(SUB)),
    cond2_uppCI = vector("list", length(SUB)),
    cond2_lowCI = vector("list", length(SUB)),
    diff_av = vector("list", length(SUB)),
    diff_lowCI = vector("list", length(SUB)),
    diff_uppCI = vector("list", length(SUB)),
    spec_tcrit = vector("list", length(SUB)),
    spec_tval = vector("list", length(SUB)),
    spec_t2_max_thresh = vector("list", length(SUB)),
    max_t2_max_thresh = vector("list", length(SUB)),
    max_t2_cs_thresh = vector("list", length(SUB)),
    max_t2_cs_test = vector("list", length(SUB))
  )

  for (i in seq_along(SUB)) {
    
    sub <- SUB[i]
    if (sub %in% exclude) {
      cat(sprintf("Skipping subject %s for component %s\n", sub, comp))
      next
    }
    
    svm_path <- file.path(svmDIR, paste0(sub, "_SVM.mvpc"))
    mvpc <- readMat(svm_path)
    av_decoding <- as.numeric(mvpc[[1]][[25]])
    
    filepath <- file.path(resultsDIR, paste0(sub, "_results.mat"))
    res <- readMat(filepath)$res
    field_names <- unlist(attr(res, "dimnames")[[1]])
    names(res) <- field_names
    
    # Extract data
    Xf <- as.numeric(res$Xf)
    Nf <- length(Xf)
    spec_t2 <- as.numeric(res$spec.t2)   # only in new_perm_res files
    spec_HT2 <- as.numeric(res$spec.HT2.gradient) # only in new_perm_res files
    spec_t2_bn <- as.numeric(res$spec.t2.bn)
    spec_HT2_bn <- as.numeric(res$spec.HT2.bn)
    
    spec_t2_perm <- as.data.frame(res$spec.t2.perm)
    spec_HT2_perm <- as.data.frame(res$spec.HT2.perm)
    spec_t2_perm_bn <- as.data.frame(res$spec.t2.perm.bn)
    spec_HT2_perm_bn <- as.data.frame(res$spec.HT2.perm.bn)
    
    spec_tval <- as.numeric(res$spec.tval)
    spec_tcrit <- as.numeric(res$spec.tcrit)
    
    # Extract and save average wave forms for plotting
    average_res$cond1_av[[i]] <- as.numeric(res$av.c1.spec)
    average_res$cond2_av[[i]] <- as.numeric(res$av.c2.spec)
    average_res$diff_av[[i]] <- as.numeric(res$diff.spec)
    if (comp != "N2pc" & comp != "LRP" ) {
      average_res$cond2_uppCI[[i]] <-  as.numeric(res$ci.upper.c2av[chanIdx, ])
      average_res$cond2_lowCI[[i]] <- as.numeric(res$ci.lower.c2av[chanIdx, ])
      average_res$diff_uppCI[[i]] <- as.numeric(res$ci.upper.diff[chanIdx, ])
      average_res$diff_lowCI[[i]] <- as.numeric(res$ci.lower.diff[chanIdx, ])
      average_res$cond1_uppCI[[i]] <- as.numeric(res$ci.upper.c1av[chanIdx, ])
      average_res$cond1_lowCI[[i]] <- as.numeric(res$ci.lower.c1av[chanIdx, ])
    } else {
      average_res$cond2_uppCI[[i]] <-  as.numeric(res$ci.upper.c2av)
      average_res$cond2_lowCI[[i]] <- as.numeric(res$ci.lower.c2av)
      average_res$diff_uppCI[[i]] <- as.numeric(res$ci.upper.diff)
      average_res$diff_lowCI[[i]] <- as.numeric(res$ci.lower.diff)
      average_res$cond1_uppCI[[i]] <- as.numeric(res$ci.upper.c1av)
      average_res$cond1_lowCI[[i]] <- as.numeric(res$ci.lower.c1av)
    }
    
    if (comp != "N2pc" & comp != "LRP" ) {
      max_tval = as.numeric(res$max.tval)
      max_tcrit = as.numeric(res$max.tcrit)
      
      maxt2 <- as.numeric(res$maxt2)
      maxHT2 <- as.numeric(res$maxHT2.gradient)
      maxt2_bn <- as.numeric(res$maxt2.bn)
      maxHT2_bn <- as.numeric(res$maxHT2.gradient.bn)
      maxt2_perm_bn <- as.data.frame(res$maxt2.perm.bn) 
      maxHT2_perm_bn <- as.data.frame(res$maxHT2.perm.bn) 
    }
    
    # CPD on t2 including the baseline period
    if (comp != "N2pc" & comp != "LRP" ) {
      onset_res$maxt2_cpd_prestim[[i]] <- get_cpts(maxt2_bn, Xf)
      onset_res$maxHT2_cpd_prestim[[i]] <- get_cpts(maxHT2_bn, Xf)
    }
    onset_res$spec_t2_cpd_prestim[[i]] <- get_cpts(spec_t2_bn, Xf)
    onset_res$spec_HT2_cpd_prestim[[i]] <- get_cpts(spec_HT2_bn, Xf)
    
    # CPD on SVM decoding accuracy
    onset_res$svm_cpd_prestim[[i]] <- get_cpts(av_decoding, Xf)
    
    # Disregard the baseline period
    if (comp == "LRP") {
      postStimIdx <- Xf >= -600
    } else if (comp == "ERN") {
      postStimIdx <- Xf >= -200
    } else {
      postStimIdx <- Xf >= 0}
    
    Xf_post <- Xf[postStimIdx]
    Nf_post <- length(Xf_post)
    spec_t2_bn_post <- spec_t2_bn[postStimIdx]
    spec_HT2_bn_post <- spec_HT2_bn[postStimIdx]
    spec_t2_perm_bn_post <- spec_t2_perm_bn[postStimIdx, ]
    spec_HT2_perm_bn_post <- spec_HT2_perm_bn[postStimIdx, ]
    
    onset_res$spec_t2_cpd[[i]] <- get_cpts(spec_t2_bn_post, Xf_post)
    onset_res$spec_HT2_cpd[[i]] <- get_cpts(spec_HT2_bn_post, Xf_post)
    
    av_decoding_post <- av_decoding[postStimIdx]
    onset_res$svm_cpd[[i]] <- get_cpts(av_decoding_post, Xf_post)
    
    if (comp != "N2pc" & comp != "LRP" ) {
      maxt2_bn_post <- maxt2_bn[postStimIdx]
      maxHT2_bn_post <- maxHT2_bn[postStimIdx]
      maxt2_perm_bn_post <- maxt2_perm_bn[postStimIdx, ]
      maxHT2_perm_bn_post <- maxHT2_perm_bn[postStimIdx, ]
      
      onset_res$maxt2_cpd[[i]] <- get_cpts(maxt2_bn_post, Xf_post)
      onset_res$maxHT2_cpd[[i]] <- get_cpts(maxHT2_bn_post, Xf_post)
    }
    
    # MAX Correction
    # get max from each of Nperm permutation columns, across all time-points
    # get the pth quantile of the of the obtained max values distribution
    # extract latency of first time point using find_onset.R custom function
    spec_t2_thresh <- quantile(apply(spec_t2_perm_bn_post, 2, max), pth)
    spec_HT2_thresh <- quantile(apply(spec_HT2_perm_bn_post, 2, max), pth)
    
    onset_res$spec_t2_MAX[i] <- find_onset(spec_t2_bn_post >= spec_t2_thresh, Xf_post)
    onset_res$spec_HT2_MAX[i] <- find_onset(spec_HT2_bn_post >= spec_HT2_thresh, Xf_post)
    
    # save MAX correction threshold for plotting
    average_res$spec_t2_max_thresh[[i]] <- spec_t2_thresh
    
    if (comp != "N2pc" & comp != "LRP" ) {
      max_t2_thresh <- quantile(apply(maxt2_perm_bn_post, 2, max), pth) 
      max_HT2_thresh <- quantile(apply(maxHT2_perm_bn_post, 2, max), pth)
      onset_res$maxt2_MAX[i] <- find_onset(maxt2_bn_post >= max_t2_thresh, Xf_post)
      onset_res$maxHT2_MAX[i] <- find_onset(maxHT2_bn_post >= max_HT2_thresh, Xf_post)
      
      average_res$max_t2_max_thresh[[i]] <- max_t2_thresh
    }
    
    # Cluster-Based Permutation Testing (cluster sum of t2 max electrode)
    if (comp == "N2pc" | comp == "LRP" ) {
      Nperm <- length(spec_t2_perm_bn_post)
      
      # get p-values (p-value = (perm max t2 ≥ max t2) / Nperm)
      perm_pvals <- rowSums(spec_t2_perm_bn_post >= spec_t2_bn_post) + 1
      perm_pvals <- perm_pvals / (Nperm + 1)
      
      # get clusters
      cluster_map <- cluster.make(perm_pvals <= alpha)
      perm_max_sums <- vector(mode = "numeric", length = Nperm)
      perm_t2_transposed <-  t(spec_t2_perm_bn_post) # because functions work on perm as rows
      perm_thresh <- apply(perm_t2_transposed, 2, quantile, probs = pth)
      
    } else { # use max electrode for all other components
      Nperm <- length(maxt2_perm_bn_post)
      perm_pvals <- rowSums(maxt2_perm_bn_post >= maxt2_bn_post) + 1
      perm_pvals <- perm_pvals / (Nperm + 1)
      cluster_map <- cluster.make(perm_pvals <= alpha)
      perm_max_sums <- vector(mode = "numeric", length = Nperm)
      perm_t2_transposed <-  t(maxt2_perm_bn_post) 
      perm_thresh <- apply(perm_t2_transposed, 2, quantile, probs = pth)
    }
    
    for(p in 1:Nperm){
      # threshold permutation t2 values and form clusters
      perm_cmap <- cluster.make(perm_t2_transposed[p,] <= perm_thresh)  
      perm_max_sums[p] <- max(cluster.sum(values = perm_t2_transposed[p,], cmap = perm_cmap))
    }
    
    # cluster-sum threshold
    cs_thresh <- quantile(perm_max_sums, probs = pth)
    # cluster test
    cs_test <- cluster.test(values = perm_t2_transposed, cmap = cluster_map, cs_thresh)
    # onset.cs <- Xf[cs.test][1]
    onset_res$cluster_sum[i] <- find_onset(cs_test, Xf_post)
    
    average_res$max_t2_cs_thresh[[i]] <- cs_thresh
    average_res$max_t2_cs_test[[i]] <- as.numeric(cs_test)
    
    
    # Function to plot onset estimation from each method 
    plot_onset <- function(time_vec, max_vec, spec_vec, perm_df, onset_time,title_prefix, y_label, file_suffix) {
      
      df_max <- data.frame(Time = time_vec, Value = max_vec)
      df_spec <- data.frame(Time = time_vec, Value = spec_vec)
      # Convert perm_df to long format for ggplot
      perm_df$Time <- time_vec
      perm_df_long <- tidyr::pivot_longer(perm_df, cols = -Time, names_to = "Perm", values_to = "Value")
      
      p <- ggplot() +
        # Add all permutation lines in light grey
        geom_line(data = perm_df_long, aes(x = Time, y = Value, group = Perm, color = "Permutated max statistic"), size = 0.2) +
        # Add main line
        geom_line(data = df_max, aes(x = Time, y = Value, color = "Max Electrode"), size = 0.5) +
        # Add spec line
        geom_line(data = df_spec, aes(x = Time, y = Value, color = "Electrode of Interest"), size = 0.5) +
        scale_color_manual(
          values = c("Permutated max statistic" = "lightgrey", "Max Electrode" = "black", "Electrode of Interest" = "green3"),
          name = NULL) +
        # Axes and theme
        geom_vline(xintercept = 0, color = "black") +
        geom_hline(yintercept = 0, color = "black") +
        xlab("Time (ms)") +
        ylab(paste("Baseline normalized", y_label)) +
        theme_light(base_size = 12) + # font size of all text in subplots
        theme(
          legend.position = c(0.85, 0.85),  # top right within plot
          legend.background = element_rect(fill = "white", color = "gray70"),
          legend.key = element_rect(fill = "transparent"),
          legend.key.size = unit(0.4, "cm"),
          legend.spacing.y = unit(0.2, "cm"))
      
      # Add dashed line and title only if onset is detected
      if (!is.na(onset_time)) {
        p <- p +
          geom_vline(xintercept = onset_time, linetype = "dashed", color = "red") +
          ggtitle(paste0(comp," | P", SUB[i], ": ", title_prefix, " at ", onset_time, " ms"))
      } else {
        p <- p +
          ggtitle(paste0(comp, " | P", SUB[i], ": no ", title_prefix, " detected"))
      }
      ggsave(filename = file.path(figureDIR, paste0(sub, "_", file_suffix, ".png")), 
             plot = p, width = 6, height = 4, dpi = 300) # default unit = "in"
    }
    
    # p1: max t2 CPD
    onset <- round(onset_res$maxt2_cpd[[i]][1], 2)
    plot_onset(Xf, maxt2_bn, spec_t2_bn,maxt2_perm_bn, onset, "first CP in max t²", "max t²", "max_t2_cpd")
    
    # p2: max HT2 CPD
    onset <- round(onset_res$maxHT2_cpd[[i]][1], 2)
    plot_onset(Xf, maxHT2_bn, spec_HT2_bn, maxHT2_perm_bn, onset,"first CP in max HT²", "max HT²", "max_HT2_cpd")
    
    # p3: spec t2 CPD
    onset <- round(onset_res$spec_t2_cpd[[i]][1], 2)
    plot_onset(Xf, maxt2_bn, spec_t2_bn,maxt2_perm_bn, onset, "first CP in t²", "single electrode t²", "spec_t2_cpd")

    # p4: spec HT2 CPD
    onset <- round(onset_res$spec_HT2_cpd[[i]][1], 2)
    plot_onset(Xf, maxHT2_bn, spec_HT2_bn, maxHT2_perm_bn, onset,"first CP in HT²", "single electrode HT²", "spec_HT2_cpd")

    # p5: max t2 MAX
    onset <- round(onset_res$maxt2_MAX[[i]][1], 2)
    plot_onset(Xf, maxt2_bn, spec_t2_bn, maxt2_perm_bn, onset,"onset of MAX corrected t²", "max t²", "maxt2_MAX")

    # p6: max HT2 MAX
    onset <- round(onset_res$maxHT2_MAX[[i]][1], 2)
    plot_onset(Xf, maxHT2_bn, spec_HT2_bn, maxHT2_perm_bn, onset,"onset of MAX corrected HT²", "max HT²", "maxHT2_MAX")

    # p7: spec t2 MAX
    onset <- round(onset_res$spec_t2_MAX[[i]][1], 2)
    plot_onset(Xf, maxt2_bn, spec_t2_bn, spec_t2_perm_bn, onset,"onset of MAX corrected t²", "spec t²", "spec_t2_MAX")

    # p8: spec HT2 MAX
    onset <- round(onset_res$spec_HT2_MAX[[i]][1], 2)
    plot_onset(Xf, maxHT2_bn, spec_HT2_bn, spec_HT2_perm_bn, onset,"onset of MAX corrected HT²", "spec HT²", "HT2_MAX")

    # p9: SVM decoding CPD
    onset <- round(onset_res$svm_cpd[[i]][1], 2)
    plot_onset(Xf, av_decoding, , , onset,"onset of SVM decoding accuracy", "average accuracy", "svm_cpd")
    

    # Save group_averages (nested list with 40 subjects x 9 vectors)
    export_df <- tibble(
      time = Xf,
      cond1_av = average_res$cond1_av[[i]],
      cond1_uppCI = average_res$cond1_uppCI[[i]],
      cond1_lowCI = average_res$cond1_lowCI[[i]],
      cond2_av = average_res$cond2_av[[i]],
      cond2_uppCI = average_res$cond2_uppCI[[i]],
      cond2_lowCI = average_res$cond2_lowCI[[i]],
      diff_av = average_res$diff_av[[i]],
      diff_uppCI = average_res$diff_uppCI[[i]],
      diff_lowCI = average_res$diff_lowCI[[i]],
      #spec_tval = spec_tval,
      #spec_tcrit = spec_tcrit,
      spec_t2_max_thresh = average_res$spec_t2_max_thresh[[i]]
    )
    
    if (comp != "N2pc" & comp != "LRP") {
      export_df <- export_df |> mutate(
        #max_tval = max_tval,
        #max_tcrit = max_tcrit,
        max_t2_max_thresh = average_res$max_t2_max_thresh[[i]],
        max_t2_cs_thresh =  average_res$max_t2_cs_thresh[[i]],
        max_t2_cs_test = c(rep(NA, Nf - Nf_post), average_res$max_t2_cs_test[[i]])
      )
    }
    
    write_xlsx(export_df, path = file.path(averagesDIR, paste0("average_", average_res$participant[[i]], ".xlsx")))
    
    cat("Finished participant ", SUB[i], " of component ", comp, "\n")
    
  } # end participant loop
  
  # Unlist the data onset_res table to prepare for export to xlsx
  onset_res_export <- onset_res %>%
    mutate(across(where(is.list), ~ map_dbl(., ~ if (is.null(.x)) NA_real_ else .x[1])))
  write_xlsx(onset_res_export, path = file.path(currentDIR, "r_onsets_all_methods.xlsx"))
  
} # end component loop

print("Script finished running.")


