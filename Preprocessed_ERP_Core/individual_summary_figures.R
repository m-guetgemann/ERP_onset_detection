# Subject Summary Figure
 # Create a summary figure of all Components for individual Participants

# Set up
library(tidyverse) 
library(ggplot2)
library(readxl)
library(R.matlab)
library(patchwork)

# Specify directories
DIR <- 'E:/Preprocessed_ERP_Core' 
figDIR <- file.path(DIR,'r_indv_summary_figures') # where plots will save to
dir.create(figDIR, showWarnings = FALSE)

COMP <- c('N170', 'MMN', 'N2pc', 'N400', 'P3', 'LRP', 'ERN')
SUB <- as.character(1:40)
# ERP Core defined electrodes of interest for each component
chanIdx <- list(
  "N170" = 26,
  "MMN" = 20,
  "N2pc" = 9, # 8 and 26 but ERP Core uses 9 for analysis
  "N400" = 14,
  "P3" = 13,
  "LRP" = 5, # 5 and 22 but ERP Core uses 5
  "ERN" = 20
)
elecIdx <- list(
  "N170" = "PO8",
  "MMN" = "FCz",
  "N2pc" = "PO7", # 9 and 26 but ERP Core uses 9 for analysis
  "N400" = "CPz",
  "P3" = "Pz",
  "LRP" = "C3", # 5 and 22 but ERP Core uses 5
  "ERN" = "FCz"
)

# Load onset estimation results of all components
df_all_onsets <- setNames(vector("list", length(COMP)), COMP)
for (c in seq_along(COMP)) {
  comp <- COMP[c]
  path <- file.path(DIR, comp, "r_onsets_all_methods.xlsx")
  df_all_onsets[[comp]] <- read_excel(path)
}

# Custom function to plot both t2 and HT2 CPD in one frame
dual_plot_onset <- function(time_vec, t2_vec, HT2_vec, onset_t2, onset_HT2, title_prefix, electrode) {
  df <- data.frame(Time = time_vec, t2 = t2_vec, HT2 = HT2_vec)
  onset_t2 <- round(onset_t2[[i]][1], 2)
  onset_HT2 <- round(onset_HT2[[i]][1], 2)
  
  p <- ggplot(df, aes(x = Time)) +
    geom_line(aes(y = t2, color = "t²")) +
    geom_line(aes(y = HT2, color = "HT²")) +
    geom_vline(xintercept = 0, color = "black") +
    geom_hline(yintercept = 0, color = "black") +
    xlab("Time (ms)") +
    ylab("Baseline normalized value") +
    scale_color_manual(values = c("t²" = "blue2", "HT²" = "green3"), name = NULL,
                         labels = c("t²" = paste(electrode, "t²"), "HT²" = paste(electrode, "HT²"))) +
    theme_light(base_size = 10) + # font size of all text in subplots
    theme(
      legend.position = c(0.85, 0.85),  # top right within plot
      legend.background = element_rect(fill = "white", color = "gray70"),
      legend.key = element_rect(fill = "transparent"),
      legend.key.size = unit(0.4, "cm"),
      legend.spacing.y = unit(0.2, "cm")) +
   ggtitle(paste0(comp, " | ", " t² CP ", onset_HT2, " ms | HT² CP ", onset_t2, " ms")) 
  
  # Plot onset estimates if it is given
  if (!is.na(onset_t2)) {
    p <- p + geom_vline(xintercept = onset_t2, linetype = "dashed", color = "blue2")
  }
  if (!is.na(onset_HT2)) {
    p <- p + geom_vline(xintercept = onset_HT2, linetype = "dashed", color = "green3")
  }
  return(p)
}

#######################################################################

# Subject Loop - Creates and Saves Summary Figure for each subject
for (i in seq_along(SUB)) {
  # i <- 1 # for debugging purposes
  sub <- SUB[i]
  p_subject_comp <- list()

for (c in seq_along(COMP)) {
  comp <- COMP[c]
  chan <- chanIdx[[comp]] # get channel index for current component
  elec <- elecIdx[[comp]] # get electrode name for current component
  
  # load analysis results
  resultsDIR <- file.path(DIR, comp, 'perm_results') # where results are stored
  filepath <- file.path(resultsDIR, paste0(sub, "_results.mat"))
  res <- readMat(filepath)$res
  field_names <- unlist(attr(res, "dimnames")[[1]])
  names(res) <- field_names
  
  # Extract relevant data
  Xf <- as.numeric(res$Xf)
  maxt2 <- as.numeric(res$maxt2)
  maxHT2 <- as.numeric(res$maxHT2.gradient)
  t2_bn <- as.numeric(res$t2.bn)
  maxt2_bn <- as.numeric(res$maxt2.bn)
  HT2_bn <- as.numeric(res$HT2.subsetNe.bn)
  maxHT2_bn <- as.numeric(res$maxHT2.gradient.bn)
  spec_t2_bn <- as.numeric(res$t2.bn[chan,]) # extract specified channel
  spec_HT2_bn <- as.numeric(res$HT2.gradient.bn[chan,]) # extract specified channel
  
  # Extract relevant onsets
  onset_res <- df_all_onsets[[comp]]
  
  # Create a plot of CPD onsets using max t2 and max HT2 
  p1 <- dual_plot_onset(Xf, maxt2_bn, maxHT2_bn, onset_res$maxt2_cpd, onset_res$maxHT2_cpd, comp, "max")
  p1
  # Create a plot of CPD onsets using t2 and HT2 from electrode of interest
  p2 <- dual_plot_onset(Xf, spec_t2_bn, spec_HT2_bn, onset_res$spec_t2_cpd, onset_res$spec_HT2_cpd, comp, elec)
  p2
  p_subject_comp[[comp]] <- p1 + p2
  
  } # end component loop


plot_subject_comp <- wrap_plots(plotlist = p_subject_comp, ncol = 1, heights = rep(1, 7)) +
plot_annotation(title = sprintf("S%s | Individual onset from CPD | max t² and  max HT² | single electrode t² and HT²", sub))

ggsave(paste0(figDIR,"/", sub, "_fig_summary.png"), plot_subject_comp, width = 8, height = 14, dpi = 300)

} # end subject loop




