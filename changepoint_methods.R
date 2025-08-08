# Compare Change Point Detection Algorithms
  # obtain onset results using algorithms provided in the changepoint R package 
  # implement CPD both using Segment Neighbor (SN method) and PELT method

  # runs only on P3 component data

library(R.matlab)
library(changepoint)
library(tidyverse) # includes rlang and dplyr
library(writexl)     


SUB <- as.character(1:40)
alpha <- 0.05 # specify alpha level of significance
pth <-  1 - alpha # percentile threshold for MAX Correction

# Run script on data of ERP Component P3
comp <- 'P3'
chanIdx <- 13
resultsDIR <- "E:/Preprocessed_ERP_Core/P3/perm_results"

# Create functions  for running CPD

# Binary Segmentation (BS) method proposed by Auger and Lawrence (1989)
get_BS_cpts <- function(data, Xf) {
  cp <- cpt.meanvar(data, method = "BinSeg", Q = 2, penalty = "Asymptotic", pen.value = alpha)
  # Q specifies maximum of Q change points in "BinSeg"
  cps <- cpts(cp) # extract indices of CPs from object "cp" returned by cpt.meanvar()
  return(Xf[cps[1]]) # return only the first CP
}

# Segment Neighbourhoods (SN) method proposed by Auger and Lawrence (1989)
get_SN_cpts <- function(data, Xf) {
  cp <- cpt.meanvar(data, method = "SegNeigh", Q = 3, penalty = "Asymptotic", pen.value = alpha)
  # Q specifies maximum of Q segments in "SegNeigh"
  cps <- cpts(cp) # extract indices of CPs from object "cp" returned by cpt.meanvar()
  return(Xf[cps[1]]) # return only the first CP
}

# Pruned Exact Linear Time (PELT) method 
# R. Killick, P. Fearnhead & I. A. Eckley (2012) Optimal Detection of Changepoints With a Linear Computational Cost, Journal of the American Statistical Association, 107:500, 1590-1598, DOI: 10.1080/01621459.2012.737745()
get_asymp_PELT_cpts <- function(data, Xf) {
  cp <- cpt.meanvar(data, method = "PELT", penalty = "Asymptotic", pen.value = alpha)
  # if "Asymptotic" is specified, the theoretical type I error is contained in pen.value
  # pen.value = the theoretical type I error e.g.,0.05 when using the "Asymptotic" penalty
  cps <- cpts(cp) # extract indices of CPs from object "cp" returned by cpt.meanvar()
  return(Xf[cps[1]]) # return only the first CP
}

get_manual_PELT_cpts <- function(data, Xf) {
  cp <- cpt.meanvar(data, method = "PELT", penalty = "Manual", pen.value = 50)
  # Pruned Exact Linear Time (PELT) allows to control manual penalty value
  # the higher the penalty value, the greater the cost adding another change-point
  # lower pen.value for more sensitive detection of more change-points
  cps <- cpts(cp) # extract indices of CPs from object "cp" returned by cpt.meanvar()
  return(Xf[cps[1]]) # return only the first CP
}

# Initialize result tables
onset_res <- tibble(
  participant = SUB,
  t2_cpd_prestim = vector("list", length(SUB)),
  maxt2_cpd = vector("list", length(SUB)),
  spec_t2_cpd = vector("list", length(SUB)),)

for (i in seq_along(SUB)) {
  # i <- 1 # for testing purposes
  sub <- SUB[i]
  filepath <- file.path(resultsDIR, paste0(sub, "_results.mat"))
  res <- readMat(filepath)$res
  field_names <- unlist(attr(res, "dimnames")[[1]])
  names(res) <- field_names
  
  # Extract data
  Xf <- as.numeric(res$Xf)
  maxt2_bn <- as.numeric(res$maxt2.bn)
  spec_t2_bn <- res$t2.bn[chanIdx, ] # from electrode of interest
  
  # Select Post-stimulus data
  postStimIdx <- Xf >= 0
  Xf_post <- Xf[postStimIdx]
  maxt2_bn_post <- maxt2_bn[postStimIdx]
  spec_t2_bn_post <- spec_t2_bn[postStimIdx]
  
  # Pre-stim CPD
  onset_res$t2_cpd_prestim[[i]] <- get_BS_cpts(maxt2_bn, Xf)
  
  # CPD excluding baseline period
  onset_res$maxt2_cpd[[i]] <- get_BS_cpts(maxt2_bn_post, Xf_post)
  onset_res$spec_t2_cpd[[i]] <- get_BS_cpts(spec_t2_bn_post, Xf_post)
  
  
end # subject loop
