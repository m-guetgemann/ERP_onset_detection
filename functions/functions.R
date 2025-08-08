# Generate permutation distributions of t values for independent groups
## cond1 and cond2 are matrices (Nf time frames x Nt trials)
permtdist <- function(cond1, cond2, Nf, Nt, nboot = 2000){
  
  st2 <- Nt+1 # index to start condition 2 
  en2 <- Nt+Nt # index to end condition 2
  all <- transpose(cbind(cond1, cond2)) # combine trials, then transpose -> trials x time points
  
  perm.tvals <- matrix(data = 0, nrow = nboot, ncol = Nf)
  # perm.pvals <- matrix(data = 0, nrow = nboot, ncol = Nf)

  for(B in 1:nboot){
    
    permdata <- all[sample(1:en2, size = en2, replace = FALSE),]
    permcond1 <- permdata[1:Nt,]
    permcond2 <- permdata[st2:en2,]

    perm.tvals[B,] <- Rfast::ttests(permcond1, permcond2, paired=FALSE)[,1]    
  }
  # list(tvals = perm.tvals, pvals = perm.pvals)
  perm.tvals
}

# Generate permutation distributions of t values for independent groups
## cond1 and cond2 are matrices (Ne electrodes x Nf time frames x Nt trials)
permtdist.elec <- function(cond1, cond2, Ne, Nf, Nt, nboot = 2000){
  
  st2 <- Nt+1 # index to start condition 2 
  en2 <- Nt+Nt # index to end condition 2
  
  # combine data
  all.data <- array(0, dim = c(Ne,Nf,Nt*2))
  all.data[,,1:Nt] <- cond1
  all.data[,,st2:en2] <- cond2
    # transpose(cbind(cond1, cond2)) # combine trials, then transpose -> trials x time points
  
  perm.tvals <- array(data = 0, dim = c(Ne, Nf, nboot))
  
  for(B in 1:nboot){
    
    # permutation: apply same indices to every electrode and time point
    permdata <- all.data[,,sample(1:en2, size = en2, replace = FALSE)]
    permcond1 <- permdata[,,1:Nt]
    permcond2 <- permdata[,,st2:en2]
    
    for(E in 1:Ne){
    perm.tvals[E,,B] <- Rfast::ttests(transpose(permcond1[E,,]), transpose(permcond2[E,,]), paired=FALSE)[,1]    
    }
  }
  perm.tvals
}

# Generate permutation distributions of t and F values for independent groups:
## t-test on amplitudes and MANOVA on amplitudes and gradients.
## cond1 and cond2 are matrices (Nf time frames x Nt trials)
## cond1.grad and cond2.grad are the same size as cond1 and cond2.
permtfdist <- function(cond1, cond2, cond1.grad, cond2.grad, Nf, Nt, nboot = 2000){
  
  st2 <- Nt+1 # index to start condition 2 
  en2 <- Nt+Nt # index to end condition 2
  all.amp <- transpose(cbind(cond1, cond2)) # combine trials, then transpose -> trials x time points
  all.grad <- cbind(cond1.grad, cond2.grad) # combine trials only
  
  perm.t.amp <- matrix(data = 0, nrow = nboot, ncol = Nf)
  perm.f <- matrix(data = 0, nrow = nboot, ncol = Nf)
  
  for(B in 1:nboot){
    
    # permutation indices
    perm.ind <- sample(1:en2, size = en2, replace = FALSE)
    
    # t-test on amplitudes
    perm.amp <- all.amp[perm.ind,]
    perm.amp1 <- perm.amp[1:Nt,]
    perm.amp2 <- perm.amp[st2:en2,]
    perm.t.amp[B,] <- Rfast::ttests(perm.amp1, perm.amp2, paired=FALSE)[,1]    
    
    # MANOVA on amplitudes and gradients
    perm.grad <- all.grad[,perm.ind]
    perm.grad1 <- perm.grad[,1:Nt]
    perm.grad2 <- perm.grad[,st2:en2]
    for(F in 2:Nf){
      perm.f[B,F] <- hotelling.stat(cbind(perm.amp1[,F],perm.grad1[F,]), cbind(perm.amp2[,F],perm.grad2[F,]),
                                 shrinkage=FALSE, var.equal=FALSE)$statistic
    }
  }
  list(t.amp = perm.t.amp, f.amp.grad = perm.f)
}



sim.counter <- function(S, nsim, inc){
  if(S == 1){
    # print(paste(nsim,"iterations:",S))
    cat(nsim,"iterations:",S)
    beep(2)
  }
  if(S %% inc == 0){
    # print(paste("iteration",S,"/",nsim))
    cat(" /",S)
    beep(2)
  }
}

# INPUTS:
#   sigmask = a vector of 0s and 1s indicating points below and above threshold.
#   Xf      = a vector of time points, matching sigmask.
#   rmzero  = option to discard an onset if it belongs to a cluster that starts in the baseline, including zero; default = TRUE.
#
# OUTPUT:
#   onset   = latency of the first cluster in the units of Xf.
find_onset <- function(sigmask, Xf, rmzero = TRUE){
  onset <- NA
  onset <- try(Xf[which(sigmask)[1]], silent = TRUE)
  if(rmzero){ # remove un-realistic cluster that includes zero
    if(is.finite(onset)){
      if(onset == 0){
        cmap <- cluster.make(sigmask)
        cmap[cmap==1] <- 0
        onset <- Xf[which(sigmask*(cmap>0)>0)[1]]
      }
    }
  }
  onset
}

# https://www.statology.org/mode-in-r/
find_mode <- function(x, na.rm = TRUE) {
  if(na.rm){
    x <- na.omit(x)
  }
  u <- unique(x)
  tab <- tabulate(match(x, u))
  u[tab == max(tab)]
}

keeporder <- function(x){
  x <- as.character(x)
  x <- factor(x, levels=unique(x))
  x
}

# code from Rand Wilcox https://osf.io/xhe8u/
hd<-function(x,q=.5,na.rm=TRUE){
  #
  #  Compute the Harrell-Davis estimate of the qth quantile
  #
  #  The vector x contains the data,
  #  and the desired quantile is q
  #  The default value for q is .5.
  #
  if(na.rm)x=elimna(x)
  n<-length(x)
  m1<-(n+1)*q
  m2<-(n+1)*(1-q)
  vec<-seq(along=x)
  w<-pbeta(vec/n,m1,m2)-pbeta((vec-1)/n,m1,m2)  # W sub i values
  y<-sort(x)
  hd<-sum(w*y)
  hd
}

elimna<-function(m){
  #
  # remove any rows of data having missing values
  #
  DONE=FALSE
  if(is.list(m) && is.matrix(m)){
    z=pool.a.list(m)
    m=matrix(z,ncol=ncol(m))
    DONE=TRUE
  }
  if(!DONE){
    if(is.list(m) && is.matrix(m[[1]])){
      for(j in 1:length(m))m[[j]]=na.omit(m[[j]])
      e=m
      DONE=TRUE
    }}
  if(!DONE){
    if(is.list(m) && is.null(dim(m))){ #!is.matrix(m))
      for(j in 1:length(m))m[[j]]=as.vector(na.omit(m[[j]]))
      e=m
      DONE=TRUE
    }}
  if(!DONE){
    m<-as.matrix(m)
    ikeep<-c(1:nrow(m))
    for(i in 1:nrow(m))if(sum(is.na(m[i,])>=1))ikeep[i]<-0
    e<-m[ikeep[ikeep>=1],]
  }
  e
}


