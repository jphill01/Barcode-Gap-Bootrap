################################################################################
# Computes the sampling distribution of the DNA barcode gap, intraspecific     #
# distance and minimum interspecific distance via the m-out-of-n resampling    #
# using either with-replacement resampling (bootstrapping) or without-         # 
# replacement sampling (subsampling)                                           #
#                                                                              #
# Coded by Jarrett D. Phillips and Robert G. Young, DoIB, University of Guelph #       
################################################################################

library(pegas) # to calculate pairwise genetic distance matrix
library(spider) # to calculate intra- and interspecific genetic distances

boot.mn <- function(intra, inter, statistic = c("barcode.gap", 
                    "min.inter", "max.intra"), m, B = 10000, replacement = TRUE, 
                    conf.level = 0.95) {
  
  #############################################################################
  # intra - intraspecific pairwise genetic distances                          #                                                                   #
  # inter - interspecific pairwise genetic distances                                                                         #
  # model - DNA substitution model - K2P by default                           #
  # statistic - statistic of interest to resample - barcode.gap by default    #                            #``
  # m - resample size - should be less than the number of data observations   #
  # such that both m and n approach infinity, but m/n approaches 0            #
  # B - resample size - 10000 by default                                      #
  # replacement - with- or without- replacement resampling - TRUE by default  #
  # conf.level - confidence level for interval estimation                     #
  #############################################################################
  
  # pre-allocate storage vector of bootstrap resamples
  boot.samples <- numeric(B) 
  # remove singleton specimens (those with NAs), if any
  genetic.dists <- na.omit(cbind(intra, inter)) 
  # number of specimens
  N <- nrow(genetic.dists) 
  # select desired statistic
  statistic <- match.arg(statistic)
  
  for (i in 1:B) {
    # sample m genetic distances with or without replacement
    if (replacement == TRUE) { # bootstrapping
      intra.boot <- sample(genetic.dists[, "intra"], size = m, replace = TRUE)
      inter.boot <- sample(genetic.dists[, "inter"], size = m, replace = TRUE)
    } else { # subsampling
      intra.boot <- sample(genetic.dists[, "intra"], size = m, replace = FALSE)
      inter.boot <- sample(genetic.dists[, "inter"], size = m, replace = FALSE)
    }
    
    if (statistic == "barcode.gap") {
      # bootstrapped barcode gap
      boot.samples[i] <- min(inter.boot) - max(intra.boot)
      # observed sample barcode gap
      stat.obs <- min(genetic.dists[, "inter"]) - max(genetic.dists[, "intra"]) 
    } else if (statistic == "min.inter") {
      # bootstrapped minimum intraspecific distance
      boot.samples[i] <- min(inter.boot) 
      # observed sample minimum interspecfic distance
      stat.obs <- min(genetic.dists[, "inter"]) 
    } else { # max.intra
      # bootstrapped minimum intraspecific distance
      boot.samples[i] <- max(intra.boot) 
      # observed sample maximum intraspecific distance
      stat.obs <- max(genetic.dists[, "intra"]) 
    }
  }
  
  # bootstrap mean
  stat.boot.est <- mean(boot.samples)
  # bootstrap bias
  stat.boot.bias <- stat.boot.est - stat.obs 
  # bootstrap standard error
  stat.boot.se <- sd(boot.samples) 
  
  # Percentile confidence interval
  stat.boot.perc.ci <- quantile(boot.samples, c((1 - conf.level) / 2, 
                                            (1 + conf.level) / 2)) 

  ### Output ###
  
  par(mfrow = c(1, 2))
  
  # histogram
  hist(boot.samples) 
  abline(v = stat.boot.est, lty = 2)
  # QQ plot
  qqnorm(boot.samples) 
  qqline(boot.samples)
  
  return(list(N = N,
              genetic.dists = genetic.dists,
              boot.samples = boot.samples, 
              stat.obs = stat.obs,
              stat.boot.est = stat.boot.est,
              stat.boot.bias = stat.boot.bias,
              stat.boot.se = stat.boot.se,
              stat.boot.perc.ci = stat.boot.perc.ci))
}

# Real data

data(anoteropsis)
# calculate distance matrix with complete deletion of missing/ambiguous data 
dist.mat <- dist.dna(anoteropsis, model = "raw", pairwise.deletion = FALSE) 
# retrieve taxon names
spp <- sapply(strsplit(dimnames(anoteropsis)[[1]], split = "_"), 
              function(x) paste(x[1], x[2], sep = "_")) 

# calculate max intraspecific distance
intra <- maxInDist(dist.mat, spp) 
# calculate interspecfic (nearest neighbour) distance
inter <- nonConDist(dist.mat, spp)  


(out <- boot.mn(intra = intra, inter = inter, statistic = "barcode.gap", 
                m = ceiling(log(length(intra))), B = 10000, replacement = TRUE, 
                conf.level = 0.95))

summary(out$boot.samples)
mean(out$boot.samples >= out$stat.obs)

# Compare to standard bootstrap - does not work

library(boot)

f <- function(x, i) {
  return(min(x[i, 2]) - max(x[i, 1])) # barcode gap
}
y <- boot(out$genetic.dists, f, R = 10000)
plot(y)
boot.ci(y)

summary(y$t0) # all equal to observed barcode gap
length((which(y$t == y$t0))) / 10000 # very high

