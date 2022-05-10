setwd("/Users/jarrettphillips/desktop")
set.seed(0673227)

boot.mn <- function(intra, inter, statistic = c("barcode.gap", "min.inter", "max.intra"), m, B = 10000, 
                    replacement = TRUE, conf.level = 0.95) {
  
  boot.samples <- numeric(B) # pre-allocate storage vector of bootstrap resamples
  genetic.dists <- na.omit(cbind(intra, inter)) # remove singleton specimens (those with NAs), if any
  N <- nrow(genetic.dists) # number of specimens
  
  for (i in 1:B) {
    # sample m genetic distances with or without replacement
    if (replacement == TRUE) { # bootstrapping
      intra.boot <- sample(genetic.dists[, "intra"], size = m, replace = TRUE)
      inter.boot <- sample(genetic.dists[, "inter"], size = m, replace = TRUE)
    } else { # subsampling
      intra.boot <- sample(genetic.dists[, "intra"], size = m, replace = FALSE)
      inter.boot <- sample(genetic.dists[, "inter"], size = m, replace = FALSE)
    }
    
    statistic <- match.arg(statistic)
    
    if (statistic == "barcode.gap") {
      boot.samples[i] <- min(inter.boot) - max(intra.boot) # bootstrapped barcode gap
      stat.obs <- min(genetic.dists[, "inter"]) - max(genetic.dists[, "intra"]) # observed sample barcode gap
    } else if (statistic == "min.inter") {
      boot.samples[i] <- min(inter.boot) # bootstrapped minimum intraspecific distance
      stat.obs <- min(genetic.dists[, "inter"]) # observed sample minimum interspecfic distance
    } else {
      boot.samples[i] <- max(intra.boot) # bootstrapped minimum intraspecific distance
      stat.obs <- max(genetic.dists[, "intra"]) # observed sample maximum intraspecific distance
    }
  }
  
  stat.boot.est <- mean(boot.samples) # bootstrap mean
  stat.boot.bias <- stat.boot.est - stat.obs # bootstrap bias
  stat.boot.se <- sd(boot.samples)
  
  ### Bootstrap confidence intervals ###
  
  idx <- trunc((B + 1) * c((1 - conf.level) / 2, (1 + conf.level) / 2))
  z.crit <- qnorm(c((1 - conf.level) / 2, (1 + conf.level) / 2)) # z critical values
  
  ## BCa interval ##
  
  z0 <- qnorm(mean(boot.samples <= stat.obs))
  
  I <- rep(NA, N)
  for (i in 1:N) {
    # Remove ith data point
    intra.new <- genetic.dists[, "intra"][-i]
    inter.new <- genetic.dists[, "inter"][-i]
    # Estimate parameter
    if (statistic == "max.intra") {
      jack.est <- max(intra.new)
      jack.mean <- mean(intra.new)
    } else if (statistic == "min.inter") {
      jack.est <- min(inter.new)
      jack.mean <- mean(inter.new)
    } else {
      jack.est <- min(inter.new) - max(intra.new)
      jack.mean <- mean(inter.new - intra.new)
    }
    I[i] <- jack.mean - jack.est
  } 
  
  # Estimate acceleration constant based on jackknife samples
  a.hat <- (sum(I^3) / sum(I^2)^(3/2)) / 6 # sample skewness 
  
  p.adjusted <- pnorm(z0 + (z0 + z.crit) / (1 - a.hat * (z0 + z.crit))) # adjusted z critical value
 
  stat.boot.norm.ci <- (stat.obs - stat.boot.bias) + z.crit * sd(boot.samples) # Normal
  stat.boot.basic.ci <- rev(2 * stat.obs - sort(boot.samples)[idx]) # Basic
  stat.boot.perc.ci <- quantile(boot.samples, c((1 - conf.level) / 2, (1 + conf.level) / 2)) # Percentile
  stat.boot.bca.ci <- quantile(boot.samples, round(p.adjusted, 1)) # BCa
  
  ### Output ###
  
  par(mfrow = c(1, 2))
  
  hist(boot.samples) # histogram
  abline(v = stat.boot.est, lty = 2)
  qqnorm(boot.samples) # QQ plot
  qqline(boot.samples)
  
  return(list(N = N,
              genetic.dists = genetic.dists,
              boot.samples = boot.samples, 
              stat.obs = stat.obs,
              stat.boot.est = stat.boot.est,
              stat.boot.bias = stat.boot.bias,
              stat.boot.se = stat.boot.se,
              stat.boot.norm.ci = stat.boot.norm.ci,
              stat.boot.basic.ci = stat.boot.basic.ci,
              stat.boot.perc.ci = stat.boot.perc.ci,
              stat.boot.bca.ci = stat.boot.bca.ci))
}

# Real data

library(pegas)
library(spider)

data(anoteropsis)
anoDist <- dist.dna(anoteropsis)
anoSpp <- sapply(strsplit(dimnames(anoteropsis)[[1]], split = "_"), 
                 function(x) paste(x[1], x[2], sep = "_"))

intra <- maxInDist(anoDist, anoSpp)
inter <- nonConDist(anoDist, anoSpp)

(out <- boot.mn(intra = intra, inter = inter, statistic = "barcode.gap", m = ceiling(log(nrow(genetic.dists))), B = 10000, replacement = TRUE, conf.level = 0.95))

summary(out$boot.samples)
length(which(out$boot.samples == out$stat.obs)) / 10000

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

