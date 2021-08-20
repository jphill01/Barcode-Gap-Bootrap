set.seed(0673227)

boot.mn <- function(intra, inter, statistic = c("barcode.gap", "min.inter", "max.intra"), m, B = 10000, 
                    replacement = TRUE, conf.level = 0.95) {
  
  boot.samples <- numeric(B) # pre-allocate storage vector
  genetic.dists <- na.omit(cbind(intra, inter)) # remove singleton specimens (those with NAs), if any
  N <- length(genetic.dists)
  
  for (i in 1:B) {
    # sample m genetic distances with or without replacement
    if (replacement == TRUE) {
      intra.boot <- sample(genetic.dists[, "intra"], size = m, replace = TRUE)
      inter.boot <- sample(genetic.dists[, "inter"], size = m, replace = TRUE)
    } else {
      intra.boot <- sample(genetic.dists[, "intra"], size = m, replace = FALSE)
      inter.boot <- sample(genetic.dists[, "inter"], size = m, replace = FALSE)
    }
  
  statistic <- match.arg(statistic
                         )
  if (statistic == "barcode.gap") {
    boot.samples[i] <- min(inter.boot) - max(intra.boot) # barcode gap
  } else if (statistic == "min.inter") {
    boot.samples[i] <- min(inter.boot) # minimum intraspecific distance
  } else {
    boot.samples[i] <- max(intra.boot) # minimum intraspecific distance
    }
  }
  
  stat.obs <- min(genetic.dists[, "inter"]) - max(genetic.dists[, "intra"]) # observed statistic
  stat.boot.est <- mean(boot.samples) # bootstrap mean
  stat.boot.bias <- stat.boot.est - stat.obs # bootstrap bias
  stat.boot.se <- sqrt(sum((boot.samples - stat.boot.est)^2) / (B - 1)) # bootstrap standard error

  ### Bootstrap confidence intervals ###
  
  idx <- trunc((B + 1) * c((1 - conf.level) / 2, (1 + conf.level) / 2))
  z.crit <- qnorm(c((1 - conf.level) / 2, (1 + conf.level) / 2)) # z critical values
  
  stat.boot.norm.ci <- (stat.obs - stat.boot.bias) + z.crit * sd(boot.samples) # Normal
  stat.boot.basic.ci <- rev(2*stat.obs - sort(boot.samples)[idx]) # Basic
  stat.boot.perc.ci <- sort(boot.samples)[idx] # Percentile
  
  ### Output ###
  
  par(mfrow = c(1, 2))
  
  hist(boot.samples) # histogram
  abline(v = stat.boot.est, lty = 2)
  qqnorm(boot.samples) # QQ plot
  qqline(boot.samples)
  
  return(list(genetic.dists = genetic.dists,
              boot.samples = boot.samples, 
              stat.obs = stat.obs,
              stat.boot.est = stat.boot.est,
              stat.boot.bias = stat.boot.bias,
              stat.boot.se = stat.boot.se,
              stat.boot.norm.ci = stat.boot.norm.ci,
              stat.boot.basic.ci = stat.boot.basic.ci,
              stat.boot.perc.ci = stat.boot.perc.ci))
}

# Real data

library(pegas)
library(spider)

data(anoteropsis)
anoDist <- dist.dna(anoteropsis)
anoSpp <- sapply(strsplit(dimnames(anoteropsis)[[1]], split = "_"), 
                 function(x) paste(x[1], x[2], sep = "_"))

data(dolomedes)
doloDist <- dist.dna(dolomedes)
doloSpp <- substr(dimnames(dolomedes)[[1]], 1, 5)


intra <- maxInDist(doloDist, doloSpp)
inter <- nonConDist(doloDist, doloSpp)

intra <- maxInDist(anoDist, anoSpp)
inter <- nonConDist(anoDist, anoSpp)

(out <- boot.mn(intra = intra, inter = inter, statistic = "max.intra", m = 8, B = 10000, replacement = TRUE, conf.level = 0.95))
plot(density(out$boot.samples))


### Compare to standard bootstrapping ###

library(boot)

f <- function(x, i) {
  return(min(x[i, 2]) - max(x[i, 1]))
}

(y <- boot(out$genetic.dists, f, R = 10000))
plot(y)
boot.ci(y, type = c("norm", "basic", "perc"))
 
        
