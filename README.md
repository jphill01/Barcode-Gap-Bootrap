# Barcode-Gap-Bootstrap

Contains a script to compute the estimated sampling distribution of the DNA barcode gap, maximum intraspecific distance and minimum interspecific distance using the *m*-out-of-*n* bootstrap. Resampling can be done either with or without replacement.

The function **boot.mn(intra, inter, statistic = c("barcode.gap", "min.inter", "max.intra"), m, B = 10000, replacement = TRUE, conf.level = 0.95)** returns a histogram and QQ plot of the estimated sampling distribution for the statistic of interest. Nonparametric confidence intervals are also computed. 
