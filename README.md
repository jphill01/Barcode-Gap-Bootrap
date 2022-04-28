# Barcode-Gap-Bootstrap

Contains an R script to compute the estimated sampling distribution of the DNA barcode gap, maximum intraspecific distance and minimum interspecific distance using the *m*-out-of-*n* bootstrap. Resampling can be done either with or without replacement.

The function 

    boot.mn(intra, inter, statistic = c("barcode.gap", "min.inter", "max.intra"), m, B = 10000, replacement = TRUE, conf.level = 0.95) 

returns a histogram and QQ plot of the estimated sampling distribution for the statistic of interest. Nonparametric confidence intervals are also computed. 

## Choosing *m*

The value of *m* should be chosen such that *m* and *n* each approach 0 and <img src="https://render.githubusercontent.com/render/math?math=\frac{m}{n}"> approaches &#x221e;. Functions that satisfy this property are $m$ = &radic; $n$.
