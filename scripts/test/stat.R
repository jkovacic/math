#!/usr/bin/env Rscript

# A R script that reproduces expected results of the test module
# 'test/statTest.cpp'.
#
# From a shell, run the script as:
#   Rscript /path/to/stat.R
#
# From R/RStudio, the script can be run as:
#   source('/path/to/stat.R')

# Note: the script requires the library "moments".
# It can be installed in R via
#   install.packages('moments')


statTest <- function()
{
  library(moments)
  #data(mtcars)

  cat('Sum of all mpgs: ', sum(mtcars$mpg), '\n')
  cat('Mean mpg: ', mean(mtcars$mpg), '\n')
  cat('Variance of mpg: ', var(mtcars$mpg), '\n')
  cat('Standard deviation of mpg: ', sd(mtcars$mpg), '\n')
  n <- length(mtcars$mpg)
  cat('Variance of mpg w/o Bessel\'s correction: ', (n-1)/n * var(mtcars$mpg), '\n')
  cat('Standard deviation of mpg w/o Bessel\'s correction: ', sqrt((n-1)/n) * sd(mtcars$mpg), '\n')
  cat('4th central moment of mpg about the mean: ', moment(mtcars$mpg, order=4, central=TRUE), '\n')
  cat('4th moment of mpg about the origin:', moment(mtcars$mpg, order=4, central=FALSE), '\n')
  cat('Skewness of mpg: ', skewness(mtcars$mpg), '\n')
  cat('Skewness of mpg (w/ sample sd): ', skewness(mtcars$mpg)*((n-1)/n)^(3/2), '\n')
  cat('Kurtosis of mpg: ', kurtosis(mtcars$mpg), '\n')
  cat('Excess kurtosis of mpg: ', kurtosis(mtcars$mpg) - 3, '\n')
  cat('Proportion of mpgs <=26: ', sum(mtcars$mpg <= 26) / length(mtcars$mpg), '\n')
  cat('Proportion of mpgs <26: ', sum(mtcars$mpg < 26) / length(mtcars$mpg), '\n')

  cat('\n')
  cat('Covariance of mpg and wt: ', cov(mtcars$mpg, mtcars$wt), '\n')
  cat('Covariance of mpg and wt w/o Bessel\'s correction: ', (n-1)/n * cov(mtcars$mpg, mtcars$wt), '\n')
  cat('Pearson\'s r (correlation) of mpg and wt: ', cor(mtcars$mpg, mtcars$wt), '\n')
  cat('r^2 of mpg and wt: ', cor(mtcars$mpg, mtcars$wt)^2, '\n')

  cat('Sum-product of mpg and wt: ', sum(mtcars$mpg * mtcars$wt), '\n')

  cat('\n')
  cat('Median of mpg: ', median(mtcars$mpg), '\n')
  cat('Approx. median of mpg: ', sort(mtcars$mpg)[(n %/% 2L) + 1L], '\n')
  cat('1st and 3rd quartile of mpg:\n')
  print(quantile(mtcars$mpg, c(0.25, 0.75)))
  cat('IQR of mpg: ', IQR(mtcars$mpg), '\n')
  cat('63th percentile of mpg: ', quantile(mtcars$mpg, 0.63), '\n')
  cat('Proportion of mpgs <=25: ', sum(mtcars$mpg <= 25) / length(mtcars$mpg), '\n')
  cat('Proportion of mpgs <=14.5: ', sum(mtcars$mpg <= 14.5) / length(mtcars$mpg), '\n')

  cat('\nQuantiles by various algorithms:\n')
  x <- mtcars$mpg;
  types <- 1:9;
  p <- c(0.01, 0.1, 0.25, 0.375, 0.5, 0.625, 0.75, 0.9, 0.99);
  alg <- 1L
  for (t in types)
  {
    cat(alg, ': ')
    print(quantile(x, p, type=t, names=FALSE))
    alg <- alg + 1L
  }

  cat('\n')
  cat('Min. mpg: ', min(mtcars$mpg), '\n')
  cat('Max. mpg: ', max(mtcars$mpg), '\n')
  N <- length(mtcars$mpg);
  cat('5th smallest mpg: ', sort(mtcars$mpg, partial=4+1)[4+1], '\n')
  cat('5th largest mpg: ', sort(mtcars$mpg, partial=N-4)[N-4], '\n')
  cat('3rd smallest mpg: ', sort(mtcars$mpg, partial=3)[3], '\n')
  cat('6th largest mpg: ', sort(mtcars$mpg, partial=N-6+1)[N-6+1], '\n')

  cat('\n')
  q <- quantile(mtcars$mpg, c(0.25, 0.75), name=FALSE);
  d <- IQR(mtcars$mpg);
  cat('Is each mpg in range between ', q[1]-1.0*d, ' and ', q[2]+1.0*d, ':\n')
  print(mtcars$mpg < (q[1]-1.0*d) | mtcars$mpg > (q[2]+1.0*d))
  cat('Sorted list of outliers in mpg:\n')
  print(sort(unique(mtcars$mpg[mtcars$mpg<(q[1]-0.5*d) | mtcars$mpg>(q[2]+0.5*d)])))
}


statTest()
