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
  cat('Variance of mpg w/o Bessel\'s correction: ', sqrt((n-1)/n) * sd(mtcars$mpg), '\n')
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

  cat('\nNormal distribution:\n')
  cat('N(2,3): z for x=6.7: ', (6.7 - 2) / 3, '\n')
  cat('N(2,3): x for z = -1.3: ', 2 - 1.3 * 3, '\n')
  cat('N(2,3) at x=4.5: ', dnorm(4.5, mean=2, sd=3), '\n')
  cat('N(2,3): P(X<1.72): ', pnorm(1.72, mean=2, sd=3), '\n')
  cat('N(2,3): P(X>2.48): ', pnorm(2.48, mean=2, sd=3, lower.tail=FALSE), '\n')
  cat('N(2,3): P(1<X<3):', pnorm(3, mean=2, sd=3) - pnorm(1, mean=2, sd=3), '\n')
  cat('N(2,3): q(p>0.75): ', qnorm(0.75, mean=2, sd=3), '\n')
  cat('N(2,3): q(p<0.52):', qnorm(0.52, mean=2, sd=3, lower.tail=FALSE), '\n')

  cat('\nStudent\'s (a.k.a. T) distribution:\n')
  cat('T(n=10, mu=2, s=1.5): t for x=3: ', (3-2) / (1.5/sqrt(10)), '\n')
  cat('T(n=10, mu=2, s=1.5): x for t=-1.2:', 2 - 1.2 * 1.5/sqrt(10), '\n')
  cat('T(df=5):  pdf at x=2: ', dt(2, df=5), '\n')
  cat('T(df=12): P(X<2): ', pt(2, df=12), '\n')
  cat('T(df=12): P(X>1.1): ', pt(1.1, df=12, lower.tail=FALSE), '\n')
  cat('T(df=12): P(-0.5<X<1):', pt(1, df=12) - pt(-0.5, df=12), '\n')
  cat('T(df=12): q(p>0.75): ', qt(0.75, df=12), '\n')
  cat('T(df=12): q(p<0.52): ', qt(0.52, df=12, lower.tail=FALSE), '\n')

  cat('\nChi^2 distribution:\n')
  cat('ChiSq(df=2) : pdf at x=1.2: ', dchisq(1.2, df=2), '\n')
  cat('ChiSq(df=7) : pdf at x=3.1: ', dchisq(3.1, df=7), '\n')
  cat('ChiSq(df=1): P(X<2.7): ', pchisq(2.7, df=1), '\n')
  cat('ChiSq(df=4): P(X<1.8): ', pchisq(1.8, df=4), '\n')
  cat('ChiSq(df=0.3): P(X>3.4): ', pchisq(3.4, df=0.3, lower.tail=FALSE), '\n')
  cat('ChiSq(df=5):   P(X>1.7): ', pchisq(1.7, df=5, lower.tail=FALSE), '\n')
  cat('ChiSq(df=1.3): P(2<X<3): ', pchisq(3, df=1.3) - pchisq(2, df=1.3), '\n')
  cat('ChiSq(df=4.2): P(2<X<3): ', pchisq(3, df=4.2) - pchisq(2, df=4.2), '\n')
  cat('ChiSq(df=0.75): q(p>0.25): ', qchisq(0.25, df=0.75), '\n')
  cat('ChiSq(df=3.8):  q(p>0.25): ', qchisq(0.25, df=3.8), '\n')
  cat('ChiSq(df=0.8):  q(p<0.25): ', qchisq(0.25, df=0.8, lower.tail=FALSE), '\n')
  cat('ChiSq(df=6):    q(p<0.25): ', qchisq(0.25, df=6, lower.tail=FALSE), '\n')

  cat('\nF distribution:\n')
  cat('F(1,3): pdf at x=2.1: ', df(2.1, df1=1, df2=3), '\n')
  cat('F(4,3): pdf at x=3.5: ', df(3.5, df1=4, df2=3), '\n')
  cat('F(0.7, 2.5): P(X<4): ', pf(4, df1=0.7, df2=2.5), '\n')
  cat('F(2.5, 0.7): P(X<4): ', pf(4, df1=2.5, df2=0.7), '\n')
  cat('F(0.8, 3.5): P(X>3): ', pf(3, df1=0.8, df2=3.5, lower.tail=FALSE), '\n')
  cat('F(3.5, 1.5): P(X>3): ', pf(3, df1=3.5, df2=1.5, lower.tail=FALSE), '\n')
  cat('F(0.5, 0.5): P(1<X<3): ', pf(3, df1=0.5, df2=0.5) - pf(1, df1=0.5, df2=0.5), '\n')
  cat('F(4, 0.2):   P(1<X<3): ', pf(3, df1=4, df2=0.2) - pf(1, df1=4, df2=0.2), '\n')
  cat('F(0.7, 0.3): q(p>0.63): ', qf(0.63, df1=0.7, df2=0.3), '\n')
  cat('F(5, 6):     q(p>0.63): ', qf(0.63, df1=5, df2=6), '\n')
  cat('F(0.3, 0.7): q(p<0.72): ', qf(0.72, df1=0.3, df2=0.7, lower.tail=FALSE), '\n')
  cat('F(6, 5):     q(p<0.72): ', qf(0.72, df1=6, df2=5, lower.tail=FALSE), '\n')

  cat('\nUniform distribution:\n')
  cat('U(1,3): pdf at x=0: ', dunif(0, min=1, max=3), '\n')
  cat('U(1,3): pdf at x=2: ', dunif(2, min=1, max=3), '\n')
  cat('U(1,3): pdf at x=4.5: ', dunif(4.5, min=1, max=3), '\n')
  cat('U(1,3): P(X<0): ', punif(0, min=1, max=3), '\n')
  cat('U(1,3): P(X<1.5): ', punif(1.5, min=1, max=3), '\n')
  cat('U(1,3): P(X<3.5): ', punif(3.5, min=1, max=3), '\n')
  cat('U(1,3): P(X>-2): ', punif(-2, min=1, max=3, lower.tail=FALSE), '\n')
  cat('U(1,3): P(X>2.7): ', punif(2.7, min=1, max=3, lower.tail=FALSE), '\n')
  cat('U(1,3): P(X>4): ', punif(4, min=1, max=3, lower.tail=FALSE), '\n')
  cat('U(1,3): P(0.25<X<1.4): ', punif(1.4, min=1, max=3) - punif(0.25, min=1, max=3), '\n')
  cat('U(1,3): P(1.9<X<3.7): ', punif(3.7, min=1, max=3) - punif(1.9, min=1, max=3), '\n')
  cat('U(1,3): q(p>0.42): ', qunif(0.42, min=1, max=3), '\n')
  cat('U(1,3): q(p>0.57): ', qunif(0.57, min=1, max=3), '\n')
  cat('U(1,3): q(p<0.34): ', qunif(0.34, min=1, max=3, lower.tail=FALSE), '\n')
  cat('U(1,3): q(p<0.81): ', qunif(0.81, min=1, max=3, lower.tail=FALSE), '\n')

  cat('\nBinomial distribution:\n')
  cat('Binom(5, 0.6): pmf at k=2: ', dbinom(2, 5, 0.6), '\n')
  cat('Binom(5, 0.6): exp. value: ', 5 * 0.6, '\n')
  cat('Binom(5, 0.6): variance: ', 5 * 0.6 * (1-0.6), '\n')
  cat('Binom(5, 0.6): std. dev.: ', sqrt( 5 * 0.6 * (1-0.6) ), '\n')
  cat('Binom(20, 0.6): normal approx: ', 20*0.6>=10 && 20*(1-0.6)>=10, '\n')
  cat('Binom(30, 0.6): normal approx: ', 30*0.6>=10 && 30*(1-0.6)>=10, '\n')
  cat('Binom(10, 0.6): P(X<=7): ', pbinom(7, size=10, prob=0.6), '\n')
  cat('Binom(10, 0.6): P(X<7): ', pbinom(6, size=10, prob=0.6), '\n')
  cat('Binom(10, 0.6): P(X>=6): ', pbinom(5, size=10, prob=0.6, lower.tail=FALSE), '\n')
  cat('Binom(10, 0.6): P(X>6): ', pbinom(6, size=10, prob=0.6, lower.tail=FALSE), '\n')
  cat('Binom(10, 0.6): P(X<0): ', pbinom(-1, size=10, prob=0.6), '\n')
  cat('Binom(10, 0.6): P(X<=0): ', pbinom(0, size=10, prob=0.6), '\n')
  cat('Binom(10, 0.6): P(X>0): ', pbinom(0, size=10, prob=0.6, lower.tail=FALSE), '\n')
  cat('Binom(10, 0.6): P(X>=0): ', pbinom(-1, size=10, prob=0.6, lower.tail=FALSE), '\n')
  cat('Binom(10, 0.6): P(X<10): ', pbinom(9, size=10, prob=0.6), '\n')
  cat('Binom(10, 0.6): P(X<=10): ', pbinom(10, size=10, prob=0.6), '\n')
  cat('Binom(10, 0.6): P(X>10): ', pbinom(10, size=10, prob=0.6, lower.tail=FALSE), '\n')
  cat('Binom(10, 0.6): P(X>=10): ', pbinom(9, size=10, prob=0.6, lower.tail=FALSE), '\n')
  cat('Binom(10, 0.6): P(5<=X<=7): ', pbinom(7, size=10, prob=0.6) - pbinom(4, size=10, prob=0.6), '\n')
  cat('Binom(10, 0.6): P(5<X<=7): ', pbinom(7, size=10, prob=0.6) - pbinom(5, size=10, prob=0.6), '\n')
  cat('Binom(10, 0.6): P(4<=X<9): ', pbinom(8, size=10, prob=0.6) - pbinom(3, size=10, prob=0.6), '\n')
  cat('Binom(10, 0.6): P(3<X<8): ', pbinom(7, size=10, prob=0.6) - pbinom(3, size=10, prob=0.6), '\n')
  cat('Binom(10, 0.6): q(p<0.4): ', qbinom(0.4, size=10, prob=0.6), '\n')
  cat('Binom(10, 0.6): q(p<=0.15): ', qbinom(0.15, size=10, prob=0.6) - 1, '\n')
  cat('Binom(10, 0.6): q(p>0.3): ', qbinom(0.3, size=10, prob=0.6, lower.tail=FALSE) + 1, '\n')
  cat('Binom(10, 0.6): q(p>=0.3): ', qbinom(0.3, size=10, prob=0.6, lower.tail=FALSE), '\n')

  cat('\nPoisson\'s distribution:\n')
  cat('Poisson(4): pmf at k=3: ', dpois(3, lambda=4), '\n')
  cat('Poisson(4): P(X<=3): ', ppois(3, lambda=4), '\n')
  cat('Poisson(4): P(X<5): ', ppois(5-1, lambda=4), '\n')
  cat('Poisson(4): P(X>=6): ', ppois(6-1, lambda=4, lower.tail=FALSE), '\n')
  cat('Poisson(4): P(X>2): ', ppois(2, lambda=4, lower.tail=FALSE), '\n')
  cat('Poisson(4): P(X<=0): ', ppois(0, lambda=4), '\n')
  cat('Poisson(4): P(X<0): ', ppois(0-1, lambda=4), '\n')
  cat('Poisson(4): P(X>=0): ', ppois(0-1, lambda=4, lower.tail=FALSE), '\n')
  cat('Poisson(4): P(X>0): ', ppois(0, lambda=4, lower.tail=FALSE), '\n')
  cat('Poisson(4): P(5<=X<=7): ', ppois(7, lambda=4) - ppois(5-1, lambda=4), '\n')
  cat('Poisson(4): P(5<X<=7): ', ppois(7, lambda=4) - ppois(5, lambda=4), '\n')
  cat('Poisson(4): P(5<=X<7): ', ppois(7-1, lambda=4) - ppois(5-1, lambda=4), '\n')
  cat('Poisson(4): P(5<X<7): ', ppois(7-1, lambda=4) - ppois(5, lambda=4), '\n')
  cat('Poisson(4): q(p<0.3): ', qpois(0.3, lambda=4), '\n')
  cat('Poisson(4): q(p<=0.7): ', qpois(0.7, lambda=4) - 1, '\n')
  cat('Poisson(4): q(p>0.4): ', qpois(0.4, lambda=4, lower.tail=FALSE) + 1, '\n')
  cat('Poisson(4): q(p>0.6): ', qpois(0.6, lambda=4, lower.tail=FALSE), '\n')

  cat('\n')
  cat('Median of mpg: ', median(mtcars$mpg), '\n')
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
  alg <- 1
  for (t in types)
  {
    cat(alg, ': ')
    print(quantile(x, p, type=t, names=FALSE))
    alg <- alg + 1
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
  
  return(0)
}

statTest()
