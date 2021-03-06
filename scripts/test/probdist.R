#!/usr/bin/env Rscript

# A R script that reproduces expected results of the test module
# 'test/probdistTest.cpp'.
#
# From a shell, run the script as:
#   Rscript /path/to/probdist.R
#
# From R/RStudio, the script can be run as:
#   source('/path/to/probdist.R')


# Note: the script requires the library "triangle".
# For more details about the package, see:
#   https://cran.r-project.org/package=triangle
# It can be installed in R via
#   install.packages('triangle')


probdistTest <- function()
{
  library(triangle)

  cat('Normal distribution:\n')
  cat('N(2,3): z for x=6.7: ', (6.7 - 2) / 3, '\n')
  cat('N(2,3): x for z = -1.3: ', 2 - 1.3 * 3, '\n')
  cat('N(2,3) pdf at x=4.5: ', dnorm(4.5, mean=2, sd=3), '\n')
  cat('N(2,3): P(X<1.72): ', pnorm(1.72, mean=2, sd=3), '\n')
  cat('N(2,3): P(X>2.48): ', pnorm(2.48, mean=2, sd=3, lower.tail=FALSE), '\n')
  cat('N(2,3): P(1<X<3):', pnorm(3, mean=2, sd=3) - pnorm(1, mean=2, sd=3), '\n')
  cat('N(2,3): q(p>0.75): ', qnorm(0.75, mean=2, sd=3), '\n')
  cat('N(2,3): q(p<0.52): ', qnorm(0.52, mean=2, sd=3, lower.tail=FALSE), '\n')
  
  cat('\nLog-normal distribution:\n')
  cat('LN(2,3): z for x=6.7: ', (log(6.7) - 2) / 3, '\n')
  cat('LN(2,3): x for z = -1.3: ', exp(2 - 1.3 * 3), '\n')
  cat('LN(2,3) pdf at x=4.5: ', dlnorm(4.5, mean=2, sd=3), '\n')
  cat('LN(2,3): P(X<1.72): ', plnorm(1.72, mean=2, sd=3), '\n')
  cat('LN(2,3): P(X>2.48): ', plnorm(2.48, mean=2, sd=3, lower.tail=FALSE), '\n')
  cat('LN(2,3): P(1<X<3):', plnorm(3, mean=2, sd=3) - plnorm(1, mean=2, sd=3), '\n')
  cat('LN(2,3): q(p>0.75): ', qlnorm(0.75, mean=2, sd=3), '\n')
  cat('LN(2,3): q(p<0.52): ', qlnorm(0.52, mean=2, sd=3, lower.tail=FALSE), '\n')
  
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

  cat('\nTriangular distribution:\n')
  cat('Tr(12, 68, 30): pdf at x = 7.5: ', dtriangle(7.5, a=12, b=68, c=30), '\n')
  cat('Tr(12, 68, 30): pdf at x = 25.3: ', dtriangle(25.3, a=12, b=68, c=30), '\n')
  cat('Tr(12, 68, 30): pdf at x = 52.4: ', dtriangle(52.4, a=12, b=68, c=30), '\n')
  cat('Tr(12, 68, 30): pdf at x = 72: ', dtriangle(72, a=12, b=68, c=30), '\n')
  cat('Tr(12, 68, 30): P(X<6): ', ptriangle(6, a=12, b=68, c=30), '\n')
  cat('Tr(12, 68, 30): P(X<17.8): ', ptriangle(17.8, a=12, b=68, c=30), '\n')
  cat('Tr(12, 68, 30): P(X<39.2): ', ptriangle(39.2, a=12, b=68, c=30), '\n')
  cat('Tr(12, 68, 30): P(X<70): ', ptriangle(70, a=12, b=68, c=30), '\n')
  cat('Tr(12, 68, 30): P(X>4): ', 1 - ptriangle(4, a=12, b=68, c=30), '\n')
  cat('Tr(12, 68, 30): P(X>22.1): ', 1 - ptriangle(22.1, a=12, b=68, c=30), '\n')
  cat('Tr(12, 68, 30): P(X>42.7): ', 1 - ptriangle(42.7, a=12, b=68, c=30), '\n')
  cat('Tr(12, 68, 30): P(X>69): ', 1 - ptriangle(69, a=12, b=68, c=30), '\n')
  cat('Tr(12, 68, 30): P(3<X<5): ', ptriangle(5, a=12, b=68, c=30) - ptriangle(3, a=12, b=68, c=30), '\n')
  cat('Tr(12, 68, 30): P(8<X<14.2): ', ptriangle(14.2, a=12, b=68, c=30) - ptriangle(8, a=12, b=68, c=30), '\n')
  cat('Tr(12, 68, 30): P(7<X<40): ', ptriangle(40, a=12, b=68, c=30) - ptriangle(7, a=12, b=68, c=30), '\n')
  cat('Tr(12, 68, 30): P(9<X<80): ', ptriangle(80, a=12, b=68, c=30) - ptriangle(9, a=12, b=68, c=30), '\n')
  cat('Tr(12, 68, 30): P(22<X<28): ', ptriangle(28, a=12, b=68, c=30) - ptriangle(22, a=12, b=68, c=30), '\n')
  cat('Tr(12, 68, 30): P(27<X<50): ', ptriangle(50, a=12, b=68, c=30) - ptriangle(27, a=12, b=68, c=30), '\n')
  cat('Tr(12, 68, 30): P(25<X<85): ', ptriangle(85, a=12, b=68, c=30) - ptriangle(25, a=12, b=68, c=30), '\n')
  cat('Tr(12, 68, 30): P(35<X<54): ', ptriangle(54, a=12, b=68, c=30) - ptriangle(35, a=12, b=68, c=30), '\n')
  cat('Tr(12, 68, 30): P(31<X<73): ', ptriangle(73, a=12, b=68, c=30) - ptriangle(31, a=12, b=68, c=30), '\n')
  cat('Tr(12, 68, 30): P(69<X<100): ', ptriangle(100, a=12, b=68, c=30) - ptriangle(69, a=12, b=68, c=30), '\n')
  cat('Tr(12, 68, 30): q(p>0): ', qtriangle(0, a=12, b=68, c=30), '\n')
  cat('Tr(12, 68, 30): q(p>0.23): ', qtriangle(0.23, a=12, b=68, c=30), '\n')
  cat('Tr(12, 68, 30): q(p>0.42): ', qtriangle(0.42, a=12, b=68, c=30), '\n')
  cat('Tr(12, 68, 30): q(p>=1): ', qtriangle(1, a=12, b=68, c=30), '\n')
  cat('Tr(12, 68, 30): q(p<=0): ', qtriangle(1-0, a=12, b=68, c=30), '\n')
  cat('Tr(12, 68, 30): q(p<0.15): ', qtriangle(1-0.15, a=12, b=68, c=30), '\n')
  cat('Tr(12, 68, 30): q(p<0.87): ', qtriangle(1-0.87, a=12, b=68, c=30), '\n')
  cat('Tr(12, 68, 30): q(p<1): ', qtriangle(1-1, a=12, b=68, c=30), '\n')

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
  
  cat('\nExponential distribution:\n')
  cat('Exp(0.6): pdf at x=3: ', dexp(3, rate=0.6), '\n')
  cat('Exp(4): pdf at x=0.7: ', dexp(0.7, rate=4), '\n')
  cat('Exp(0.6): P(X<3): ', pexp(3, rate=0.6), '\n')
  cat('Exp(0.6): P(X>1.5): ', pexp(1.5, rate=0.6, lower.tail=FALSE), '\n')
  cat('Exp(4): P(X<2.5): ', pexp(2.5, rate=4), '\n')
  cat('Exp(4): P(X>3.7): ', pexp(3.7, rate=4, lower.tail=FALSE), '\n')
  cat('Exp(0.6): P(0.2<X<4.5): ', pexp(4.5, rate=0.6) - pexp(0.2, rate=0.6), '\n')
  cat('Exp(4): P(0.12<X<1.8): ', pexp(1.8, rate=4) - pexp(0.12, rate=4), '\n')
  cat('Exp(0.6): q(p<0.3): ', qexp(0.3, rate=0.6), '\n')
  cat('Exp(0.6): q(p>0.4): ', qexp(0.4, rate=0.6, lower.tail=FALSE), '\n')
  cat('Exp(4): q(p<0.6): ', qexp(0.6, rate=4), '\n')
  cat('Exp(4): q(p>0.2): ', qexp(0.2, rate=4, lower.tail=FALSE), '\n')
}


probdistTest()
