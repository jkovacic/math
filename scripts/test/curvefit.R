#!/usr/bin/env Rscript

# A R script that reproduces expected results of the test module
# 'test/curvefitTest.cpp'.
#
# From a shell, run the script as:
#   Rscript /path/to/curvefit.R
#
# From R/RStudio, the script can be run as:
#   source('/path/to/curvefit.R')


curvefitTest <- function()
{
  x <- seq(0, 5)
  y <- exp(-x)

  p1 <- lm(y ~ x)
  p2 <- lm(y ~ x + I(x^2))
  p3 <- lm(y ~ x + I(x^2) + I(x^3))

  cat("1st degree regression polynomial:\n")
  print(p1)
  cat("\n2nd degree regression polynomial:\n")
  print(p2)
  cat("\n3rd degree regression polynomial:\n")
  print(p3)


  xt = c(-1.5, -0.75, 0, 0.75, 1.5)
  yt <- tan(xt)

  intp <- lm(yt ~ xt + I(xt^2) + I(xt^3))
  cat("\n3rd degree interpolation polynmial for tan(x):\n")
  print(intp)


  p5 <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5))
  cat("\n5th degree interpolation polynomial for exp(-x):\n")
  print(p5)
}

curvefitTest()
