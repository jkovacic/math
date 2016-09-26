#!/usr/bin/env Rscript

# A R script that reproduces expected results of the test module
# 'test/selectionTest.cpp'.
#
# From a shell, run the script as:
#   Rscript /path/to/selection.R
#
# From R/RStudio, the script can be run as:
#   source('/path/to/selection.R')

selectionTest <- function()
{
  #data(mtcars)
  mpg <- mtcars$mpg

  cat("Order of mpgs' indices in ascending order:\n")
  print(order(mpg) - 1)

  cat("\nOrder of mpgs' indices in descending order:\n")
  print(order(mpg, decreasing=TRUE) - 1)
  
  cat("Rank table of 'mpg' in ascending order:\n")
  print(rank(mpg, ties.method="min") - 1)
  
  cat("Rank table of 'mpg' in descending order:\n")
  print(rank(-mpg, ties.method="min") - 1)

  cat('\n')
  cat('Min. mpg: ', min(mpg), '\n')
  cat('Max. mpg: ', max(mpg), '\n')
  cat('Index of min. mpg: ', which.min(mpg) - 1, '\n')
  cat('Index of max. mpg: ', which.max(mpg) - 1, '\n')
  
  cat('\n')
  srt <- sort(mpg)
  for (i in 1:length(srt))
  {
    cat(i, '. smallest mpg: ', srt[i], '\n', sep='')
  }

  cat('\n')
  srt <- sort(mpg, decreasing=TRUE)
  for (i in 1:length(srt))
  {
    cat(i, '. largest mpg: ', srt[i], '\n', sep='')
  }
  
  return(0)
}

selectionTest()
