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

  cat("Order of mpgs' indices in ascending order:\n")
  print(order(mtcars$mpg) - 1)

  cat("\nOrder of mpgs' indices in descending order:\n")
  print(order(mtcars$mpg, decreasing=TRUE) - 1)
  
  cat("Rank table of 'mpg' in ascending order:\n")
  print(rank(mtcars$mpg, ties.method="min") - 1)
  
  cat("Rank table of 'mpg' in descending order:\n")
  print(rank(-mtcars$mpg, ties.method="min") - 1)

  cat('\n')
  cat('Min. mpg: ', min(mtcars$mpg), '\n')
  cat('Max. mpg: ', max(mtcars$mpg), '\n')
  cat('Index of min. mpg: ', which.min(mtcars$mpg) - 1, '\n')
  cat('Index of max. mpg: ', which.max(mtcars$mpg) - 1, '\n')
  
  return(0)
}

selectionTest()
