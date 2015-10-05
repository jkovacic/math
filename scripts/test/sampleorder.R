# A R script that reproduces expected results of the test module
# 'test/sampleorderTest.cpp'.
#
# From a shell, run the script as:
#   Rscript /path/to/sampleorder.R
#
# From R/RStudio, the script can be run as:
#   source('/path/to/sampleorder.R')

data(mtcars)

cat("Order of mpgs' indices in ascending order:\n")
print(order(mtcars$mpg) - 1)

cat("\nOrder of mpgs' indices in descending order:\n")
print(order(mtcars$mpg, decreasing=TRUE) - 1)
