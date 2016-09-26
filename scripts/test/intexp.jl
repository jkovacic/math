#!/usr/bin/env julia

# A Julia script that reproduces the second part of expected results of
# the test module 'test/intexpTest.cpp'.
# The other parts are performed in 'scripts/test/intexp.m' and
# 'scripts/test/intexp.mac'.
#
# From a shell, run the script as:
#   julia /path/to/intexp.jl
#
# From Julia, the script can be run as:
#   include("path/to/intexp.jl")


# The test cases are performed in three scripts because each mathematical
# structure is (natively) supported by one programming language/mathematical
# software better than by the others.


# continued from test cases in 'intexp.m'

function intexpTest()
    f = -3//4
    println("f = ", f)
    println("f^0 = ", f^0)
    println("f^7 = ", f^7)
end


intexpTest()

# test cases continue in 'intexp.mac'.
