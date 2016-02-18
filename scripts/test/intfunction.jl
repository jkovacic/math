# A Julia script that reproduces expected results of the test module
# 'test/intfunctionTest.cpp'.
#
# From a shell, run the script as:
#   julia /path/to/intfunction.jl
#
# From Julia, the script can be run as:
#   include("path/to/intfunction.jl")


function intfunctionTest()
    nums = [ 12, 100, 37423 ]
    for n = nums
       println("floor(sqrt(", n, ")) = ", convert(Integer, floor(sqrt(n))) )
    end
    
    println()
    for n = nums
       println("ceil(sqrt(", n, ")) = ", convert(Integer, ceil(sqrt(n))) )
    end

    println()
    nums = [ 1, 8, 63, 64, 1024, 1025, 65000 ]
    for n = nums
       println("floor(log2(", n, ")) = ", convert(Integer, floor(log2(n))) )
    end

    println()
    for n = nums
       println("ceil(log2(", n, ")) = ", convert(Integer, ceil(log2(n))) )
    end

    println()
    nums = [ 1, 16, 65 ]
    for n = nums
       println("floor(log2(", n, ")) = ", convert(Integer, floor(log2(n))) )
    end

    println()
    for n = nums
       println("ceil(log2(", n, ")) = ", convert(Integer, ceil(log2(n))) )
    end
end


intfunctionTest()
