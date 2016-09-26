#!/usr/bin/env julia

# A Julia script that reproduces expected results of the test module
# 'test/ratmatTest.cpp'.
#
# From a shell, run the script as:
#   julia /path/to/ratmat.jl
#
# From Julia, the script can be run as:
#   include("path/to/ratmat.jl")


function ratmatTest()
    a = [ -1//2  2//5; 3//4  -1//3 ]
    println("a:")
    println(a);  println()

    invv = inv(a)
    println("inverse of a")
    println(invv);  println()

    a += invv
    println("a += invv:")
    println(a);  println()

    invv *= a
    println("determinant of inv * (a+inv): ", det(invv));  println()

    r = 3//4;
    b = a + r
    println("a + 3/4 :");
    println(b);  println()
    b = r + a
    println("3/4 + a :");
    println(b);  println()

    b = a - r
    println("a - 3/4 :");
    println(b);  println()
    b = r - a
    println("3/4 - a :");
    println(b);  println()

    r = 1//2;
    b = a;
    b += r
    println("a + 1/2 :");
    println(b);  println()
    b = a;
    b -= r
    println("a - 1/2 :");
    println(b);  println()
    b = a / r
    println("a / 1/2 :");
    println(b);  println()
    b /= r
    println("a / (1/2)^2 :");
    println(b);  println()

    b = a .* (a+r)
    println("a .* (a+1/2) :");
    println(b);  println()
    b = a;
    b .*= a
    println("a .* a :");
    println(b);  println()
    b = (a-r) ./ a
    println("(a-1/2) ./ a :");
    println(b);  println()
    b = a - r;
    b ./= (a + r)
    println("(a-1/2) ./ (a+1/2) :");
    println(b);  println()
end


ratmatTest()
