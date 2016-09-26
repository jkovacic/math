#!/usr/bin/env julia

# A Julia script that reproduces expected results of the test module
# 'test/rationalTest.cpp'.
#
# From a shell, run the script as:
#   julia /path/to/rational.jl
#
# From Julia, the script can be run as:
#   include("path/to/rational.jl")


function rationalTest()
    s = Rational(-12)
    println("Rational(-12): ", s)
    s = 34824//1000
    println("Rational('34.824'): ", s)
    s = (-4285714+4) // 999999
    println("-4.|285714 = ", s)
    println("abs(", s, ") = ", abs(s))
    s = (35167-35) // (10000-10)
    println("3.5|167 = ", s)
    println("abs(", s, ") = ", abs(s))

    println()
    a = 15//20
    println("a = ", a, " = ", num(a)*20, "//", den(a)*20)
    b = a
    c = Rational(0)
    print("b = ", b, "    c = ", c)
    c = a
    c1 = Rational(-7)
    println("    c = ", c, "   c = ", c1)

    println()
    b = inv(a)
    c = 2//5
    println("a = ", a, "    b = ", b, "    c = ", c, "    c^(-1) = ", inv(c))

    a = 1//3
    b = 1//2
    c = +a + b
    println(a, " + ", b, " = ", c, " = ", float(c))
    b += a
    a = 1//4
    c = b - a
    println(b, " - ", a, " = ", c, " = ", float(c))
    a -= c
    c = a * b
    println(a, " * ", b, " = ", c)
    c /= b
    a = b / c
    println(b, " / ", c, " = ", a)

    c = 0//1
    println("a=", a, "   isZero: ", a==0, "   isPositive: ", a>0, "    isNegative: ", a<0)
    println("b=", b, "   isZero: ", b==0, "   isPositive: ", b>0, "    isNegative: ", b<0)
    println("c=", c, "   isZero: ", c==0, "   isPositive: ", c>0, "    isNegative: ", c<0)

    println(a, " + 2 = ", a + 2)
    println("5 + ", a, " = ", 5 + a)
    println(a, " - 3 = ", a - 3)
    println("1 - ", a, " = ", 1 - a)
    println(a, " * (-4) = ", a * (-4))
    println("2 * ", a, " = ", 2 * a)
    println(a, " / 8 = ", a / 8)
    println("10 / ", a, " = ", 10 / a)

    c = a
    println("a = ", a, "    b = ", b, "    c = ", c)

    if a==c
        println("a == c")
      else
        println("a != c")
    end

    if a==b
        println("a == b")
      else
        println("a != b")
    end

    if a!=c
        println("a != c")
      else
        println("a == c")
    end

    if a!=b
        println("a != b")
      else
        println("a == b")
    end

    println(a, " == 3: ", a == 3)
    println("2 == ", b, ": ", 2 == b)
    println(a, " != 1: ", a != 1)
    println("-2 != ", b, ": ", -2 != b)
    println(a, " > 1: ", a > 1)
    println("3 > ", b, ": ", 3 > b)
    println(a, " >= 2: ", a >= 2)
    println("5 >= ", b, ": ", 5 >= b)
    println(a, " < 6: ", a < 6)
    println("-2 < ", b, ": ", -2 < b)
    println(a, " <= 4: ", a <= 4)
    println("3 <= ", b, ": ", 3 <= b)

    print(a, " + 3 = ")
    a += 3
    println(a)

    print(a, " - 1 = ")
    a -= 1
    println(a)

    print(a, " * 5 = ")
    a *= 5
    println(a)

    print(a, " / 2 = ")
    a /= 2
    println(a)
end


rationalTest()
