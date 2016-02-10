# A Python script that reproduces the second part of expected results of the
# test module 'test/specfunTest.cpp'.
# The first part is performed in 'scripts/test/specfun.mac'.
#
# From a shell, run the script as:
#   python /path/to/specfun.py
#
# From Python 2.x, the script can be run as:
#   execfile("/path/to/specfun.py")
# From Python 3.x, the script can be run as:
#   exec(open("/path/to/specfun.py").read())

# The test cases are performed in two scripts because neither Python
# (SciPy) nor Maxima fully support special functions.


# Note: the script requires the library 'scipy'.


from __future__ import print_function


def specfunTest() :
    import scipy.special as sp

    # continued from test cases in 'specfun.mac'

    print("Inverse of reg. lower inc. gamma(0.2, 0.3):  ", sp.gammaincinv(0.2, 0.3))
    print("Inverse of reg. lower inc. gamma(3, 0.7):    ", sp.gammaincinv(3.0, 0.7))
    print("Inverse of reg. upper inc. gamma(0.3, 0.4):  ", sp.gammainccinv(0.3, 0.4))
    print("Inverse of reg. upper inc. gamma(5.2, 0.82): ", sp.gammainccinv(5.2, 0.82))
    print("Inverse of lower inc. gamma(0.24, 0.94):     ", sp.gammaincinv(0.24, 0.94/sp.gamma(0.24)))
    print("Inverse of lower inc. gamma(2.8, 0.17):      ", sp.gammaincinv(2.8, 0.17/sp.gamma(2.8)))
    print("Inverse of upper inc. gamma(0.65, 0.86):     ", sp.gammainccinv(0.65, 0.86/sp.gamma(0.65)))
    print("Inverse of upper inc. gamma(3.5, 0.43):      ", sp.gammainccinv(3.5, 0.43/sp.gamma(3.5)))

    print()
    print("Inverse of reg. lower inc. beta(0.3, 0.2, 0.7):  ", sp.betaincinv(0.3, 0.2, 0.7))
    print("nverse of reg. lower inc. beta(2.4, 3.5, 0.6):   ", sp.betaincinv(2.4, 3.5, 0.6))
    print("Inverse of reg. upper inc. beta(0.9, 1.5, 0.7):  ", sp.betaincinv(0.9, 1.5, 1-0.7))
    print("Inverse of reg. upper inc. beta(1.9, 2.7, 0.25): ", sp.betaincinv(1.9, 2.7, 1-0.25))
    print("Inverse of lower inc. beta(2.8, 0.3, 2):         ", sp.betaincinv(2.8, 0.3, 2/sp.beta(2.8, 0.3)))
    print("Inverse of lower inc. beta(1.1, 1.3, 0.4):       ", sp.betaincinv(1.1, 1.3, 0.4/sp.beta(1.1, 1.3)))
    print("Inverse of upper inc. beta(0.4, 0.5, 1.8):       ", sp.betaincinv(0.4, 0.5, 1-1.8/sp.beta(0.4, 0.5)))
    print("Inverse of upper inc. beta(1.7, 1.1, 0.2):       ", sp.betaincinv(1.7, 1.1, 1-0.2/sp.beta(1.7, 1.1)))

    print()
    print("erf(-1.2): ", sp.erf(-1.2))
    print("erf(0.7):  ", sp.erf(0.7))
    print("erfc(0.2): ", sp.erfc(0.2))

    print()
    print("erfInv(-0.12): ", sp.erfinv(-0.12))
    print("erfInv(0.34):  ", sp.erfinv(0.34))
    print("erfcInv(1.7):  ", sp.erfcinv(1.7))
    print("erfcInv(0.65): ", sp.erfcinv(0.65))


specfunTest()
