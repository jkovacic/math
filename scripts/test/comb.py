# A Python script that reproduces expected results of the test module
# 'test/combTest.cpp'.
#
# From a shell, run the script as:
#   python /path/to/comb.py
#
# From Python, the script can be run as:
#   import comb


import itertools as it

perm = it.permutations("abcde")

print "Permutations...\n"
cnt = 1
for i in list(perm) :
    print cnt, ": ", ''.join(i)
    cnt += 1

print "\nCombinations...\n"
for k in range(1, 6) :
    print "\nK = ", k, "\n"
    comb = it.combinations("abcde", k)
    cnt = 1
    for i in list(comb) :
        print cnt, ": ", ''.join(i)
        cnt += 1
