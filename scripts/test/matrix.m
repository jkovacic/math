% An Octave script that reproduces expected results of the test module
% 'test/matrixTest.cpp'.
%
% From a shell, run the script as:
%   octave /path/to/matrix.m
%
% From Octave, the script can be run as:
%   matrix


function matrixTest
    f = [1, 0.5, 4.5; 0, 1, 0.4]
    printf('f(0,1) = %f\n', f(0+1, 1+1));
    printf('f(4) = %f\n', f(4+1));

    printf('\nf multiplied by 3:\n');  3 * f
    printf('f conjugated:\n');  conj(f)
    printf('f multiplied by 0.5:\n');  f * 0.5

    printf('t = f transposed:\n');
    t = f'
    printf('rank(t): %f\n', rank(t));

    printf('t * f:\n');  t * f
    printf('f * t:\n');  f * t

    printf('3x3 unit matrix:\n')
    a = eye(3)
    m1 = [1, 2.4, -1.4, 0; 4.5, 1, 0, -0.5; 0, 1.75, 1, 2]
    printf('rank(m1): %f\n', rank(m1));
    printf('removed the 3rd column from m1:\n');
    m1(:, 3) = []
    printf('inserted a column of zeros between the 1st and the 2nd column of m1:\n');
    m1 = [m1(:, 1:1), zeros(3, 1), m1(:, 2:end)]
    printf('removed the 2nd row of m1\n');
    m1(2, :) = []
    printf('\inserted a row (zeroes) between the 1st and 2nd row of m1\n');
    m1 = [m1(1:1, :); zeros(1, 4); m1(2:end, :)]
    printf('rank: %f\n', rank(m1));
    printf('previous matrix transposed:\n');
    m1'

    sq = [1, 2; 3, 4]; 
    kv = sq
    printf('-kv:\n');  -kv
    m3 = -kv;

    a1 = [1, 2, 3; 4, 5, 6; 7, 9, 8]
    printf('Determinant of a1: %f\n', det(a1));
    printf('Rank of a1: %f\n', rank(a1));
    printf('\ninverse of a1:\n');
    invv = inv(a1)
    printf('a1 * inv   (must be a unit matrix):\n');  a1 * invv

    printf('Upper triangular part (incl. diag) of a1:\n');
    triu(a1)
    printf('Lower triangular part (excl. diag) of a1:\n');
    tril(a1, -1)
    printf('Diagonal part of a1:\n');
    m1 = diag(diag(a1))

    printf('invv transposed:\n');
    invv = invv'
    printf('Add 0.5 to inv(2, 0):\n');
    invv(2+1, 0+1) = invv(2+1, 0+1) + 0.5
    printf('Swap the 2nd and the 3rd row:\n');
    invv([2, 3], :) = invv([3, 2], :)
    printf('Swap the 2nd and the 3rd column:\n');
    invv(:, [2, 3]) = invv(:, [3, 2])
    printf('Maximums of each row of invv:');
    max(invv, [], 2)
    printf('Absolute maximums of each row of invv:');
    max(abs(invv), [], 2)
    printf('Minimums of each row of invv:');
    min(invv, [], 2)
    printf('Absolute minimums of each row of invv:');
    min(abs(invv), [], 2)
    printf('Maximums of each column of invv:');
    max(invv, [], 1)
    printf('Absolute maximums of each column of invv:');
    max(abs(invv), [], 1)
    printf('Minimums of each column of invv:');
    min(invv, [], 1)
    printf('Absolute minimums of each column of invv:');
    min(abs(invv), [], 1)
    a1 = 5.5
    printf('\n');

    c1 = [1+i, 1-2i; 2-3i, 2+4i]
    printf('c1 conjugated:\n');
    conj(c1)
    fm = [1.1, 2.2, 3.3, 4.4]'
    printf('fm transposed:\n');
    fmt = fm'
    printf('fm transposed transposed:\n');
    fmt = fmt'
end


matrixTest
