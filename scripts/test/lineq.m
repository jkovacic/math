% An Octave script that reproduces expected results of the test module
% 'test/lineqTest.cpp'.
%
% From a shell, run the script as:
%   octave /path/to/lineq.m
%
% From Octave, the script can be run as:
%   lineq

A = [1+i, 2-i, -1+0.5i; i, 0.25+3i, -2-0.5i; 3, 2-3i, -0.5+i]
b = [1+0.2i; -2-i; 1]
x = A \ b
