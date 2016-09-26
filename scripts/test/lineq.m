#!/usr/bin/env octave

% An Octave script that reproduces expected results of the test module
% 'test/lineqTest.cpp'.
%
% From a shell, run the script as:
%   octave /path/to/lineq.m
%
% From Octave, the script can be run as:
%   lineq


function lineqTest
    A = [1+i, 2-i, -1+0.5i; i, 0.25+3i, -2-0.5i; 3, 2-3i, -0.5+i]
    b = [1+0.2i; -2-i; 1]
    x = A \ b

    A1 = [4, -1, -1; -1, 1, 7; -2, 6, 1]
    b1 = [3, 1; -6, 7; 9, -3]
    x1 = A1 \ b1

    A2 = [0.4+0.4i, 1-0.7i, 5+3i, 2-i; ...
          6-i, 0.5+0.3i, -0.1+2i, 1-0.7i; ...
          0.5+0.4i, 0.2-0.3i, -1+0.3i, 6-4i; ...
          -0.5+i, -7-2i, 0.7+0.3i, 2-i]
    b2 = [2-i, 1.6+2.4i; ...
          3+3i, 2.7-1.2i; ...
          -2-4i, -0.8-3.1i; ...
          -1+5i, -1.1+0.4i]
    x2 = A2 \ b2 
end


lineqTest
