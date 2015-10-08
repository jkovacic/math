% An Octave script that reproduces the first part of expected results
% of the test module 'test/intexpTest.cpp'.
% The other parts are performed in 'scripts/test/intexp.jl' and
% 'scripts/test/intexp.mac'.
%
% From a shell, run the script as:
%   octave /path/to/intexp.m
%
% From Octave, the script can be run as:
%   intexp

% Note: the script requires the library 'quaternion'.
% It is available at:
% http://octave.sourceforge.net/quaternion/

pkg load quaternion

printf('5^0 = %f\n', 5^0);
printf('2^10 = %f\n', 2^10);
printf('(-3)^17 = %f\n', (-3)^17);
printf('5^12 = %f\n', 5^12);
printf('(-7)^4 = %f\n', (-7)^4);
printf('sqrt(2)^30 = %f\n', sqrt(2)^30);

printf('\n');
m = [0.1, -0.2, 0.3; -0.4, 0.5, -0.6; 0.7, -0.8, 0.9]
printf('m^0 =\n');
m^0
printf('m^5 =\n');
m^5

printf('\n');
q = quaternion(1, -0.8, 1.2, -1.5)
printf('q^0 = ');
q^0

q10 = q^10; 
printf('q^10 = %f %+fi %+fj %+fk\n', ...
    get(q10, 'w'), get(q10, 'x'), get(q10, 'y'), get(q10, 'z'));

% test cases continue in 'intexp.jl'.
