% An Octave script that reproduces expected results of the test module
% 'test/quatTest.cpp'.
%
% From a shell, run the script as:
%   octave /path/to/quat.m
%
% From Octave, the script can be run as:
%   quat

% Note: the script requires the library 'quaternion'.
% It is available at:
% http://octave.sourceforge.net/quaternion/

pkg load quaternion

zeroq = quaternion(0);
printf('Zero: '); zeroq
q = quaternion(1.0, -0.5, 1.2, -2.3)
printf('q: 1:%f   i:%f   j:%f   k:%f\n', ...
    get(q, 'w'), get(q, 'x'), get(q, 'y'), get(q, 'z') );
printf('abs(q): %f\n', abs(q));
qu = unit(q);
printf('Unit quaternion of q: ');  qu
printf('norm(qu): %f\n', norm(qu));
qu = quaternion(2.3)

o = quaternion(1);
i = qi();
j = qj();
k = qk();
printf('i*j: ');  i * j
printf('j*i: ');  j * i
printf('j*k: ');  j * k
printf('k*j: ');  k * j
printf('k*i: ');  k * i
printf('i*k: ');  i * k

p = quaternion(1, 2, 3, 4)
q
printf('p+q: ');  p + q
printf('p-q: ');  p - q
printf('p*q: ');  p * q
printf('q*p: ');  q * p
printf('p+5: ');  p + 5
printf('3+p: ');  3 + p
printf('p-4: ');  p - 4
printf('1-p: ');  1 - p
printf('p+2: ');  p += 2
printf('p+2-2: ');  p -= 2

qc = -0.5 * (q +i*q*i + j*q*j + k*q*k);
printf('-0.5*(q+i*q*i+j*q*j+k*q*k): ');  qc
printf('conj(q): ');  conj(q)

qrec = inv(q)
printf('q * qrec: ');  q * qrec
printf('qrec * q: ');  qrec * q
