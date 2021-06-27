addpath('../.')

%        [   x    y    z  phi theta psi]
w = diag([  25  100   45    1   1     1]);
z = diag([ 0.5  1.5  1.5    0   0     0]);

Md = diag([2 2 2 2 2 2]);
KP = w^2;
KD = 2*z*w;

g = 9.81;