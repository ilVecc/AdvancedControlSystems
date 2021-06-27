addpath('../.')

traj_on = 1;
%        [    x    y    z   phi theta psi]
w = diag([   15   30   30    10   1     1]);
z = diag([ 0.95 0.95 0.95  0.95   0     0]);

KP = w^2;
KD = 2*z*w;

g = 9.81;