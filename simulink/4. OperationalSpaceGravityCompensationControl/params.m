addpath('../.')

%         [   x    y    z  phi theta psi]
KP = diag([ 250  350  100   50     1   1]);
KD = diag([  35   45   15    1     0   0]);

g = 9.81;