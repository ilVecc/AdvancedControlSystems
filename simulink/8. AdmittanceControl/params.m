addpath('../.')

% motion [   x    y    z   phi theta psi]
w = diag([  25  100   45     1     1   1]);
z = diag([ 0.5  1.5  1.5     0     0   0]);

Md = diag([2 2 2 2 2 2]);
KP = w^2;
KD = 2*z*w;


% compliance  [x y z phi theta psi]
wt =  10*diag([1 1 1   1     1   1]);
zt = 0.9*diag([1 1 1   1     1   1]);

Mt  = 0.75*diag([1 1 1 1 1 1]);
KPt = wt^2;
KDt = 2*zt*wt;


g = 9.81;