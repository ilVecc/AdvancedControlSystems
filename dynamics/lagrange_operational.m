function [eq, dyn_op] = lagrange_operational(q, qD, dynq, dynqD, tau, kin, dyn)
%LAGRANGE_OPERATIONAL Compute the operational space dynamic model for the
%given dynamic matrices and kinematics

J = kin.J;
JA = kin.JA;
TA = kin.TA;

% get derivative of JA (JAD) using the functions (dynq/dynqD) and then
% going back to the symbolics (q/qD)
JAdyn = subs(JA, q, dynq);
JAdynD = diff(JAdyn);
JAD = simplify(subs(JAdynD, [dynqD; dynq], [qD; q]));

B = dyn.B;
C = dyn.C;
G = dyn.G;

% the robot is non-redundant, thus these simplified formulas
BA = simplify(pinv(JA')*B*pinv(JA));
CA_XD = simplify((pinv(JA')*C - BA*JAD)*qD);
GA = simplify(pinv(JA')*G);
yA = TA'*pinv(J')*tau; % yA = TA'*ye  where  tau = J'*ye

%syms x y z phi theta psi;
%syms xD yD zD phiD thetaD psiD;
syms xDD yDD zDD phDD thDD psDD real;

%X = [x, y, z, phi, theta, psi]';
%XD = [xD, yD, zD, phiD, thetaD, psiD]';
XDD = [xDD, yDD, zDD, phDD, thDD, psDD]';

eq = BA*XDD + CA_XD + GA == yA;

dyn_op.BA = BA;
dyn_op.CA_XD = CA_XD;
dyn_op.GA = GA;

end

