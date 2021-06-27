function [eq, dyn] = lagrange(q, qD, qDD, kin)
%LAGRANGE get robot dynamics using the Lagrange formulation

% everything is calculated w.r.t. frame 0

% EQUATION OF MOTION
% inertial matrix
[B, U] = inertial_matrix(kin);
% Coriolis matrix
C = coriolis_matrix(q, qD, B);
% gravitational moment vector
G = gravity_matrix(kin);

% equation of dynamics (without tau)
eq = B * qDD + C * qD + kin.P.Fv * qD + kin.P.Fs * sign(qD) + G;

dyn.B = B;
dyn.C = C;
dyn.G = G;

% ENERGY
% kinetic energy
K = simplify(1/2 * qD' * B * qD);
% potential energy
U = U;

dyn.K = K;
dyn.U = U;

end

