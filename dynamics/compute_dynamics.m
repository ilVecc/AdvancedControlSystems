function [tau_val] = compute_dynamics(P, eq, q, qD, qDD, he, qd, qDd, qDDd, hed, tau)
%COMPUTE_DYNAMICS calculate numeric tau for Lagrange and Newton-Euler 
%formulations.

eq = simplify(subs(P.numeric(eq), [q; qD; qDD], [qd; qDd; qDDd]));
eq = simplify(subs(eq, he, hed));
tau_val = solve(eq, tau);
tau_val = structfun(@double, tau_val, 'uniformoutput', 0);

end

