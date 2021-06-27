function [] = simulate(eqn, init_cond)
%SIMULATE run the simulation of the given equations
[ode, Y] = odeToVectorField(eqn);
odeFun = matlabFunction(ode, 'Vars', {'t', 'Y'});
time_span = [0 10];
[t, Q] = ode45(odeFun, time_span, init_cond);
Q = [Q(:,3) Q(:,1) Q(:,5) Q(:,4) Q(:,2) Q(:,6)];
plot(t, Q);
legend({'$t_1(t)$', '$t_2(t)$', '$d_3(t)$', ...
        '$\dot{t}_1(t)$', '$\dot{t}_2(t)$', '$\dot{d}_3(t)$'}, ...
       'Location', 'best', 'Interpreter', 'latex', 'FontSize', 15)
end

