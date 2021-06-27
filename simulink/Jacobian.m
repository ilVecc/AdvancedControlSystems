function [J] = Jacobian(q)

t1 = q(1);
t2 = q(2);
d3 = q(3);

% wrt frame 0
J = [[ - (3*sin(t1 + t2))/10 - (2*sin(t1))/5, -(3*sin(t1 + t2))/10, 0];
     [   (3*cos(t1 + t2))/10 + (2*cos(t1))/5,  (3*cos(t1 + t2))/10, 0];
     [                                     0,                    0, 1];
     [                                     0,                    0, 0];
     [                                     0,                    0, 0];
     [                                     1,                    1, 0]];
 
end

