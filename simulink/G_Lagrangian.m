function [G] = G_Lagrangian(q, g)

t1 = q(1);
t2 = q(2);
d3 = q(3);

G = [ -(271*g*pi*((9*cos(t1 + t2))/10 + 2*cos(t1)))/1250;
                       -(65853*pi*g*cos(t1 + t2))/800000;
                                                       0];
 
end

