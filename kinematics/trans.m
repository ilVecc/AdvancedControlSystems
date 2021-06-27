function [P] = trans(H, p)
%TRANSFORM Apply an homogeneous transformation on a 3D point/vector
P = H * [p; 1];
P = P(1:3);
end

