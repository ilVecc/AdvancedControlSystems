function [C] = coriolis_matrix(q, qD, B)
%CORIOLIS_MATRIX compute the Coriolis matrix of the system given the
%inertial matrix and the joint parameters

dB = cat(3, diff(B, q(1)), diff(B, q(2)), diff(B, q(3)));

C = zeros(3,'sym');
for i=1:3
    for j=1:3
        for k=1:3
            cijk = 1/2 * (dB(i,j,k) + dB(i,k,j) - dB(j,k,i)) * qD(k);
            C(i, j) = C(i, j) + cijk;
        end
    end
end
C = simplify(C);
end

