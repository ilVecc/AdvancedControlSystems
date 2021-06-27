function [g] = gravity_matrix(kin)
%GRAVITY_MATRIX computes the gravity contributions to the dynamics

g = sym(zeros([kin.dofs, 1]));
for i=1:kin.dofs
    % gravity contributions
    gj = sym(zeros([3, 1]));
    for j=1:kin.dofs
        % gravity contributions for each link
        JLjPi = kin.JL(1:3,i,j);
        gj = gj + kin.P.m(i) * JLjPi;
    end
    g(i) = simplify(kin.P.g0' * gj);
end

end