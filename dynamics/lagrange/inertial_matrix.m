function [B,U] = inertial_matrix(kin)
%INERTIAL_MATRIX compute partial Jacobians and then the inertial matrix 
%given the mechanical/joint parameters and the homogeneous matrices

% compute partial jacobians
pL = sym(zeros(3, kin.dofs));
for i=1:kin.dofs
    % CoM w.r.t. local frame
    pLi_i = kin.P.pCoM(i);
    % CoM w.r.t. base frame
    pL(:,i) = simplify(trans(kin.H(0, i), pLi_i));
end
kin.JL = jacobian_partial(kin, pL);

% compute inertial tensors
IL = sym(zeros(3, 3, kin.dofs));
for i=1:kin.dofs
    H0_i = kin.H(0, i);
    R0_i = H0_i(1:3,1:3);
    % inertial tensor w.r.t. local frame
    ILi_i = kin.P.I(i);
    % inertial tensor w.r.t. base frame
    IL(:,:,i) = simplify(R0_i * ILi_i * R0_i');
end


% KINETIC AND POTENTIAL ENERGY
% compute inertial matrix and potential scalar
B = sym(zeros(3));
U = sym(zeros(1));
for i=1:kin.dofs
    m = kin.P.m(i);
    JLiP = kin.JL(1:3,:,i);
    JLiO = kin.JL(4:6,:,i);
    ILi = IL(:,:,i);
    % local ineratial matrix
    B = B + simplify(m * (JLiP' * JLiP) + (JLiO' * ILi * JLiO));
    % local potential energy
    U = U + m * kin.P.g0' * pL(:, i);
end
B = simplify(B);
U = simplify(U);


end

