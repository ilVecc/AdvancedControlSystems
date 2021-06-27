function [q, qD, qDD, dynq, dynqD, dynqDD, tau] = params_joint()
%PARAMS_JOINT Summary of this function goes here

% KINEMATICS
% joint position parameters
q = [sym('t1', 'real'); 
     sym('t2', 'real'); 
     sym('d3', 'real')];

% joint velocity parameters
qD = [sym('Dt1', 'real');
      sym('Dt2', 'real');
      sym('Dd3', 'real')];

% joint acceleration parameters
qDD = [sym('DDt1', 'real');
       sym('DDt2', 'real');
       sym('DDd3', 'real')];


% DYNAMICS
% joint position parameters
t = sym('t', 'real');
dynq = [symfun('t1(t)', t);
        symfun('t2(t)', t);
        symfun('d3(t)', t)];

% joint velocity parameters
dynqD = diff(dynq);

% joint acceleration parameters
dynqDD = diff(dynqD);

% join imposed torques
tau = [sym('Tt1', 'real');
       sym('Tt2', 'real');
       sym('Td3', 'real')];

end

