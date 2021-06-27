%% IMPORT
addpath(genpath('.'));

% !!! doesn't follow the DH convention (prismatic frame is at CoM)
robot = importrobot('RRP_2.urdf');
robot.DataFormat = 'col';  % d vectors are given as rows

showdetails(robot)
%show(robot, homeConfiguration(robot));


% mechanical parameters
P = MechanicalParameters();
% numerical evaluation variables and functions for the dynamic model
[q, qD, qDD, dynq, dynqD, dynqDD, tau] = params_joint();



%% KINEMATICS DIRECT
% Denavit-Hartenberger table (d, T, a, A)
DH = [[    P.L(0)       0      0      pi/2 ];
      [      0         q(1)  P.L(1)    0   ];
      [      0         q(2)  P.L(2)    0   ];
      [ P.L(3) + q(3)   0      0       0   ];
      [      0          0      0       0   ]];
joint_type = {'R', 'R', 'P'};

kin = Kinematics(DH, joint_type, P);
P.g0 = kin.P.g0;

clear joint_type


%% KINEMATICS INVERSE
% robot has only 3 DoF, and actually orientation is not needed
qd = kin.inv([0.5; -0.2; 0.3; 0; 0; 0])
show(robot, qd);


%% GEOMETRIC JACOBIAN
J = kin.J;
% here we obtain the Jacobian from frame 0 to ee frame, but the robotic 
% systems toolbox computes it from base frame, so we introduce the necessary
% transformation
Hb_0 = kin.H_num('b', 0, [], []);
Rb_0 = Hb_0(1:3,1:3);
TJb_0 = [    Rb_0  zeros(3); 
         zeros(3)      Rb_0];

clear Hb_0 Rb_0


%% ANALYTICAL JACOBIAN
JA = kin.JA;



%% TEST KINEMATICS
test = true;
if test
    qd = [pi/2; pi/2; -0.1]

    % MANUAL DIRECT obtain homogeneous matrix for specific configuration
    H_manual = kin.H_num('b', 'e', q, qd)
    % AUTO DIRECT computed from ee to base
    H_auto = getTransform(robot, qd, 'ee', 'base_link')

    % MANUAL INVERSE
    x = kin.pose(H_manual)
    qd_manual = kin.inv(x);
    % AUTO INVERSE
    ik = inverseKinematics('RigidBodyTree', robot);
    weights = [1 1 1 1 1 1];  % for the tolerances of the errors
    qd_auto = ik('ee', H_manual, weights, robot.homeConfiguration);
    

    % MANUAL GEOMETRIC JACOBIAN
    % obtain Jacobian matrix for specific configuration
    % computed from frame 0 to ee frame
    J_manual = kin.J_num(q, qd)
    J_manual_base = TJb_0 * J_manual
    % AUTO GEOMETRIC JACOBIAN
    % computed from the base frame to ee frame
    % (be aware that here linear and angular velocities are swapped, so 
    % angular is above and linear below, thus the second line)
    J_auto_base = geometricJacobian(robot, qd, 'ee');
    J_auto_base = [J_auto_base(4:6,:); J_auto_base(1:3,:)]


    % MANUAL ANALYTICAL JACOBIAN
    JA_manual = kin.JA_num(q, qd);
    % AUTO ANALYTICAL JACOBIAN
    % i don't think this exists
    
    clear qd weights x ik
    clear H_manual J_manual J_manual_base qd_manual JA_manual
    clear H_auto J_auto_base qd_auto
end



%% LAGRANGE FORMULATION
% generalized forces on the ee (forces [position] and torques [orientation])
syms f1 f2 f3 mu1 mu2 mu3 real
he = [f1; f2; f3; mu1; mu2; mu3];
clear f1 f2 f3 mu1 mu2 mu3

[eq, dyn] = lagrange(q, qD, qDD, kin);
eq_lag = (eq == tau - kin.J'*he);
clear eq

%% RECURSIVE NEWTON-EULER
eq = rne(qD, qDD, kin);
eq_rne = (eq == tau);
clear eq

%% TEST DYNAMICS (ROTTO)
qd = [-pi/2; pi/2; -0.1];
qDd = [0; 0; 0];
qDDd = [0; 0; 0];
hed = [0; 0; 0; 0; 1; 0];

% manual
tau_lag_manual = compute_dynamics(P, eq_lag, q, qD, qDD, he, qd, qDd, qDDd, hed, tau)
tau_rne_manual = compute_dynamics(P, eq_rne, q, qD, qDD, he, qd, qDd, qDDd, hed, tau)

fext = externalForce(robot,'ee',TJb_0'*hed,qd);
tau_auto = inverseDynamics(robot, qd, qDd, qDDd, fext)

%clear qd qDd qDDd hed i



%% DYNAMIC MODEL IN OPERATIONAL SPACE (Lagrange)
[eq_lag_op, dyn_op] = lagrange_operational(q, qD, dynq, dynqD, tau, kin, dyn);



%% SIMULATION  (ignore, go to simulink)
test = false;
if test
    % insert gravity, mass, inertia in the robot representation
    setup_robot_dynamics(robot, kin);

    % imposed torque
    taud = gravityTorque(robot, homeConfiguration(robot));

    % dynamic substitution (symbolic variables -> symbolic functions)
    eq = subs(eq, q, dynq);
    eq = subs(eq, qD, dynqD);
    eq = subs(eq, qDD, dynqDD);
    eq = subs(eq, tau, taud);

    % compute simulation
    q.d = homeConfiguration(robot);
    q.D.d = [0 0 0]; % which means it's not initially moving
    init_cond = [
        q.d(2);  
        q.D.d(2);
        q.d(1); 
        q.D.d(1);
        q.d(3); 
        q.D.d(3);
    ];

    simulate(eq, init_cond);
end




