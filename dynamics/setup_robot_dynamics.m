function [] = setup_robot_dynamics(robot, kin)
%SETUP_ROBOT 

% init robot
robot.Gravity           = kin.P.g0;
IL1_1 = kin.P.I(1);
robot.Bodies{1}.Mass    = kin.P.L1.m;
robot.Bodies{1}.Inertia = [IL1_1(1,1) IL1_1(2,2) IL1_1(3,3) ...
                           IL1_1(2,3) IL1_1(1,3) IL1_1(1,2)];
IL2_2 = kin.P.I(2);
robot.Bodies{2}.Mass    = kin.P.L2.m;
robot.Bodies{2}.Inertia = [IL2_2(1,1) IL2_2(2,2) IL2_2(3,3) ...
                           IL2_2(2,3) IL2_2(1,3) IL2_2(1,2)];                    
IL3_3 = kin.P.I(3);
robot.Bodies{3}.Mass    = kin.P.L3.m;
robot.Bodies{3}.Inertia = [IL3_3(1,1) IL3_3(2,2) IL3_3(3,3) ...
                           IL3_3(2,3) IL3_3(1,3) IL3_3(1,2)];

end

