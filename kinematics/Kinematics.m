classdef Kinematics < handle
    %KINEMATICS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        DH(:,4) sym 
        dofs double 
        joint_type cell
        J(6,:) sym
        JA(6,:) sym
        TA sym
        JL(6,:,:) sym
        P MechanicalParameters
    end
    
    properties (SetAccess = private)
        cache_H containers.Map
        idx cell
    end
    
    methods
        function this = Kinematics(DH, joint_type, P)
            %KINEMATICS Store a DH table and prepare for caching
            this.DH = DH;
            this.dofs = size(DH, 1) - 2;
            this.cache_H = containers.Map();
            this.joint_type = joint_type;
            this.idx = [{'b'}, {'0'}, arrayfun(@num2str, 1:this.dofs, 'UniformOutput', 0), {'e'}];
            
            % compute all [i-1 -> i] homogenous matrices from base to ee
            for from=1:(this.dofs + 2)
                key = strcat(this.idx{from}, '_', this.idx{from+1});
                this.cache_H(key) = this.homogeneous(this.DH(from,:));
            end
            
            % compute jacobians
            this.J = this.jacobian_geometric();
            [this.JA, this.TA] = this.jacobian_analytical();
            
            % add gravity w.r.t. frame 0
            P.g0 = trans(this.H('b', 0)', P.gbase);
            this.P = P;
        end
        
        function Hi = H(this, from, to)
            %DIRECT computes direct kinematics given the starting and
            %ending reference frames
            
            if from == to
                Hi = sym(eye(4));
                return
            end
            
            if ~ischar(from)
                from = int2str(from);
            end
            if ~ischar(to)
                to = int2str(to);
            end
            
            key = strcat(from, '_', to);
            if this.cache_H.isKey(key)
                % get from cache
                Hi = this.cache_H(key);
            else
                % compute and cache recursively
                if strcmp(to, 'e')
                    to_prev = this.dofs;
                else
                    to_prev = str2double(to) - 1;
                end
                H_prev = this.H(from, to_prev);
                H_this = this.H(to_prev, to);
                Hi = simplify(H_prev * H_this);
                this.cache_H(key) = Hi;
            end
        end
        
        function Hi = H_num(this, from, to, q, qd)
            %DIRECT_NUMERIC computes direct kinematics given the starting 
            %and ending reference frames and convert it to numeric
            Hi = double(subs(this.numeric(this.H(from, to)), q, qd));
        end
        
        function J_num = J_num(this, q, qd)
            J_num = double(subs(this.numeric(this.J), q, qd));
        end
        
        function JA_num = JA_num(this, q, qd)
            Hb_e = this.H_num('b', 'e', q, qd);
            Phi = rotm2eul(Hb_e(1:3,1:3), 'ZYZ');
            val.phi   = Phi(1);
            val.theta = Phi(2);
            val.psi   = Phi(3);
            JA_num = double(subs(this.numeric(this.JA), q, qd));
            JA_num = double(subs(JA_num, val));
        end
        
        function [q] = inv(this, X)
            % INV inverse kinematics of the robot
            L0 = this.numeric(this.P.L(0));
            L1 = this.numeric(this.P.L(1));
            L2 = this.numeric(this.P.L(2));
            L3 = this.numeric(this.P.L(3));
            wx = X(1);
            wz = X(3) - L0;
            c2 = (wx^2 + wz^2 - L1^2 - L2^2) / (2*L1*L2);
            % we use abs so to ignore precision error leading to undesirable imaginary part
            s2 = abs(sqrt(1 - c2^2));
            t2 = atan2(s2, c2);

            c1 = ((L1 + L2*c2) * wx + L2*s2*wz) / (wx^2 + wz^2);
            s1 = ((L1 + L2*c2) * wz - L2*s2*wx) / (wx^2 + wz^2);
            t1 = atan2(s1, c1);

            d3 = -X(2) - L3;

            q = [t1, t2, d3]';
        end
        
        function [matrix] = numeric(this, matrix)
            matrix = this.P.numeric(matrix);
        end
        
    end
    
    methods (Static)
        
        function [X] = pose(H)
            %H_TO_POSE Extracts the pose from an homogeneous matrix
            R = H(1:3,1:3);
            Phi = rotm2eul(R, 'ZYZ');
            X = [H(1:3,4); Phi'];
        end
        
    end
    
    methods (Access = private)
        
        function [J] = jacobian_geometric(this)
            %JACOBIAN_GEOMETRIC computes the geometric Jacobian given the 
            %homogeneous matrices
            z0 = [0; 0; 1];
            
            H0_e = this.H(0, 'e');
            d0_e = H0_e(1:3,4);

            J = sym(zeros(6, this.dofs));
            for i=1:this.dofs
                % velocity of the i-th joint
                H0_i = this.H(0, i-1);
                R0_i = H0_i(1:3,1:3);
                if strcmp(this.joint_type{i}, 'R')
                    % is rotational joint
                    d0_i = H0_i(1:3,4);
                    J(:,i) = [cross(R0_i * z0, d0_e - d0_i);
                                                  R0_i * z0];
                else
                    % is prismatic joint
                    J(:,i) = [R0_i * z0;
                              zeros(3,1)];
                end
            end
            J = simplify(J);
        end
        
        function [JA, TA] = jacobian_analytical(this)
            %JACOBIAN_ANALYTICAL computes the analytical Jacobian given the 
            %homogeneous matrices
            Phi = [sym('ph', 'real'); 
                   sym('th', 'real'); 
                   sym('ps', 'real')];
            Tphi = [0 -sin(Phi(1)) cos(Phi(1))*sin(Phi(2));
                    0  cos(Phi(1)) sin(Phi(1))*sin(Phi(2));
                    1       0             cos(Phi(2))     ];
            TA = [eye(3)   zeros(3);
                  zeros(3)   Tphi  ];
            JA = TA \ this.J;
        end
        
    end
    
    methods (Static, Access = private)
        
        function [Hi] = homogeneous(DH_row)
            %HOMOGENEOUS Computes the homogeneous transofrmation matrix for the given
            % parameters (z displacement, z rotation, x displacement, x rotation)
            Di = DH_row(1);
            ti = DH_row(2);
            Ai = DH_row(3);
            ai = DH_row(4);
            Hi = [[cos(ti)  -sin(ti)*cos(ai)   sin(ti)*sin(ai)  Ai*cos(ti)]; 
                  [sin(ti)   cos(ti)*cos(ai)  -cos(ti)*sin(ai)  Ai*sin(ti)]; 
                  [   0          sin(ai)           cos(ai)          Di    ]; 
                  [   0             0                 0              1    ]];
        end
        
    end
    
end

