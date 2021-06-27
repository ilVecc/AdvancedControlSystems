classdef MechanicalParameters
    %MECHANICALPARAMETERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dofs int32
        Fv(:,:) sym
        Fs(:,:) sym
        gbase sym
        g0 sym
        Im(1,:) sym
        Kr(1,:) sym
        zm(3,:) sym
    end
    
    properties (SetAccess = private)
        cache_m(1,:) sym
        cache_L(1,:) sym
        cache_I(3,3,:) sym
        cache_pCoM(3,:) sym
        cache_Iaug(3,3,:) sym
    end
    
    methods
        
        function this = MechanicalParameters()
            
            this.dofs = 3;
            
            % KINEMATICS
            % mass [kg]
            this.cache_m = [sym('L0m', 'positive')
                            sym('L1m', 'positive')
                            sym('L2m', 'positive')
                            sym('L3m', 'positive')];
            this.cache_L = [sym('L0h', 'positive')   %0.4;   % [m]
                            sym('L1h', 'positive')   %0.4;   % [m]
                            sym('L2h', 'positive')   %0.3;   % [m]
                            sym('L3h', 'positive')]; %0.4;   % [m]
            
            %base link
            L0m = this.m(0);
            L0h = this.L(0);
            L0r = sym('L0r', 'positive'); %0.025; % [m]
            % link 1
            L1m = this.m(1);
            L1h = this.L(1);
            L1r = sym('L1r', 'positive'); %0.020; % [m]
            % link 2
            L2m = this.m(2);
            L2h = this.L(2);
            L2r = sym('L2r', 'positive'); %0.015; % [m]
            % link 3
            L3m = this.m(3);
            L3h = this.L(3);
            L3w = sym('L3w', 'positive'); %0.025; % [m]
            L3d = sym('L3d', 'positive'); %0.025; % [m]
            
            % DYNAMICS
            % centers of mass w.r.t. local frame
            this.cache_pCoM = cat(2, [-L1h/2;     0;      0], ...
                                     [-L2h/2;     0;      0], ...
                                     [     0;     0; -L3h/2]);
            % inertial tensors w.r.t. local frame
            L1I = this.I_cylin(L1m, L1r, L1h)      + this.steiner(L1m, this.pCoM(1));
            L2I = this.I_cylin(L2m, L2r, L2h)      + this.steiner(L2m, this.pCoM(2));
            L3I = this.I_paral(L3m, L3w, L3d, L3h) + this.steiner(L3m, this.pCoM(3));
            this.cache_I = cat(this.dofs, L1I, L2I, L3I);
            % viscous and static frictions of each link
            this.Fv = zeros(this.dofs, 'sym'); %diag(sym('Fv', [this.dofs,1], 'positive'));
            this.Fs = zeros(this.dofs, 'sym'); %diag(sym('Fv', [this.dofs,1], 'positive'));
            % gravity vector w.r.t frame 0
            this.gbase = [0; 0; -sym('g', 'real')];

            % MOTORS
            % inertial value of i-th motor over the z-axis (assuming symmetry of motor)
            this.Im = zeros(1, this.dofs, 'sym'); % = sym('Im', [1, this.dofs], 'positive');
            % gear rateos
            this.Kr = zeros(1, this.dofs, 'sym'); % = sym('Kr', [1, this.dofs], 'positive');
            % motor rotation axis w.r.t. its link (motor i is attached to link i-1)
            this.zm = sym(repmat([0; 0; 1], [1, this.dofs])); % = sym('zm', [3, this.dofs], 'real');

            % AUGMENTED LINK
            % motors have no inertia
            this.cache_Iaug = cat(this.dofs, L1I, L2I, L3I);
            
        end
        
        function mi = m(this,idx)
            mi = this.cache_m(idx+1);
        end
        
        function Li = L(this,idx)
            Li = this.cache_L(idx+1);
        end
        
        function pCoMi = pCoM(this,idx)
            pCoMi = this.cache_pCoM(:,idx);
        end
        
        function Ii = I(this,idx)
            Ii = this.cache_I(:,:,idx);
        end
        
        function Iaugi = Iaug(this,idx)
            Iaugi = this.cache_Iaug(:,:,idx);
        end

    end
    
    methods (Static)
        
        function matrix = numeric(matrix)
            % density (aluminum)
            density = 2710; % [kg/m^3]
            %base link
            val.L0h = 0.4;   % [m]
            val.L0r = 0.025; % [m]
            L0V = pi * val.L0r^2 * val.L0h; % [m^3]
            val.L0m = density * L0V;        % [kg]
            % link 1
            val.L1h = 0.4;   % [m]
            val.L1r = 0.020; % [m]
            L1V = pi * val.L1r^2 * val.L1h; % [m^3]
            val.L1m = density * L1V;        % [kg]
            % link 2
            val.L2h = 0.3;   % [m]
            val.L2r = 0.015; % [m]
            L2V = pi * val.L2r^2 * val.L2h; % [m^3]
            val.L2m = density * L2V;        % [kg]
            % link 3
            val.L3h = 0.4;   % [m]
            val.L3w = 0.025; % [m]
            val.L3d = 0.025; % [m]
            L3V = val.L3w * val.L3d * val.L3h; % [m^3]
            val.L3m = density * L3V;   % [kg]
            % gravity
            val.g = 9.81;
            
            matrix = subs(matrix, val);
            
            if ~hasSymType(matrix, 'variable')
                matrix = double(matrix);
            end
        end
        
    end
    
    methods (Static, Access = private)
        
        function [I] = I_cylin(m, R, h) % h lies on the x-axis
            %I_CYLD Creates the inertial tensor of a cylinder with reference frame in 
            % the CoM and its x-axis aligned with the main axis of inertia 
            I = m * [[1/2*R^2          0                  0        ]; 
                     [   0    1/12*(3*R^2 + h^2)          0        ]; 
                     [   0             0         1/12*(3*R^2 + h^2)]];
        end
        
        function [I] = I_paral(m, x, y, z)
            %I_PRISM Creates the inertial tensor of a parallelepiped prism with 
            % reference frame in the CoM and its x-axis aligned with the main axis 
            % of inertia 
            I = m * [[1/12*(y^2 + z^2)        0                0        ]; 
                     [       0         1/12*(x^2 + z^2)        0        ]; 
                     [       0                0         1/12*(x^2 + y^2)]];
        end
        
        function [I] = steiner(m, r)
            %STEINER Gives the inertial contribution on a rigid body due to the 
            % displacement of its reference frame from the CoM
            I = m * (r'*r*eye(3) - r*r');
        end
        
    end
end

