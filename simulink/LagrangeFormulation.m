function [sys,x0,str,ts,simStateCompliance] = LagrangeFormulation(t,x,u,flag,g)

    switch flag

      %%%%%%%%%%%%%%%%%%
      % Initialization %
      %%%%%%%%%%%%%%%%%%
      case 0
        [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes;

      %%%%%%%%%%%%%%%
      % Derivatives %
      %%%%%%%%%%%%%%%
      case 1
        sys=mdlDerivatives(t,x,u,g);

      %%%%%%%%%%
      % Update %
      %%%%%%%%%%
      case 2
        sys=mdlUpdate(t,x,u);

      %%%%%%%%%%%
      % Outputs %
      %%%%%%%%%%%
      case 3
        sys=mdlOutputs(t,x,u);

      %%%%%%%%%%%%%%%%%%%%%%%
      % GetTimeOfNextVarHit %
      %%%%%%%%%%%%%%%%%%%%%%%
      case 4
        sys=mdlGetTimeOfNextVarHit(t,x,u);

      %%%%%%%%%%%%%
      % Terminate %
      %%%%%%%%%%%%%
      case 9
        sys=mdlTerminate(t,x,u);

      %%%%%%%%%%%%%%%%%%%%
      % Unexpected flags %
      %%%%%%%%%%%%%%%%%%%%
      otherwise
        DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));

    end

end

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes()

    %
    % call simsizes for a sizes structure, fill it in and convert it to a
    % sizes array.
    %
    % Note that in this example, the values are hard coded.  This is not a
    % recommended practice as the characteristics of the block are typically
    % defined by the S-function parameters.
    %
    sizes = simsizes;

    sizes.NumContStates  = 6;  % x = [q; qD]
    sizes.NumDiscStates  = 0;
    sizes.NumOutputs     = 6;  % y = [qD; q]
    sizes.NumInputs      = 3;  % u = [tau]
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 1;   % at least one sample time is needed

    sys = simsizes(sizes);

    %
    % initialize the initial conditions
    %
    x0  = [0; 0; 0; 0; 0; 0];  % [q; qD]

    %
    % str is always an empty matrix
    %
    str = [];

    %
    % initialize the array of sample times
    %
    ts  = [0 0];

    % Specify the block simStateCompliance. The allowed values are:
    %    'UnknownSimState', < The default setting; warn and assume DefaultSimState
    %    'DefaultSimState', < Same sim state as a built-in block
    %    'HasNoSimState',   < No sim state
    %    'DisallowSimState' < Error out when saving or restoring the model sim state
    simStateCompliance = 'UnknownSimState';

end


%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys=mdlDerivatives(t,x,u,g)
    % x = [q; qD]
    q = x(1:3);
    qD = x(4:6);
    tau = u;

    % B(q)qDD + C(q,qD)qD + F(q)qD + G(q) = tau - J'he
    % all these are just pre computed matrices
    B = double(B_Lagrangian(q));
    C = double(C_Lagrangian(q, qD));
    G = double(G_Lagrangian(q, g));
    J = double(Jacobian(q));
    
    qDD = B \ (tau - C*qD - G);  % - J'*he
    
    sys = [qD; qDD];
end


%
%=============================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=============================================================================
%
function sys=mdlUpdate(t,x,u)

    sys = [];

end


%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)
    % x = [q; qD]
    q = x(1:3);
    qD = x(4:6);

    % output is simply the state
    sys = [qD; q];
end


%
%=============================================================================
% mdlGetTimeOfNextVarHit
% Return the time of the next hit for this block.  Note that the result is
% absolute time.  Note that this function is only used when you specify a
% variable discrete-time sample time [-2 0] in the sample time array in
% mdlInitializeSizes.
%=============================================================================
%
function sys=mdlGetTimeOfNextVarHit(t,x,u)

    sampleTime = 1;
    sys = t + sampleTime;

end


%
%=============================================================================
% mdlTerminate
% Perform any end of simulation tasks.
%=============================================================================
%
function sys=mdlTerminate(t,x,u)

    sys = [];

end

