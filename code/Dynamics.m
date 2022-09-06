function Xdot = Dynamics(t,X,Var)

    %% Input Variables

    % The following quantities are in the N frame
    R = X(1:3);                     % Position 
    nu = X(7:9);                    % Translational Velocity
    
    % Angular velcoity is in B frame
    w = X(4:6);                     % Angular Velocity

    OB = reshape(X(10:18),3,3);     % DCM mapping from B to O frame
    NO = NODCM(t, Var);             % DCM mapping from O to N frame
    NB = NO * OB; % B to N; from B to O then from O to N

    % Convert position and attitude to SE(3)
    g = SE3(NB, R); % R is in N frame in SE(3)

    % Planetary positions in N frame
    [R_E, R_L, ~] = SysInert(t, Var);

    %% Velocity dynamics

    v = [w;  NB'*nu]; % v is in body-fixed frame
    
    F_G_N = GravF_Inert(R, R_E, R_L, Var);
    F_G_B = NB' * F_G_N';

    phi = [0;0;0;F_G_B];
    
    vdot = Var.I*(Coadj(v)*Var.I*v + phi); % vdot is in body-fixed frame
    wdot = vdot(1:3); % wdot kept in body-fixed frame
    
    %nudot = F_G_N';
    nudot = NB*vdot(4:6);

    %% Pos/Attitude Dynamics
    
    vvee = se3alg(v);
    gdot = g*vvee;

    [NBdot, Rdot] = invSE3(gdot); % Rdot is in N frame by definition.

    %% Output

    Xdot = [Rdot; wdot; nudot; reshape(NBdot,[],1)];

end