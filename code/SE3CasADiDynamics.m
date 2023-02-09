function [statedot, veldot] = SE3CasADiDynamics(states, vels, controls, Var)
    
    BS = reshape(states(1:9), 3,3);
    R_BS = states(10:12);

    g_BS = SE3(BS, R_BS);
    v_BS = vels;

    mu = Var.mu;

    %% Synodic Frame relative to Inertial Frame Info
    v_S = [0 0 1 0 0 0]';

    %% v_B
    v_B = v_BS + AdSE3(g_BS)\v_S;

    %% I_S
    I = Var.I;

    I_S = AdSE3(g_BS)\I*inv(AdSE3(g_BS));
    
    J = BS*I(1:3,1:3)*BS';
    m = I(4,4);

    %% Relative Pose to the Primary Bodies
    
    d = R_BS - [-mu 0 0]';
    r = R_BS - [1-mu 0 0]';

    %% Dynamics using SE(3)

    % Gravity Gradient Torque
    %Mg = 3*(1-mu)/norm(d)^5 * CrossProd(d)*J*d ...
    %    + 3*mu/norm(r)^5 * CrossProd(r)*J*r;
    Mg = [0 0 0]';

    % Gravitational Force
    Fg = -m*(1-mu)/norm(d)^3 * d...
        -m*mu/norm(r)^3 * r;
        
    %Fg = -m*(1-mu)/norm(d)^3 * (eye(3) + 3/m/norm(d)^2*(J+1/2*(trace(J) - 5*(d/norm(d))' * J * d/norm(d))*eye(3))) *d...
    %    -m*mu/norm(r)^3 * (eye(3) + 3/m/norm(r)^2*(J+1/2*(trace(J) - 5*(r/norm(r))' * J * r/norm(r))*eye(3))) * r;

    PHI = [Mg; Fg];

    %% SE(3) Kinematics
    
    g_BS_dot = g_BS * vee(v_B) - vee(v_S) * g_BS;
    statedot = [reshape(g_BS_dot(1:3,1:3),9,1); g_BS_dot(1:3,4)];

    %g_BS_dot = g_BS * vee(v_BS);
    w_S = [0 0 1]';
    w = v_BS(1:3);
    nu = v_BS(4:6);
    
    %v_BS_dot = I\(Coadj(v_BS)*I*v_BS + PHI + u);

    wdot = -CrossProd(w_S)*w - J\CrossProd(w)*J*w + J\Mg;
    nudot = -2*CrossProd(w_S)*nu - CrossProd(w_S)^2*R_BS + 1/m * Fg;

    veldot = [wdot; nudot] + controls;
end