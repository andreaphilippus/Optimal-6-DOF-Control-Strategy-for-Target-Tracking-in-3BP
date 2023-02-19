function Xdot = SE3CR3BPEoM(t, X, I, mu, u)
    
    g_BS = reshape(X(1:16),4,4);
    v_BS = X(17:22);

    %% Synodic Frame relative to Inertial Frame Info
    v_S = [0 0 1 0 0 0]';
    g_S = t * vee(v_S);
    g_S = [cos(t) sin(t) 0 0;
        -sin(t) cos(t) 0 0;
        0 0 1 0;
        0 0 0 1];
    [R_S, ~] = invSE3(g_S);
    
    %% v_B
    g_B = g_S * g_BS;
    v_B = v_BS + AdSE3(g_BS)\v_S;

    %% I_S

    I_S = AdSE3(g_BS)'\I/AdSE3(g_BS);
    [SB, R_BS] = invSE3(g_BS);
    J = SB*I(1:3,1:3)*SB';
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
    %Fg = -m*(1-mu)/norm(d)^3 * d...
    %    -m*mu/norm(r)^3 * r;
     
    Fg = R_S*m*(norm(d)^3/(1-mu) * d + norm(r)^3/mu * r);
    %Fg = -m*(1-mu)/norm(d)^3 * (eye(3) + 3/m/norm(d)^2*(J+1/2*(trace(J) - 5*(d/norm(d))' * J * d/norm(d))*eye(3))) *d...
    %    -m*mu/norm(r)^3 * (eye(3) + 3/m/norm(r)^2*(J+1/2*(trace(J) - 5*(r/norm(r))' * J * r/norm(r))*eye(3))) * r;
    
    %% Paper
    % Eq. (53~54)
    Ctau = 3*((1-mu)/norm(d)^5 * CrossProd(d)*J*d + mu/norm(r)^5 * CrossProd(r)*J*r);
    Cf = R_S'*m*((1-mu)/norm(d)^3* d + mu/norm(r)^3 * r);

    %%
    %PHI = [Mg; Fg];

    %% SE(3) Kinematics
    w_S = [0 0 1]';
    w = v_BS(1:3);
    nu = v_BS(4:6);
    
    %g_BS_dot = g_S\vee(-v_S + v_B) * g_B;
    %g_BS_dot = g_BS * vee(v_BS);
    %g_BS_dot = g_BS*vee(v_B) - vee(v_S) * g_BS;
    g_BS_dot = [CrossProd(w)*SB nu; 0 0 0 0];
    
    
    %v_BS_dot = I\(Coadj(v_BS)*I*v_BS + PHI + u);

    wdot =  - J\CrossProd(w)*J*w + J\Mg -CrossProd(w_S)*w;
    nudot = -2*CrossProd(w_S)*nu - CrossProd(w_S)^2*R_BS + 1/m * Cf;

    %% Paper
    % Eq. (52)
    %wdot = -CrossProd(w_S)*w - J\CrossProd(w)*J*w; % + J\Ctau;
    %nudot = -2*CrossProd(w_S)*nu - CrossProd(w_S)*CrossProd(w_S)*R_BS + 1/m*Cf;

    v_BS_dot = [wdot; nudot] + u;

    %% Output
    Xdot = [reshape(g_BS_dot,[],1); v_BS_dot];
end