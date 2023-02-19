function Xdot = SE3CR3BPEoM_inertial(t, X, I, mu, u)
    
    g = reshape(X(1:16),4,4); % g = [ [NB] Nr; 0 0 0 1]
    v = X(17:22);   % v = Bv

    v_N = AdSE3(g) * v;

    %% I_N

    I_N = AdSE3(g')\I/AdSE3(g);
    [NB, r] = invSE3(g);
    J_N = NB*I(1:3,1:3)*NB';
    m = I(4,4);

    %% Dynamics
    NS = [cos(t) -sin(t) 0;
        sin(t) cos(t) 0;
        0 0 1];
    p_1 = NS*[-mu 0 0]';
    p_2 = NS*[1-mu 0 0]';
    mu1 = 1-mu;
    mu2 = mu;
    
    d_N = r - p_1;
    r_N = r - p_2;

    tau_N = 3 * (mu1/norm(d_N)^5 * CrossProd(d_N)*J_N*d_N + ...
        mu2/norm(r_N)^5 * CrossProd(r_N)*J_N*r_N);

    f_N = -m * (mu1/norm(d_N)^3 * d_N + mu2/norm(r_N)^3 * r_N);

    PHI = [tau_N; f_N];

    %% SE(3) Kinematics
    %{
    g_BS_dot = g_S\vee(-v_S + v_B) * g_B;
    %g_BS_dot = g_BS * vee(v_BS);
    w_S = [0 0 1]';
    w = v_BS(1:3);
    nu = v_BS(4:6);
    
    %v_BS_dot = I\(Coadj(v_BS)*I*v_BS + PHI + u);

    wdot =  - J\CrossProd(w)*J*w + J\Mg -CrossProd(w_S)*w;
    nudot = -2*CrossProd(w_S)*nu - CrossProd(w_S)^2*R_BS + 1/m * Fg;

    %% Paper
    % Eq. (52) 
    %
    wdot = -CrossProd(w_S)*w - J\CrossProd(w)*J*w; % + J\Ctau;
    nudot = -2*CrossProd(w_S)*nu - CrossProd(w_S)*CrossProd(w_S)*R_BS + 1/m*Cf;

    v_BS_dot = [wdot; nudot] + u;
    %}
    
    g_dot = g*vee(v);
    v_dot_N = I_N\(Coadj(v_N)*I_N*v_N + PHI);
    v_dot = AdSE3(g)\v_dot_N;

    %% Output
    Xdot = [reshape(g_dot,[],1); v_dot];
end