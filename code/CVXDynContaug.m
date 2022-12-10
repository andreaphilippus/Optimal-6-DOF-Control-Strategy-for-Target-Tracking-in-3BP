function dot = CVXDynContaug(t, in, controls, Var)
    % states: 4-by-4 SE(3) state
    % costates: velocity vector
    % controls: control vector

    states = in(1:12); vels = in(13:18);
    
    %% Velocity Dynamics on SE(3)

    SB = reshape(states(1:9),3,3);
    R = states(10:12);
    g = SE3(SB, R); 
    %[SB, R] = invSE3(states);

    w = vels(1:3);
    w1 = w(1); w2 = w(2); w3 = w(3);

    nu = vels(4:6);

    x = R(1); y = R(2); z = R(3);
    xdot = nu(1); ydot = nu(2); zdot = nu(3);

    % Parameters
    m1 = Var.m1; m2 = Var.m2; mu = Var.mu;

    d13 = sqrt((x + mu)^2 + y^2 + z^2);
    d23 = sqrt((x - 1 + mu)^2 + y^2 + z^2);

    K1 = Var.K(1); K2 = Var.K(2); K3 = Var.K(3);

    E = SB' * (R - [mu 0 0]');
    E1 = E(1); E2 = E(2); E3 = E(3);

    L = SB' * (R - [1-mu 0 0]');
    L1 = L(1); L2 = L(2); L3 = L(3);

    %% Velocity Dynamics on SE(3)

    nudotS = [2*ydot + x - (1-mu)*(x+mu)/d13^3 - mu*(x-1+mu)/d23^3;
        y - 2*xdot - (1-mu)*y/d13^3 - mu*y/d23^3;
        - (1-mu)*z/d13^3 - mu*z/d23^3];

    %wdot = [K1 * (3*(1-mu)/d13^5 *E2*E3 + 3*mu/d23^5 *L2*L3 - w2*w3);
    %    K2 * (3*(1-mu)/d13^5 *E1*E3 + 3*mu/d23^5 *L1*L3 - w1*w3);
    %    K3 * (3*(1-mu)/d13^5 *E1*E2 + 3*mu/d23^5 *L1*L2 - w1*w2)];
    
    wdot = [K1 * (-w2*w3);
        K2 * (- w1*w3);
        K3 * (- w1*w2)];
    
    if in(end) ~= 1
        controls = 0;
    end

    veldot = [wdot; nudotS] + controls;

    
    v = [eye(3) zeros(3); zeros(3) SB'] * vels; % Velocities in B frame
    vvee = se3alg(v);

    [SBdot, Rdot] = invSE3(g * vvee);

    statedot = [reshape(SBdot,9,1); Rdot];

    dot = [statedot; veldot; 1];
    
    %{
    % Gravitational force in B frame
    phi = SB'*gravforceN([R(1) R(2) R(3) 0 0 0 vels(4) vels(5)], Var);   
    phi = [0;0;0;phi];

    vdot = [eye(3) zeros(3); zeros(3) SB] * (Var.I\(Coadj(v)*Var.I*v + phi));

    veldot = vdot + controls;
    %}
end