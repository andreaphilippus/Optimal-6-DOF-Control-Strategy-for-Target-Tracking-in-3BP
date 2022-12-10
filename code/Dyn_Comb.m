function Xdot = Dyn_Comb(t, X, Var, pert)
    
    %% Input Variables

    % The following quantities are in the S frame
    R = X(1:3);                     % Position 
    nu = X(7:9);                    % Translational Velocity 
    
    x = R(1); y = R(2); z = R(3);
    xdot = nu(1); ydot = nu(2); zdot = nu(3);


    % Angular velcoity is in B frame
    w = X(4:6);                     % Angular Velocity
    
    w1 = w(1); w2 = w(2); w3 = w(3);

    if nargin == 3 % When perturbation not added:
        pert = [0 0]';
    end

    % Parameters
    m1 = Var.m1; m2 = Var.m2; mu = Var.mu;

    d13 = sqrt((x + mu)^2 + y^2 + z^2);
    d23 = sqrt((x - 1 + mu)^2 + y^2 + z^2);

    K1 = Var.K(1); K2 = Var.K(2); K3 = Var.K(3);
    
    SB = reshape(X(10:18),3,3);     % DCM mapping from B to S frame

    E = SB' * (R - [mu 0 0]');
    E1 = E(1); E2 = E(2); E3 = E(3);

    L = SB' * (R - [1-mu 0 0]');
    L1 = L(1); L2 = L(2); L3 = L(3);
     
    %% Velocity Dynamics on SE(3)
    
    v = [w; SB'*nu];                % Velocity vector in B frame

    
    %{

    F_B = SB'*gravforceN(X, Var);   % Gravitational force in B frame

    % To be updated. For now, it is neglected.
    FJ2_B = [0;0;0];                % J2 Perturbation in B frame
    M_B = [0;0;0];                  % Gravity gradient in B frame

    phi = [M_B; F_B + FJ2_B];       % Gravitational effects in B frame
    
    vdot = Var.I\(Coadj(v)*Var.I*v + phi);
    wdot = vdot(1:3);
    nudot = SB*vdot(4:6);

    %}

    nudotS = [2*ydot + x - (1-mu)*(x+mu)/d13^3 - mu*(x-1+mu)/d23^3;
        y - 2*xdot - (1-mu)*y/d13^3 - mu*y/d23^3;
        - (1-mu)*z/d13^3 - mu*z/d23^3] + pert(2)*randn(3,1);
    
    %nudotB = SB' * nudotS;

    %wdot = [K1 * (3*(1-mu)/d13^5 *E2*E3 + 3*mu/d23^5 *L2*L3 - w2*w3);
    %    K2 * (3*(1-mu)/d13^5 *E1*E3 + 3*mu/d23^5 *L1*L3 - w1*w3);
    %    K3 * (3*(1-mu)/d13^5 *E1*E2 + 3*mu/d23^5 *L1*L2 - w1*w2)];

    wdot = [K1 * (-w2*w3);
        K2 * (- w1*w3);
        K3 * (- w1*w2)];
 

    %% Position, Attitude Dynamics on SE(3)

    g = SE3(SB, R);
    vvee = se3alg(v);

    gdot = g*vvee;
    [SBdot, Rdot] = invSE3(gdot); % Note that Rdot here is in S frame
    
    %% Output

    Xdot = [Rdot; wdot; nudotS; reshape(SBdot,[],1)];
end