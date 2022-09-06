function Xdot = Dyn_Comb_dim(t, X, Var)
    
    %% Input Variables

    % The following quantities are in the S frame
    R = X(1:3);                     % Position 
    nu = X(7:9);                    % Translational Velocity 
    
    x = R(1); y = R(2); z = R(3);
    xdot = nu(1); ydot = nu(2); zdot = nu(3);


    % Angular velcoity is in B frame
    w = X(4:6);                     % Angular Velocity
    
    w1 = w(1); w2 = w(2); w3 = w(3);

    % Parameters
    m1 = Var.m1; m2 = Var.m2;
    D1 = Var.D1; D2 = Var.lstar - D1;

    d13 = sqrt((x + D1)^2 + y^2 + z^2);
    d23 = sqrt((x - D2)^2 + y^2 + z^2);

    K1 = Var.K(1); K2 = Var.K(2); K3 = Var.K(3);
    
    SB = reshape(X(10:18),3,3);     % DCM mapping from B to S frame

    E = SB' * (R + [D1 0 0]');
    E1 = E(1); E2 = E(2); E3 = E(3);

    L = SB' * (R - [D2 0 0]');
    L1 = L(1); L2 = L(2); L3 = L(3);
     
    %% Velocity Dynamics on SE(3)
    
    v = [w; SB'*nu];                % Velocity vector in B frame

    nudotS = [2*ydot + x - m1*(x+D1)/d13^3 - m2*(x-D2)/d23^3;
        y - 2*xdot - m1*y/d13^3 - m2*y/d23^3;
        -m1*z/d13^3 - m2*z/d23^3];
    
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

    Xdot = [Rdot; wdot; nudotS; reshape(SBdot,[],1)]
end