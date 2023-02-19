function [statedot, veldot] = SE3CasADiDynamics(states, vels, controls, Var)
    
    SB = reshape(states(1:9), 3, 3);
    r = states(10:12);

    v = vels;

    w = SB*v(1:3);
    nu = v(4:6);

    mu = Var.mu;
    I = Var.I;

    %% Synodic frame
    wS = [0 0 1]';

    IS = SB*I(1:3,1:3)*SB';
    m = I(4,4);

    %% Gravity
    r13 = r - [-mu 0 0]';
    r23 = r - [1-mu 0 0]';

    tauS = 3*((1-mu)/norm(r13)^5 * CrossProd(r13)*IS*r13 + mu/norm(r23)^5 * CrossProd(r23)*IS*r23);

    fS = -m*((1-mu)/norm(r13)^3 * r13 + mu/norm(r23)^3 * r23);

    %% SE(3) Kinematics
    statedot = [reshape(-CrossProd(w)*SB,9,1); nu];

    %% SE(3) Kinetics
    wdot = -CrossProd(wS)*w - IS\CrossProd(w)*IS*w + IS\(tauS +controls(1:3));
    nudot = -2*CrossProd(wS)*nu - CrossProd(wS)*CrossProd(wS)*r + 1/m*(fS + controls(4:6));

    veldot = [SB'*wdot; nudot];
end