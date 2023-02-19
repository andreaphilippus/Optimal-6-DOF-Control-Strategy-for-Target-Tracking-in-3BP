function Xdot = SE3CR3BP(t, X, I, mu, u)
    
    g = reshape(X(1:16), 4, 4);
    v = X(17:end);

    w = v(1:3);
    nu = v(4:6);

    [SB, r] = invSE3(g); % R: S relative to B, [SB], r: pos in Synodic frame
    
    %% Synodic frame
    wS = [0 0 1]';

    IS = SB*I(1:3,1:3)*SB';
    m = I(4,4);

    %% Gravity
    r13 = r - [-mu 0 0]';
    r23 = r - [1-mu 0 0]';

    tauS = 3*((1-mu)/norm(r13)^5 * CrossProd(r13)*IS*r13 +mu/norm(r23)^5 * CrossProd(r23)*IS*r23);

    fS = -m*((1-mu)/norm(r13)^3 * r13 + mu/norm(r23)^3 * r23);
    

    %% 
    wdot = -CrossProd(wS)*w - IS\CrossProd(w)*IS*w + IS\tauS;
    nudot = -2*CrossProd(wS)*nu - CrossProd(wS)*CrossProd(wS)*r + 1/m*fS;

    vdot = [wdot; nudot] + u;

    gdot = [-CrossProd(w)*SB nu; 0 0 0 0];
    
    Xdot = [reshape(gdot,[],1); vdot];
end