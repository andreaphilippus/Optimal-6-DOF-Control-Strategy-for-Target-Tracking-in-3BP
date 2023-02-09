function eta = logSE3(g0, g)
    % g0 : reference state (to be targeted)
    % g : active state (to be controlled)

    [C, R] = invSE(g);
    [C0, R0] = invSE(g0);
    
    Ce = C0'*C;
    Re = C0'*(R - R0);
    invg0g = SE3(Ce, Re);

    %% thetavec
    
    thetanorm = acos(0.5*trace(Ce) - 1);

    if thetanorm == 0
        thetaX = zeros(3,3);
    else
        thetaX = thetanorm/(2*sin(thetanorm)) * (Ce - Ce');
    end

    thetavec = invSO3(thetaX);
    
    %% S(theta)
    
    S = Sfunc(thetanorm, thetaX);

    %% beta

    beta = S\Re;

    eta = [thetavec; beta];
end