function eta = logSE3(g)
    % g : usually error in this case

    [Ce, Re] = invSE3(g);

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