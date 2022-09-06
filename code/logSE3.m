function eta = logSE3(g)
    
    [NB, R] = invSE(g);

    %% thetavec
    
    thetanorm = acos(0.5*trace(NB) - 1);

    if thetanorm == 0
        thetaX = 0;
    else
        thetaX = thetanorm/(2*sin(thetanorm)) * (NB - NB');
    end

    thetavec = invSO3(thetaX);
    
    %% S(theta)
    
    S = Sfunc(thetanorm, thetaX);

    %% beta

    beta = S\R;

    eta = [thetavec; beta];
end