function PlotSys(Var)
   
    [XL,YL,ZL] = sphere;
    % Plot Luna
    rL = Var.R2 / Var.lstar; % Nondimensionalize the radius
    Xp = XL * rL; Yp = YL * rL; Zp = ZL * rL;
    Luna = surf(Xp + (1-Var.mu), Yp, Zp); hold on
    Luna.EdgeColor = 'none';
    set(Luna, 'FaceColor', [0.7 0.7 0.7])
    
    % Plot Earth
    rE = Var.R1 / Var.lstar; % Nondimensionalize the radius
    Xp = XL * rE; Yp = YL * rE; Zp = ZL * rE;
    Earth = surf(Xp - Var.mu, Yp, Zp); 
    Earth.EdgeColor = 'none';
    set(Earth, 'FaceColor', [91, 207, 244] / 255)
    
    % Plot L2
    plot(Var.L2x, 0, 'r*')
end