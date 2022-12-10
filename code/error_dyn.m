function Out = error_dyn(t, input, g_ref, v_ref, Var)

    % input: to be integrated
    % {}_ref: FIXED REFERENCE STATES

    %input: eta (6x1) + g (16x1) + v (6x1)
    % -> 28x1

    eta = input(1:6);
    g = reshape(input(7:22),4,4);
    v = input(23:end);

    global ti
    t_index = ti;

    mu = Var.mu;
    %t_index = find(t_ref == t) % Index of the timespan for REFERENCE STATES. 0 < t_index <= size(v_l)
    

    g_l = g_ref(:,:,t_index);
    v_l = v_ref(:,t_index);

    K1 = Var.K(1); K2 = Var.K(2); K3 = Var.K(3);

    
    

    [C, R] = invSE3(g);
    [C_l, R_l] = invSE3(g_l);

    %X = SE3(C'*C_l, R-R_l);
    invX = SE3(C_l'*C, R_l-R);

    x = R(1); y = R(2); z = R(3);
    x_l = R_l(1); y_l = R_l(2); z_l = R(3);

    d13 = sqrt((x + mu)^2 + y^2 + z^2);
    d23 = sqrt((x - 1 + mu)^2 + y^2 + z^2);

    d13_l = sqrt((x_l + mu)^2 + y_l^2 + z_l^2);
    d23_l = sqrt((x_l - 1 + mu)^2 + y_l^2 + z_l^2);

    w1 = v(1); w2 = v(2); w3 = v(3);
    xdot = v(4); ydot = v(5); zdot = v(6);

    w1l = v_l(1); w2l = v_l(2); w3l = v_l(3);
    xdot_l = v_l(4); ydot_l = v_l(5); zdot_l = v_l(6);

    wdot = [K1 * (-w2*w3);
        K2 * (- w1*w3);
        K3 * (- w1*w2)];

    nudot = [2*ydot + x - (1-mu)*(x+mu)/d13^3 - mu*(x-1+mu)/d23^3;
        y - 2*xdot - (1-mu)*y/d13^3 - mu*y/d23^3;
        - (1-mu)*z/d13^3 - mu*z/d23^3]+ randn(3,1);     % Perturbation
    
    vdot = [wdot; nudot];

    wdot_l = [K1 * (-w2l*w3l);
        K2 * (- w1l*w3l);
        K3 * (- w1l*w2l)];

    nudot_l = [2*ydot_l + x_l - (1-mu)*(x_l+mu)/d13_l^3 - mu*(x_l-1+mu)/d23_l^3;
        y_l - 2*xdot_l - (1-mu)*y_l/d13_l^3 - mu*y_l/d23_l^3;
        - (1-mu)*z_l/d13_l^3 - mu*z_l/d23_l^3]; % + pert(2)*randn(3,1); % Perturbation may not be considered for reference state
    
    vdot_l = [wdot_l; nudot_l];

    %v_e = v - capAd(invX, v_l);

    etadot = Gfunc(eta) * v_e;

    %v_edot = Coadj(v) * v + 
    v_edot = vdot - Adj(v) * capAd(invX) * v_l - capAd(invX, vdot_l);

    gdot = g * se3alg(v);

    
    Out = [etadot; reshape(gdot,[],1); vdot];
end