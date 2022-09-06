function F_G_N = GravF_Inert(X, R_E, R_L, Var)
    
    x = X(1); y = X(2); z = X(3); % Position
    R_s = [x; y; z];

    R_E = [R_E(1); R_E(2); 0];
    R_L = [R_L(1); R_L(2); 0];
    
    r_13 = R_s - R_E; % Position vector originated from Earth
    r_23 = R_s - R_L; % Position vector originated from Luna

    %% Force by Earth

    F_G_E = -Var.m1 / norm(r_13)^3 * r_13;

    %% Force by Luna

    F_G_L = -Var.m2 / norm(r_23)^3 * r_23;

    %% Total Force
    
    F_G_N = (F_G_E + F_G_L)';
    
end