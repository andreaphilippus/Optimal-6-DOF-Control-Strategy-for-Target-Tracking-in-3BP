function [t, g, v] = propagator_dim(X0, tspan, Var)
    opt = odeset('Reltol',Var.tol,'Abstol',Var.tol);
    
    for i = 1:Var.NumOrb-1
        tspan = [tspan; tspan(2:end) + tspan(end)];
    end

    %tspan = linspace(0,tspan(end),10);

    % Use ode45 for integrator
    [t, y] = ode45(@(t,X) Dyn_Comb_dim(t,X,Var), tspan, X0, opt);

    for i = 1:length(t)
        % The following quantities are in N frame
        R = y(i,1:3)';                     % Position, N frame 
        nu = y(i,7:9)';                    % Translational Velocity
        
        % Angular velcoity is in B frame
        w = y(i,4:6)';                     % Angular Velocity
    
        NB = reshape(y(i,10:18),3,3);     % DCM mapping from B to N frame

        g(:,:,i) = SE3(NB, R);
        v(:,i) = [w; NB' * nu];
    end
end