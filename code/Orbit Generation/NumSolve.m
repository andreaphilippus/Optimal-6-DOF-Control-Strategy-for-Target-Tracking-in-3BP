function [t,X] = NumSolve(fun, X0, tspan, tol, event)
    
    if event == 0 % no event, propagates forward
        opt = odeset('Reltol',tol,'Abstol',tol);
    else % propagates for half segment for L2
        opt = odeset('Reltol',tol,'Abstol',tol,'Events',@(t,X)L2event(t,X));
    end
    
    % Use ode113 for integrator
    [t,X] = ode113(fun, tspan, X0, opt);
end