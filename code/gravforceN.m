function GravF = gravforceN(X, Var)
    x = X(1); y = X(2); z = X(3); % Position
    xdot = X(7); ydot = X(8); % Velocity

    mu = Var.mu;
    d13 = norm([x + mu, y, z]);
    d23 = norm([x - 1 + mu, y, z]);

    %%
    
    ax = x + 2*ydot - (1-mu)*(x + mu)/(d13^3) - mu*(x - 1 + mu)/(d23^3);
    ay = y - 2*xdot - (1-mu)*y/(d13^3) - mu*y/(d23^3);
    az = - (1-mu)*z/(d13^3) - mu*z/(d23^3);
    
    GravF = [ax; ay; az];
end