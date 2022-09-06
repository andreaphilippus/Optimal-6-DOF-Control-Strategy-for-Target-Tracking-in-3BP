function Xdot = CR3BP_EoM(t, X, mu)

    x = X(1); y = X(2); z = X(3); % Position
    xdot = X(4); ydot = X(5); zdot = X(6); % Velocity

    d13 = sqrt((x + mu)^2 + y^2 + z^2); % Distance between S/C and Earth
    d23 = sqrt((x-1+mu)^2 + y^2 + z^2); % Distance between S/C and Luna

    Xdot = [xdot; ydot; zdot; % [xdot ydot zdot xddot yddot zddot]
        2*ydot + x - (1-mu)*(x+mu)/d13^3 - mu*(x-1+mu)/d23^3;
        y - 2*xdot - (1-mu)*y/d13^3 - mu*y/d23^3;
        - (1-mu)*z/d13^3 - mu*z/d23^3;];
end