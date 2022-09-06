function Xdot = StateSTMdot(t, X, mu)
    
    %% Position/Velocity State CR3BP Equations of Motion
    x = X(1); y = X(2); z = X(3); % Position
    xdot = X(4); ydot = X(5); zdot = X(6); % Velocity
    
    d13 = sqrt((x + mu)^2 + y^2 + z^2); % Distance between S/C and Earth
    d23 = sqrt((x-1+mu)^2 + y^2 + z^2); % Distance between S/C and Luna

    Statedot = [xdot; ydot; zdot; % [xdot ydot zdot xddot yddot zddot]
        2*ydot + x - (1-mu)*(x+mu)/d13^3 - mu*(x-1+mu)/d23^3;
        y - 2*xdot - (1-mu)*y/d13^3 - mu*y/d23^3;
        - (1-mu)*z/d13^3 - mu*z/d23^3;];
    
    %% State Transition Matrix
    STM = reshape(X(7:end), 6,6); % STM

    % Partial derivatives of Pseudo-potential
    Uxx = 1 - (1-mu)/d13^3 - mu/d23^3 + 3*(1-mu)*(x+mu)^2 /d13^5 + 3*mu*(x-1+mu)^2 / d23^5;
    Uyy = 1 - (1-mu)/d13^3 - mu/d23^3 + 3*(1-mu)*y^2 /d13^5 + 3*mu*y^2 / d23^5;
    Uzz = -(1-mu)/d13^3 - mu/d23^3 + 3*(1-mu)*z^2 /d13^5 + 3*mu*z^2 / d23^5;
    Uxy = 3*(1-mu)*(x+mu)*y /d13^5 + 3*mu*(x-1+mu)*y / d23^5;
    Uyx = Uxy;
    Uxz = 3*(1-mu)*(x+mu)*z /d13^5 + 3*mu*(x-1+mu)*z / d23^5;
    Uzx = Uxz;
    Uyz = 3*(1-mu)*y*z /d13^5 + 3*mu*y*z / d23^5;
    Uzy = Uyz;

    UXX = [Uxx Uxy Uxz; Uyx Uyy Uyz; Uzx Uzy Uzz];
    Omega = [0 2 0; -2 0 0; 0 0 0];
    
    % Time-varying Matrix A(t)
    A = [zeros(3) eye(3); UXX Omega];

    STMdot = reshape(A*STM, [], 1);
    
    % X = [x y z xdot ydot zdot reshape(STM)]
    % Xdot = [xdot ydot zdot xddot yddot zddot reshape(STMdot)]
    Xdot = [Statedot; STMdot];
end