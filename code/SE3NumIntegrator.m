function Out = SE3NumIntegrator(tvar, g0, v0, u, I, mu)
    
    % tvar  : time vector; shared with the reference time vector Nx1
    % g0    : Initial SE(3) state       4x4
    % v0    : Initial velocity vector   6x1
    % u     : Control vector            6xN
        % u will be changed after verifying uncontrolled dynamics works
    
    opt = odeset('RelTol',1e-12, 'AbsTol', 1e-12);
    temp = ode45(@(t,X)aug(t, X, u, I, mu), [0 tvar(end)], [reshape(g0,[],1); v0], opt);

    temp = deval(temp, tvar);
    
    Out.g = []; Out.v = [];
    for i = 1:length(tvar)
        Out.g(:,:,i) = reshape(temp(1:16, i), 4, 4);
        Out.v(:,i) = temp(17:end, i);
    end

    function Xdot = aug(t, X, u, I, mu)
        % X             : Augmented, 32x1
            % X(1:16)   : vec(g)
            % X(17:end) : v

        g = reshape(X(1:16), 4, 4);
        v = X(17:end);

        gdot = SE3Kinematics(g, v);
        vdot = SE3Dynamics(v, g, u, I, mu);

        Xdot = [reshape(gdot, [], 1); vdot]
    end
end