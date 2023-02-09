function vdot = SE3Dynamics(v, g, u, I, mu)
    
    [C, R] = invSE3(g);

    [~, UXX] = U_Partial(R(1), R(2), R(3), mu);

    % CR3BP Gravitational Force
    %d = sqrt((R(1) + mu).^2 + R(2).^2 + R(3).^2);
    %r = sqrt((R(1)-1+mu).^2 + R(2).^2 + R(3).^2);
    %
    d = [R(1) + mu; R(2); R(3)];
    r = [(R(1)-1+mu); R(2); R(3)];

    F_g_Earth = -(1-mu)/norm(d)^3 * d;

    % In ECI
    F_g_Moon = -mu * (r/norm(r)^3 + [1-mu, 0, 0]'/norm(1-mu)^3);
    F_g_Moon = F_g_Moon + CrossProd([0 0 ]')* C' *v(4:6)

    Fg = F_g_Moon + F_g_Earth;
    %}
    %Fg = UXX*R + [0 2 0; -2 0 0; 0 0 0] * C' * v(4:6);

    %Fg = -(1-mu)*R/norm(d)^3 - mu*r/norm(r)^3;

    %Fg = -((1-mu)/norm(d)^3 + mu/norm(r)^3) * R + 2*mu*(d/norm(d)^3);

    phi = [0;0;0; Fg];
    vdot = I\(Coadj(v)*I*v + u ) + phi ;
   
end