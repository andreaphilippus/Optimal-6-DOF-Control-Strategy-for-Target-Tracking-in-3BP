function G = Gfunc(eta)

    thetavec = eta(1:3); beta = eta(4:6);

    thetaX = CrossProd(thetavec);
    theta = norm(thetavec);

    S = Sfunc(thetanorm, thetaX);

    A = eye(3) + 0.5*thetaX + (1/theta^2 - 1 + cos(theta)*2*theta*sin(theta))*(thetaX'*thetaX);

    T = 0.5*CrossProd(S*beta)*A +...
        (1/theta^2 - (1+cos(theta))/(2*theta*sin(theta)))*(thetavec*beta' + (thetavec'*beta)*A) -...
        (1 + cos(theta))*(theta - sin(theta))/(2*theta*(sin(theta))^2)*S*beta*thetavec' +...
        ((1 + cos(theta))*(theta + sin(theta))/(2*theta^3 *(sin(theta))^2) - 2/theta^4)*thetavec'*beta*thetavec*thetavec';

    G = [A zeros(3); T A];
end