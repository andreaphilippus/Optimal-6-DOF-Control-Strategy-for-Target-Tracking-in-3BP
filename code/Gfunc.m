function G = Gfunc(eta)

    thetavec = eta(1:3); beta = eta(4:6);

    thetaX = CrossProd(thetavec);
    theta = norm(thetavec);
<<<<<<< Updated upstream
=======
    
    
    if theta == 0
        thetaX = zeros(3);
    else
        thetaX = CrossProd(thetavec);
    end
    
>>>>>>> Stashed changes

    S = Sfunc(thetanorm, thetaX);

    A = eye(3) + 1/2*thetaX + (1/theta^2 - (1 + cos(theta))/(2*theta*sin(theta)))*(thetaX*thetaX);

<<<<<<< Updated upstream
    T = 0.5*CrossProd(S*beta)*A +...
        (1/theta^2 - (1+cos(theta))/(2*theta*sin(theta)))*(thetavec*beta' + (thetavec'*beta)*A) -...
        (1 + cos(theta))*(theta - sin(theta))/(2*theta*(sin(theta))^2)*S*beta*thetavec' +...
        ((1 + cos(theta))*(theta + sin(theta))/(2*theta^3 *(sin(theta))^2) - 2/theta^4)*thetavec'*beta*thetavec*thetavec';
=======
    T = 1/2*CrossProd(S*beta)*A +...   
        (1/theta^2 - (1+cos(theta))/(2*theta*sin(theta)))*(thetavec*beta' + (thetavec'*beta)*A) -...
        (1 + cos(theta))*(theta - sin(theta))/(2*theta*(sin(theta))^2)*S*beta*thetavec' +...
        ((1 + cos(theta))*(theta + sin(theta))/(2*theta^3 *(sin(theta))^2) - 2/theta^4)*(thetavec'*beta)*(thetavec*thetavec');
    
    if eta == zeros(6,1)
        G = zeros(6);
    else
        G = [A zeros(3); T A];
    end
end
>>>>>>> Stashed changes

    G = [A zeros(3); T A];
end