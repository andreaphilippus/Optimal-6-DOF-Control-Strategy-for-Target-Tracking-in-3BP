function S = Sfunc(theta, thetaX)
    
    % investigate if thetaX ^2 in Eq 85 is product with transpose or not
    S = eye(3) + (1 - cos(theta))/theta^2 * thetaX + (theta - sin(theta))/theta^3 * (thetaX'*thetaX);
end