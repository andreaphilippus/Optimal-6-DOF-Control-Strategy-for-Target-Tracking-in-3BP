function [U_X, U_XX] = U_Partial(pos1, pos2, pos3, mu)
    % By default, input argument is [pos, sys]
    % [pos, pos2] can be substituted by X and Y for planar;
   
    switch nargin
        case 2
            x = pos1(1); y = pos1(2); z = pos1(3);
            sys = pos2;

        case 3
            x = pos1; y = pos2; z = zeros(size(x));
            sys = pos3;

        case 4
            x = pos1; y = pos2; z = pos3;
           
    end

    d = sqrt((x + mu).^2 + y.^2 + z.^2);
    r = sqrt((x-1+mu).^2 + y.^2 + z.^2);

    Ux = -(1-mu)*(x+mu)/d^3 - mu*(x-1+mu)/r^3 + x;
    Uy = -(1-mu)*y/d^3 - mu*y/r^3 + y;
    Uz = -(1-mu)*z/d^3 - mu*z/r^3;

    Uxx = -(1-mu)/d^3 - mu/r^3 + 3*(1-mu)*(x+mu)^2/d^5 + 3*mu*(x-1+mu)^2/r^5 + 1;
    Uyy = -(1-mu)/d^3 - mu/r^3 + 3*(1-mu)*y^2/d^5 + 3*mu*y^2/r^5 + 1;
    Uzz = -(1-mu)/d^3 - mu/r^3 + 3*(1-mu)*z^2/d^5 + 3*mu*z^2/r^5 + 1;

    Uxy = 3*(1-mu)*(x+mu)*y/d^5 + 3*mu*(x-1+mu)*y/r^5;
    Uxz = 3*(1-mu)*(x+mu)*z/d^5 + 3*mu*(x-1+mu)*z/r^5;
    Uyz = 3*(1-mu)*z*y/d^5 + 3*mu*z*y/r^5;

    U_X = [Ux; Uy; Uz];

    U_XX = [Uxx Uxy Uxz; Uxy Uyy Uyz; Uxz Uyz Uzz];
 
end