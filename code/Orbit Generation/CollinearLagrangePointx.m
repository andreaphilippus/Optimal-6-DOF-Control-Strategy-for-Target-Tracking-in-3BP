function [L1x, L2x, L3x] = CollinearLagrangePointx(mu)
   
    colli_Lagrange = @(xeq) -((1-mu).*(xeq+mu))./(abs(xeq+mu).^3)...
            - mu.*(xeq-1+mu)./(abs(xeq-1+mu).^3) + xeq; 
    L1x = fzero(colli_Lagrange,0.5);
    L2x = fzero(colli_Lagrange,1.5);
    L3x = fzero(colli_Lagrange,-2);
end