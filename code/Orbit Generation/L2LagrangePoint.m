function L2x = L2LagrangePoint(mu)
    function y = colli_Lagrange(xeq)
        y = -((1-mu).*(xeq+mu))./(abs(xeq+mu).^3)...
            - mu.*(xeq-1+mu)./(abs(xeq-1+mu).^3) + xeq;
    end
    options = optimset('Display','iter');
    L2x = fzero(@colli_Lagrange,[1,2],options);
end