function Collinear_Lagrange_Point
    m_Earth = 5.974E24;
    m_Luna = 7.348E22;
    mu = m_Luna / (m_Earth + m_Luna);
    function y = colli_Lagrange(xeq)
        y = -((1-mu).*(xeq+mu))./(abs(xeq+mu).^3)...
            - mu.*(xeq-1+mu)./(abs(xeq-1+mu).^3) + xeq;
    end
    options = optimset('Display','iter');
    L1 = fzero(@colli_Lagrange,[0.01,0.98],options);
    L2 = fzero(@colli_Lagrange,[1,2]);
    L3 = fzero(@colli_Lagrange,[-0.1,-2]);
    disp(L1)
    disp(L2)
    disp(L3)
end