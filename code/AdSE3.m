function AdSE3 = AdSE3(g)
    
    [P, b] = invSE3(g);

    AdSE3 = [P zeros(3,3); CrossProd(b)*P, P];
end