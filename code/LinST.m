function X = LinST(g)
    
    [SB, R] = invSE3(g);
    X = [reshape(SB, [], 1); R];
end