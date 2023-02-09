function out = AdSE3inv(g)
    
    [C, R] = invSE3(g);
    out = [C' zeros(3,3); -C'*CrossProd(R) C'];
end