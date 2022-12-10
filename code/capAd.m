function out = capAd(H, eps)
    if nargin == 2
        out = vee(H*se3alg(eps)*inv(H));
    else
        [P, b] = invSE3(H);

        out = [P zeros(3); CrossProd(b)*P P];
    end
    
end