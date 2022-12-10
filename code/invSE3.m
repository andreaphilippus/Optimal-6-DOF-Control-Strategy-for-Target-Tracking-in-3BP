function [C, R] = invSE3(g)
    C = g(1:3,1:3);
    R = g(1:3,4);
end