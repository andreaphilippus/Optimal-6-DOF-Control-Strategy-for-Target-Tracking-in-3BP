function Ad = Adj(a)

    Ad = [CrossProd(a(1:3)) zeros(3); CrossProd(a(4:6)) CrossProd(a(1:3))];
end