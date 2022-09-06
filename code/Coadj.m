function adstar = Coadj(a)

    adstar = [-CrossProd(a(1:3)) -CrossProd(a(4:6)); zeros(3) -CrossProd(a(1:3))];
end