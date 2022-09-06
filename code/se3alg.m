function vee = se3alg(a)

    vee = [CrossProd(a(1:3)) a(4:6); zeros(1,4)];
end