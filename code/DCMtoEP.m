function EP = DCMtoEP(C)
    e4 = sqrt(1/4 * (1 + trace(C)));

    e1 = (C(2,3) - C(3,2))/4/e4;
    e2 = (C(3,1) - C(1,3))/4/e4;
    e3 = (C(1,2) - C(2,1))/4/e4;

    EP = [e1;e2;e3;e4];
end