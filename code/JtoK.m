function K = JtoK(J)
    I1 = J(1,1); I2 = J(2,2); I3 = J(3,3);

    K = [(I3-I2)/I1; (I1-I3)/I2; (I2-I1)/I3];
end