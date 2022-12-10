function x = SE_addition(x1, x2)

    att1 = x1(1:3,1:3);
    att2 = x2(1:3,1:3);
    att = att1 * att2;
    
    R1 = x1(1:3,4);
    R2 = x2(1:3,4);
    R = R1 + R2;

    x = SE3(att, R);
end