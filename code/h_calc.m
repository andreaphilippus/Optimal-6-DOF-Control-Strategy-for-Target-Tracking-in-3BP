function h = h_calc(xx, xd)
    
    R = xx(1:3,4);
    R_ref = xd(1:3,4);

    SB = xx(1:3,1:3);
    SB_ref = xd(1:3,1:3);
    
    R_e = R_ref - R;
    SB_e = SB_ref' * SB;
    h = SE3(SB_e, R_e);
end