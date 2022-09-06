function EA323 = DCMtoEA323(C)
    disp(C)
    t1 = atan2(C(3,2), C(3,1));
    t2 = acos(C(3,3));
    t3 = atan2(C(2,3), -C(1,3));
    
    EA323 = [t1 t2 t3]';
end