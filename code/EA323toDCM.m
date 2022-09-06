function DCM = EA323toDCM(t1, t2, t3)
    
    R3_1 = [cos(t3) sin(t3) 0; -sin(t3) cos(t3) 0; 0 0 1];
    R2_2 = [cos(t2) 0 -sin(t2); 0 1 0; sin(t2) 0 cos(t2)];
    R3_3 = [cos(t1) sin(t1) 0; -sin(t1) cos(t1) 0; 0 0 1];

    DCM = R3_1*R2_2*R3_3;

end