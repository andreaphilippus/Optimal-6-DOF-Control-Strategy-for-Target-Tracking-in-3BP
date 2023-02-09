function v_BP = SE3RelVel(v_B, gPB, v_P)
    
    v_BP = v_B - AdSE3(gPB) * v_P; 

end