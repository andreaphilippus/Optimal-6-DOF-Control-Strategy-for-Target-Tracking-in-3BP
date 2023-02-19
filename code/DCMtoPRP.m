function PRP = DCMtoPRP(C)
    theta = acos(1/2*(trace(C)-1));
    lambda = 1/2/sin(theta) * [C(2,3)-C(3,2); C(3,1)-C(1,3); C(1,2)-C(2,1)];
    PRP = [theta; lambda];
end