function gdot = SE3Kinematics(g, v)

    gdot = g*se3alg(v);
end