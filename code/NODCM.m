function NO = NODCM(t, Var)

    g = Var.n * t;
    
    NO = [cos(g) -sin(g) 0; sin(g) cos(g) 0; 0 0 1];
end