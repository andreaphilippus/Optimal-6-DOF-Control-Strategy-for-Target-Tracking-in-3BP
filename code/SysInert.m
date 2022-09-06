function [R_E, R_L, R_L2] = SysInert(tspan, Var)

    temp = Var.n*tspan;
    
    % Earth Trajectory
    rE0 = -Var.D1;
    R_E = [rE0*cos(temp) rE0*sin(temp)];

    % Luna Trajectory 
    rL0 = Var.lstar - Var.D1;
    R_L = [rL0*cos(temp) rL0*sin(temp)];

    % L2 Trajectory
    L2x = Var.L2x * Var.lstar;
    R_L2 = [L2x*cos(temp) L2x*sin(temp)];
    
end