function gdot = SE3RelKin(grel, vto, vfrom)
    % grel = g_B/P
    % vto = v_B
    % vfrom = v_P
    % B relative to P

    gdot = grel*vee(vto) - vee(vfrom)*grel;
end