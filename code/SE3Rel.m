function rel = SE3Rel(to, from)
    % Relative SE3 pose: g_to/from
    [Cto, Rto] = invSE3(to);
    [Cfrom, Rfrom] = invSE3(from);

    rel = SE3(Cfrom'*Cto, Cfrom'*(Rto - Rfrom));

end