function [x0, v0] = shiftnonlinCVX(st, cste, con, T, Var)
    
    %{
    [k1,l1] = CVXDynCont(st, cste, con, Var);                   % new 
    [k2,l2] = CVXDynCont(st + T/2*k1, cste + T/2*l1, con, Var); % new
    [k3,l3] = CVXDynCont(st + T/2*k2, cste + T/2*l2, con, Var); % new
    [k4,l4] = CVXDynCont(st + T*k3, cste + T*l3, con, Var);     % new
    st = st +T/6*(k1 + 2*k2 + 2*k3 + k4); % new  
    cste = cste + T/6*(l1 + 2*l2 + 2*l3 + l4);   % new vel

    x0 = full(st);
    v0 = full(cste);
    %}

    % or use numerical 
    X0 = [st; cste; 1];
    opt = odeset('Reltol',Var.tol,'Abstol',Var.tol);
    [~, hist] = ode45(@(t,x)CVXDynContaug(t,x,con,Var),[0 T],X0,opt);
    
    hist = hist(end,:)';
    x0 = hist(1:12); v0 = hist(13:18);
end