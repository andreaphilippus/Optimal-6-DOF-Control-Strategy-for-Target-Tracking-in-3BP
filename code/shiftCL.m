function [t0, x0, v0, u0] = shiftCL(T, t0, x0, v0, u, f)
st = full(x0);
cste = full(v0);
con = u(1,:)';
%f_value = f(st, cste, con);
%st = st+ (T*f_value);

%
[k1,l1] = f(st, cste, con);                   % new 
[k2,l2] = f(st + T/2*k1, cste + T/2*l1, con); % new
[k3,l3] = f(st + T/2*k2, cste + T/2*l2, con); % new
[k4,l4] = f(st + T*k3, cste + T*l3, con);     % new
st = st +T/6*(k1 + 2*k2 + 2*k3 + k4); % new  
cste = cste + T/6*(l1 + 2*l2 + 2*l3 + l4);   % new vel

x0 = full(st);
v0 = full(cste);
t0 = t0 + T;
u0 = [u(2:size(u,1),:);u(size(u,1),:)];
%}

%{
[stdot, vldot] = f(st, cste, con);
st = st+ (T*stdot);
vl = cste + (T*vldot);
t0 = t0 + T;
x0 = full(st);
v0 = full(vl);
u0 = [u(2:size(u,1),:);u(size(u,1),:)];
%}
end