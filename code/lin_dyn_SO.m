function [A,B] = lin_dyn_SO(F,vars,var_ind)
ddot_ind = var_ind{1}; dot_ind = var_ind{2}; ind = var_ind{3}; u_ind = var_ind{4};

for j = 1:length(F)
    for k = 1:length(vars)
        PD(j,k) = diff(F(j),vars(k));
    end
end

M = PD(:,ddot_ind(1):ddot_ind(end));
D = PD(:,dot_ind(1):dot_ind(end));
S = PD(:,ind(1):ind(end));
W = PD(:,u_ind(1):u_ind(end));
A = simplify([zeros(size(M)), eye(size(M));
     M^-1*(-S), M^-1*(-D)]);
B = simplify([zeros(size(W));M^-1*(-W)]);
end