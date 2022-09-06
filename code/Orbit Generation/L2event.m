function [val, isterminal, direction] = L2event(t,X)
    isterminal = 1; % Ends after the s/c comes back to y=0 again
    val = X(2); % y = 0 needed
    direction = 1; 
end