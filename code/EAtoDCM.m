function [rot,R] = EAtoDCM(theta, seq)
    R1 = @(t) [1 0 0
        0 cos(t) sin(t)
        0 -sin(t) cos(t)];
    R2 = @(t) [cos(t) 0 -sin(t)
        0 1 0
        sin(t) 0 cos(t)];
    R3 = @(t) [cos(t) sin(t) 0
        -sin(t) cos(t) 0
        0 0 1];
    R = {R1 R2 R3};
    rot = eye(3);
    for i = 3:-1:1
        rot = rot*R{seq(i)}(theta(i));
    end
end