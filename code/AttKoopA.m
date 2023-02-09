function A = AttKoopA(n)
    
    A = zeros(9*n, 9*2);
    A(1:9,10:18) = eye(9);
    temp = zeros(9*n, 9);
    for i = 2:n-1
        temp(9*(i-1)+1:9*(i-1)+9, 1:9) = eye(9);

        A = [A temp];
        temp = zeros(9*(i-1), 9);
    end
end