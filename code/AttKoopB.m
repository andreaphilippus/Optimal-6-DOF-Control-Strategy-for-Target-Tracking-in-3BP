function B = AttKoopB(n, BS, omega, J)
    syms u [3, 1] 

    B = zeros(9,3);
    for j = 1:n-1

        temp = 0;
        for i = 1:j
            temp = temp + (CrossProd(omega))^(i - 1) * CrossProd(J * u) * (CrossProd(omega))^(n - i);
        end

        temp = reshape(BS * temp, [], 1);

        B = double([B; diff(temp, u(1)) diff(temp, u(2)) diff(temp, u(3))]);
    end

    B = [zeros(9, 9*(n-1)); eye(9*(n-1))] * B;
end