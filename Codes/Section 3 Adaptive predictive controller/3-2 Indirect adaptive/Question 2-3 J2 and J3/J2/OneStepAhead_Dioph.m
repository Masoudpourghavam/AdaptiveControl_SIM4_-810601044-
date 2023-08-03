function [F,G] = dioph(A, d, n)
    A(1, n+2:n+2+d) = zeros(1, d+1);

    F = zeros(1, d);
    F(1) = 1;
    for i = 2:d
        for j = 1:i-1
            F(i) = F(i) - A(i-j)*F(j);
        end
    end

    G = zeros(1, n);
    for i = 1:n
        for j = 1:d
            G(i) = G(i) + A(n+1+i-j)*F(j);
        end
    end
end
