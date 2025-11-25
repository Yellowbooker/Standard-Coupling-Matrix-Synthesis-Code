function M1 = NormalizeCM(M1)
    % normalized operation
    N = length(M1)-2;
    for i = 1:N+1
        if real(M1(i,i+1)) < 0
            R = eye(N+2,N+2);
            R(i+1,i+1) = -1;
            M1 = R*M1*R.';
        end
    end
    for i = 1:N+2
        for j = 1:N+2
            if abs(M1(i,j)) < 1e-5
                M1(i,j) = 0;
            end
        end
    end
end
