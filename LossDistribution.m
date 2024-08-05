%% Loss Redistribution

function M = LossDistribution(M)

% LossDistribution aims to obtain the coupling matrix with the parallel coupled
% resonator pairs having the same Q factor.

% By yellowbook, 2024-08-05

N = length(M) - 2;


if mod(N,2) == 0
    for i = 1:N/2
        k = M(2*i, 2*i + 1)/(M(2*i, 2*i) - M(2*i + 1, 2*i + 1));
        theta_r = atan(-2*k - sqrt(4*k*k + 1));
%         k = [1, M(2*i, 2*i + 1)/(M(2*i, 2*i) - M(2*i + 1, 2*i + 1)), -1];
%         zero = roots(k);
%         theta_r = atan(zero(1));
        R = eye(N + 2);
        R(2*i, 2*i) = cos(theta_r);
        R(2*i+1, 2*i+1) = cos(theta_r);
        R(2*i, 2*i+1) = -sin(theta_r);
        R(2*i+1, 2*i) = sin(theta_r);
        M = R*M*R.';
    end
end
% if mod(N,2) == 0
%     for i = 1:N/2
%         M = Rotate(M, 1, i+1 , N - i + 2, 'column');
%     end
% end

% for i = 1:N+2
%     for j = 1:N+2
%         if abs(M(i,j)) < 1e-4
%             M(i,j) = 0;
%         end
%     end
% end
end