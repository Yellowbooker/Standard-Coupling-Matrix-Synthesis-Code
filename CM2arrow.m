%% To Arrow Topology

function M = CM2arrow(M)

% CM2arrow aims to transform M to the arrow-type coupling matrix

% By yellowbook, 2024-08-02

N = length(M) - 2;

for i = 1:(N - 1)
    for j = 1:(N - i)
        M = Rotate(M, i, N - j + 2 , N + 1 - j, 'column');
    end
end

% % normalized operation
% for i = 1:N
%     if real(M(i,i+1)) < 0
%         R = eye(N+2,N+2);
%         R(i+1,i+1) = -1;
%         M = R*M*R.';
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