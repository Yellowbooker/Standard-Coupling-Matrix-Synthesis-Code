%% To Parallel Coupled Pair Topology

function M = LossyCM2PCP(M)

% lossyCM2PCP aims to transform M to the parallel coupled pair 
% coupling matrix

% By yellowbook, 2024-08-02

N = length(M) - 2;

if mod(N,2) == 0
    for i = 1:N/2
        M = Rotate(M, 1, 2*i , 2*i+1, 'column');
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