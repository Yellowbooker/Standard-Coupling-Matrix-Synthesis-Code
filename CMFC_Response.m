%% Coupling Matrix to S-parameters

function [S11, S12, S22] = CMFC_Response(M, C, G, w)

% CMFC_Response aims to calculate the S-parameters using the coupling
% matrix

% M is the coupling matrix
% C is the capacitor matrix
% G is the port admittance matrix

% By yellowbook, 2024-07-21

for k = 1:length(w)
    A = G + C.*1i*w(k) + 1i*M;
    AP = inv(A);
    S12(k) = 2*AP(end,1);
    S11(k) = -1 + 2*AP(1,1);
    S22(k) = -1 + 2*AP(end,end);
end

end
