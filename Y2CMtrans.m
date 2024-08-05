%% Y-parameters to transversal coupling matrix

function M_trans = Y2CMtrans(Y)

% Y2CMtrans aims to calculate the transversal coupling matrix (N + 2) 
% using the Y-parameters

% By yellowbook, 2024-07-21

[residueY11,Ypole,K11] = residue2(Y.Y11n,Y.Yd);
[residueY12,Ypole,K12] = residue2(Y.Y12n,Y.Yd);
[residueY22,Ypole,K22] = residue2(Y.Y22n,Y.Yd);

% ======================= Y to TCM =======================
N = length(Y.Yd) - 1;
M_trans = zeros(N+2,N+2);
for k=1:N
    M_trans(k+1,k+1) = 1i*Ypole(k,1);
    if abs(residueY11(k,1)) >= abs(residueY22(k,1))
        M_trans(1,k+1) = sqrt(residueY11(k,1));
        M_trans(N+2,k+1) = residueY12(k,1)/sqrt(residueY11(k,1));
    else if abs(residueY11(k,1)) < abs(residueY22(k,1))
        M_trans(N+2,k+1) = sqrt(residueY22(k,1));
        M_trans(1,k+1) = residueY12(k,1)/sqrt(residueY22(k,1));
        end
    end
    M_trans(k+1,1) = M_trans(1,k+1);
    M_trans(k+1,N+2) = M_trans(N+2,k+1);
end

if isempty (K11) && isempty (K22)
else
    M_trans(1,1) = K11/1i;
    M_trans(end,end) = K22/1i;
end

if length(roots(Y.Y12n)) == N
    M_trans(1,end) = K12/1i;
    M_trans(end,1) = K12/1i;
end

end