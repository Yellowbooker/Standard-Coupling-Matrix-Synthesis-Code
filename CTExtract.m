%% To CT Topology

function M2 = CTExtract(M1, Tz, Pull)

% TriExtract aims to transform the arrow-type matrix to the CT-type
% coupling matrix.
% Pull is the distance the CT section moves toward the S port.
% When the pulls of the two CT sections are the same, the two CT sections
% would becomes a CQ section.

% By yellowbook, 2024-08-04

N = length(M1) - 2;
R = eye(N+2,N+2);
theta_r = atan(M1(N,N+1)/(Tz/1i + M1(N+1,N+1)));
R(N, N) = cos(theta_r);
R(N+1, N+1) = cos(theta_r);
R(N, N+1) = -sin(theta_r);
R(N+1, N) = sin(theta_r);
M2 = R*M1*R.';

for i = 1:Pull
    M2 = Rotate(M2, N-i+2, N-i, N-i+1, 'row');
end

end