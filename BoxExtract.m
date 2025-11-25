%% To Box Topology

function M2 = BoxExtract(M1, Tz, Pull)

% BoxExtract aims to transform the arrow-type matrix to the box-type
% coupling matrix.
% Pull is the distance the box section moves toward the S port
%
%     3----4---L
%     |    |
% S---1----2

% By yellowbook, 2025-11-16

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

% annihilate M(N-pull-1,N-pull) to creat a 'box'
M2 = Rotate(M2, N-Pull-1, N-Pull-1,N-Pull, 'cross');


end