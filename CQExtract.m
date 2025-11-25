%% To CQ Topology

function M2 = CQExtract(M1, Tz, Pull, opt)

% CQExtract aims to transform the CT-type matrix to the CQ-type
% coupling matrix.
% Pull is the distance the CQ section moves toward the S port.
% opt = '13'(M24 = 0) or '24'(M13 = 0) or 'equal'(M13 = M24) to switch
% three types of CQ structures.
%       2----3
%       |\  /|
%       | \/ |
%       | /\ |
%       |/  \|
%  S----1----4----L
% By yellowbook, 2025-11-22

N = length(M1)-2;

if strcmp(opt, '13') == 1
   M1 = CTExtract(M1,Tz(1)*1i,Pull);
   M2 = CTExtract(M1,Tz(2)*1i,Pull);
end

if strcmp(opt, '24') == 1
   M1 = CTExtract(M1,Tz(1)*1i,Pull);
   M1 = CTExtract(M1,Tz(2)*1i,Pull-1);
   M2 = Rotate(M1, N-Pull-1, N-Pull+1, N-Pull, 'column');
end

if strcmp(opt, 'equal') == 1
    M1 = CTExtract(M1,Tz(1)*1i,Pull);
    M1 = CTExtract(M1,Tz(2)*1i,Pull-1);
    M2 = Rotate(M1, N-Pull-1, N-Pull+1, N-Pull, 'column');
    R = eye(N+2,N+2);
    theta_r = atan((M2(N-Pull, N-Pull+2)-M2(N-Pull-1, N-Pull+1))/...
        (M2(N-Pull-1, N-Pull) + M2(N-Pull+1, N-Pull+2)));
    R(N-Pull, N-Pull) = cos(theta_r);
    R(N-Pull+1, N-Pull+1) = cos(theta_r);
    R(N-Pull, N-Pull+1) = -sin(theta_r);
    R(N-Pull+1, N-Pull) = sin(theta_r);
    M2 = R*M2*R.';
end


end