function M2 = AfterRotate(M1, Pull, theta1)

    N = length(M1) - 2;
    R = eye(N+2,N+2);
    R(N-Pull, N-Pull) = cos(theta1);
    R(N-Pull+1, N-Pull+1) = cos(theta1);
    R(N-Pull, N-Pull+1) = -sin(theta1);
    R(N-Pull+1, N-Pull) = sin(theta1);
    M2 = R*M1*R.';

    R = eye(N+2,N+2);
    theta_r = 0.5*atan(2*M2(N-Pull+1, N-Pull+2)/(M2(N-Pull+2, N-Pull+2)...
        - M2(N-Pull+1, N-Pull+1)));
    R(N-Pull+1, N-Pull+1) = cos(theta_r);
    R(N-Pull+2, N-Pull+2) = cos(theta_r);
    R(N-Pull+1, N-Pull+2) = -sin(theta_r);
    R(N-Pull+2, N-Pull+1) = sin(theta_r);
    M2 = R*M2*R.';
    % 
    % 
    M2 = Rotate(M2, N-Pull, N-Pull-1, N-Pull, 'cross');
    M2 = Rotate(M2, N-Pull-1, N-Pull+1, N-Pull, 'column');

%     M34(i) = M2(N-Pull-1,N-Pull);
%     M45(i) = M2(N-Pull,N-Pull+1);
%     M56(i) = M2(N-Pull+1,N-Pull+2);