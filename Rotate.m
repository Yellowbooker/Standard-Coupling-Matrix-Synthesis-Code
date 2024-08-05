%% Performs A Givens Transformation on The Coupling Matrix M1

function M2 = Rotate(M1, k, i, j, opt)

% Rotate is to perform a givens transformation on the coupling matrix M1

% M1 is the input coupling matrix
% M2 is the output coupling matrix
% k is the index of the target element that needs to be transform to 0
% [i, j] is the rotate pivot
% opt is the type of the target element
% if opt = 'row' the target element is [i, k]
% if opt = 'column' the target element is [k, i]
% if opt = 'cross' the target element is [i, j]

% By yellowbook, 2024-08-02

N = length(M1) - 2;
R = eye(N + 2);
flag = 0;

if strcmp(opt, 'row') == 1
   theta_r = atan(M1(i, k)/M1(j, k));
   flag = 1;
end

if strcmp(opt, 'column') == 1
   theta_r = atan(M1(k, i)/M1(k, j));
   flag = 1;
end

if strcmp(opt, 'cross') == 1
   theta_r = 0.5*atan(2*M1(i, j)/(M1(j, j) - M1(i, i)));
   flag = 1;
end

if flag == 1
    R(i, i) = cos(theta_r);
    R(j, j) = cos(theta_r);
    R(i, j) = -sin(theta_r);
    R(j, i) = sin(theta_r);
    M2 = R*M1*R.';
else
    M2 = M1;
end


