%% Synthesis General Chebyshev Filtering Function

function [Fw, Pw] = General_Chebyshev(Pwz, N)

% General_Chebyshev aims to calculate the numerator polynomial(Pw) and the 
% denominator polynomial(Fw) of the general chebyshev fitering function(CN):
% CN = Pw/Fw

% Tz is a vector contains the transmission zeros in w-domain
% N is the order of the filter
% Fw and Pw is the normalized polynomials

% By yellowbook, 2024-07-21

Pw = [zeros(1, N - length(Pwz)) poly(Pwz)];

Pwz = [Pwz,inf(1, N-length(Pwz))];
U1 = [1, -1/Pwz(1)];
V1 = [sqrt(1 - 1/(Pwz(1)*Pwz(1)))];
U2 = U1;
V2 = V1;
for k = 2:N
    U2 = conv(U1, [1,-1/Pwz(k)]) + conv([1,0,-1],V1.*sqrt(1 - 1/(Pwz(k)*Pwz(k))));
    V2 = conv(V1, [1,-1/Pwz(k)]) + U1.*sqrt(1 - 1/(Pwz(k)*Pwz(k)));
    U1 = U2;
    V1 = V2;
end
Fw = U1./(U1(1));
% V = V1./(V1(1));
% V1 = V.*1i.^[0:1:length(V)-1];

end





