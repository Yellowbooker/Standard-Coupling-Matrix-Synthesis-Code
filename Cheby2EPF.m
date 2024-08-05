%% General Chebyshev Filtering Function to E, P, F

function Filter = Cheby2EPF(Fw, Pw, RL)

% Cheby2EPF aims to calculate the characteristic polynomials (Es, Ps, Fs, Ep,
% Epr) and of the general chebyshev fiter with the inband return loss of RL 

% Fw and Pw is the normalized polynomials of the general chebyshev fuction
% CN in w-domian (CN = Pw/Fw)
% Filter.Es: the normalized characteristic polynomial constructed by the 
% poles
% Filter.Fs: the normalized characteristic polynomial constructed by the
% zeros of S11
% Filter.F22s: the normalized characteristic polynomial constructed by the
% zeros of S22
% Filter.Ps: the normalized characteristic polynomial constructed by the 
% zeros of S12 and S21
% Filter.epsilon: the normalized constant of S12 and S21
% Filter.epsilonR: the normalized constant of S11 and S22
% S11 = Fs/Es/epsilonR; S21 = S12 = Ps/Es/epsilon; S22 = F22s/Es/epsilon; 

% By yellowbook, 2024-07-21

N = length(Fw) - 1; % order
NTz = length(roots(Pw));
Filter.epsilon = 1/sqrt(10^(RL/10)-1)*abs(polyval(Pw,1)/polyval(Fw,1));
if NTz == N
    Filter.epsilonR = Filter.epsilon/sqrt(Filter.epsilon*Filter.epsilon-1);
else
    Filter.epsilonR = 1;
end
Ew2 = conv(Pw,Pw)/(Filter.epsilon*Filter.epsilon) + ...
    conv(Fw,Fw)/(Filter.epsilonR*Filter.epsilonR);
RootE = roots(Ew2);
RootE = RootE(imag(RootE) > 0);
Ew = poly(RootE);

Filter.Es = Ew./1i.^[N:-1:0]; % from w-domain to s-domain
Filter.Es = Filter.Es./Filter.Es(1);
Filter.Fs = Fw./1i.^[N:-1:0]; % from w-domain to s-domain
Filter.Fs = Filter.Fs./Filter.Fs(1);
Filter.F22s = ((-1)^N)*conj(Filter.Fs).*(-1).^[N:-1:0]; % from w-domain to s-domain
% Filter.F22s = Filter.F22s./Filter.F22s(1);
Filter.Ps = Pw./1i.^[zeros(1, N - NTz),[NTz:-1:0]]; % from w-domain to s-domain
Filter.Ps = Filter.Ps./Filter.Ps(end - NTz);
if mod(N - NTz, 2) == 0
    Filter.Ps = Filter.Ps*1i;
end

end





