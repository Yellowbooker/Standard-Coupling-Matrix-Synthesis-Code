%% From Coupling Matrix to E, P, F

function FilterPd = Predistortion(Filter, Q, FBW, effectiveQ)

% Predistortion aims to calculate the predistorted characteristic 
% polynomials using the general Chebyshev characteristic polynomials. 

% Filter is the general Chebyshev characteristic polynomials
% Q is the actual Q factor of the filter
% FBW is the fractional bandwidth of the filter
% effectiveQ is the effective Q factor of the filter after the predistortion
% FilterPd is the predistorted characteristic polynomials

% By yellowbook, 2024-08-15

N = length(Filter.Es) - 1;
v = ones(N,1);
% v = [1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 0.4 0.4].';
Newpole = sort(roots(Filter.Es),'ComparisonMethod','real') + (1/FBW/Q - 1/FBW/effectiveQ)*v;
FilterPd.Es = poly(Newpole);
FilterPd.Ps = Filter.Ps;

Ns = 1000;
for m = 1:Ns
    k(m) = abs(polyval(FilterPd.Es,-1.05i + 2.1i/(Ns-1)*(m-1))/polyval(FilterPd.Ps,-1.05i + 2.1i/(Ns-1)*(m-1)));
end

FilterPd.epsilon = 1/min(k);
FilterPd.epsilonR = 1;

Fs2 = conv(FilterPd.Es,conj(FilterPd.Es).*((-1).^[N:-1:0])) - ...
    conv(FilterPd.Ps,conj(FilterPd.Ps).*((-1).^[N:-1:0]))/(FilterPd.epsilon*FilterPd.epsilon);
RootF = roots(Fs2);
RootF = RootF(real(RootF) < 0);
FilterPd.Fs = poly(RootF);
FilterPd.F22s = ((-1)^N)*conj(FilterPd.Fs).*(-1).^[N:-1:0];

end