%% From Coupling Matrix to E, P, F

function [FilterPd,ab_final] = AdaptivePredistortion(Filter, Q, FBW, effectiveQ)

% Predistortion aims to calculate the predistorted characteristic 
% polynomials using the general Chebyshev characteristic polynomials. 

% Filter is the general Chebyshev characteristic polynomials
% Q is the actual Q factor of the filter
% FBW is the fractional bandwidth of the filter
% effectiveQ is the effective Q factor of the filter after the predistortion
% FilterPd is the predistorted characteristic polynomials

% By yellowbook, 2024-08-16

N = length(Filter.Es) - 1;
v = ones(N,1);
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
w = [-1:0.01:1];
Num = length(w);
Before = abs(polyval(FilterPd.Ps,1i*w)./polyval(FilterPd.Es,1i*w)/FilterPd.epsilon);
Before = Before/max(Before);
a = 1:-0.1:0.1;
b = 1:0.1:1.9;

for i = 1:length(a)
    for j = 1:length(b)
        ab{i,j} = [a(i), b(j)];
        v = ones(N,1)*b(j);
        v(end) = a(i);
        v(end-1) = a(i);
        Newpole = sort(roots(Filter.Es),'ComparisonMethod','real') + (1/FBW/Q - 1/FBW/effectiveQ)*v;
        FilterPd.Es = poly(Newpole);
        FilterPd.Ps = Filter.Ps;

        Ns = 1000;
        for m = 1:Ns
            k(m) = abs(polyval(FilterPd.Es,-1.05i + 2.1i/(Ns-1)*(m-1))/polyval(FilterPd.Ps,-1.05i + 2.1i/(Ns-1)*(m-1)));
        end

        FilterPd.epsilon = 1/min(k);
%         FilterPd.epsilonR = 1;

        Fs2 = conv(FilterPd.Es,conj(FilterPd.Es).*((-1).^[N:-1:0])) - ...
            conv(FilterPd.Ps,conj(FilterPd.Ps).*((-1).^[N:-1:0]))/(FilterPd.epsilon*FilterPd.epsilon);
        RootF = roots(Fs2);
        RootF = RootF(real(RootF) < 0);
        FilterPd.Fs = poly(RootF);
        FilterPd.F22s = ((-1)^N)*conj(FilterPd.Fs).*(-1).^[N:-1:0];
        Now = abs(polyval(FilterPd.Ps,1i*w)./polyval(FilterPd.Es,1i*w)/FilterPd.epsilon);
        Now = Now/max(Now);
        After(i,j) = sum(abs(db(Now) - db(Before)))/Num;
        Gamma(i,j) = db(abs(abs(polyval(FilterPd.Fs,0)./polyval(FilterPd.Es,0))));
    end
end

Index = find(After<0.1);
[RL,X] = min(Gamma(Index));
ab_final = ab{Index(X)};

v = ones(N,1)*ab_final(2);
v(end) = ab_final(1);
v(end-1) = ab_final(1);
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