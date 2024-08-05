%% Characteristic Polynomials E, P, F to Y-parameters

function Y = EPF2Ylossy(Filter, K, Y0)

% EPF2Y aims to calculate the Y-parameters using the characteristic 
% polynomials (Es, Ps, Fs, Ep,Epr)

% Y0 is the characteristic admittance of the ports
% Y.Yd: the denominator polynomial of the Y-parameters
% Y.Y11n: the denominator polynomial of Y11
% Y.Y12n: the denominator polynomial of Y12
% Y.Y21n: the denominator polynomial of Y21
% Y.Y22n: the denominator polynomial of Y22
% Y11 = Y11n/Yd; Y21 = Y12 = Y12n/Yd; Y22 = Y22n/Yd; 

% By yellowbook, 2024-07-21
N = length(Filter.Es) - 1;
Filter.Ps = Filter.Ps./Filter.epsilon;
Filter.Fs = Filter.Fs./Filter.epsilonR;
Filter.F22s = Filter.F22s./Filter.epsilonR;
T = K*K*((-1)^N)*conj(Filter.Es).*(-1).^[N:-1:0];

Y.Yd = Filter.Es + Filter.Fs + Filter.F22s + T;
Y.Y11n = Y0*(Filter.Es - Filter.Fs + Filter.F22s - T);
Y.Y12n = Y0*(-2*Filter.Ps);
Y.Y21n = Y.Y12n;
Y.Y22n = Y0*(Filter.Es + Filter.Fs - Filter.F22s - T);

end
