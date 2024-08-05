%% Lossy Coupling Matrix Synthesis Code
% For Lossy General Chebyshev Filter

clear
close all

% N = 4; % order
% Tz = [-3.7431*1i -1.8051*1i 1.5699*1i 6.1910*1i];  % Tz in s-domain
% Tz = Tz/1i; % in w-domain
% RL = 22; % RL level within the passband

N = 10; % order
Tz = [-1.22*1i 1.22*1i -1.74*1i 1.74*1i 0.6+0.3i -0.6+0.3i 0.6-0.3i -0.6-0.3i];  % Tz in s-domain
Tz = Tz/1i; % in w-domain
RL = 22; % RL level within the passband

[Fw, Pw] = General_Chebyshev(Tz, N);
Filter = Cheby2EPF(Fw, Pw, RL);

% Attenuation of S21
attenuation = 30; % dB
K = 10^(-1*attenuation/20);
alpha = 1; % K <= alpha <= 1/K

Filter.Ps = Filter.Ps*K;
Filter.Fs = Filter.Fs*K*alpha;
Filter.F22s = Filter.F22s*K/alpha;

w = [-8:0.01:8]; % w(omega) range in the figures
figure('name','S-parameters_polynomials');
plot(w,db(polyval(Filter.Fs,w*1i)./polyval(Filter.Es,w*1i)./Filter.epsilonR)...
    ,w,db(polyval(Filter.Ps,w*1i)./polyval(Filter.Es,w*1i)./Filter.epsilon),'Linewidth',2);
ylabel('S-parameters (dB)','fontsize',18);
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',18);
set(gca,'linewidth',1.2);
xlabel('\omega (rad/s)','fontsize',18);
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',18);
set(gca,'linewidth',1.2);
ylim([-120,0]);
grid on

Y = EPF2Ylossy(Filter, K, 1); % Y0 is 1S
M_trans = Y2CMtrans(Y);

G = zeros(N+2,N+2);
G(1,1) = 1;
G(N+2,N+2) = 1;
C = eye(N+2,N+2);
C(1,1) = 0;
C(N+2,N+2) = 0;

% CM = to_foldedCM(N,M_trans);
CM = LossyCM2PCP(M_trans);
% imag(CM)
CM = LossDistribution(CM);
% imag(CM)
[S11, S12, S22] = CMFC_Response(M_trans, C, G, w);
figure('name','S-parameters_couplingmatrix');
plot(w,db(S11),w,db(S12),w,db(S22),'--','Linewidth',2);
ylabel('S-parameters (dB)','fontsize',18);
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',18);
set(gca,'linewidth',1.2);
xlabel('\omega (rad/s)','fontsize',18);
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',18);
set(gca,'linewidth',1.2);
ylim([-120,0]);
grid on