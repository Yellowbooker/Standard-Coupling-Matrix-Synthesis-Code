%% Predistortion Coupling Matrix Synthesis Code
% For Predistortion Chebyshev Filter

clear
close all

% N = 7; % order
% Tz = [-3.7431*1i -1.8051*1i 1.5699*1i 6.1910*1i];  % Tz in s-domain
% Tz = Tz/1i; % in w-domain
% RL = 22; % RL level within the passband

N = 10; % order
Tz = [-1.0912i 1.0912i -1.605389i 1.605389i 0.6173+0.34881i -0.6173+0.34881i...
    0.6173-0.34881i -0.6173-0.34881i];  % Tz in s-domain
Tz = Tz/1i; % in w-domain
RL = 22; % RL level within the passband

Q = 2000; % Q of the resonators
FBW = 0.01;
% alpha = 1/FBW/Q;
effectiveQ = 20000;

[Fw, Pw] = General_Chebyshev(Tz, N);

Filter = Cheby2EPF(Fw, Pw, RL);

% FilterPd = Predistortion(Filter, Q, FBW, effectiveQ);
[FilterPd,ab_final] = AdaptivePredistortion(Filter, Q, FBW, effectiveQ);
% display(ab_final);

w = [-5:0.005:5]; % w(omega) range in the figures
[Maxerror,S_poly] = CheckUnitary(FilterPd, w);

figure('name','S-parameters_polynomials');
plot(w,db(polyval(FilterPd.Fs,w*1i)./polyval(FilterPd.Es,w*1i)./FilterPd.epsilonR)...
    ,w,db(polyval(FilterPd.Ps,w*1i)./polyval(FilterPd.Es,w*1i)./FilterPd.epsilon),'Linewidth',2);
ylabel('S-parameters (dB)','fontsize',18);
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',18);
set(gca,'linewidth',1.2);
xlabel('\omega (rad/s)','fontsize',18);
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',18);
set(gca,'linewidth',1.2);
ylim([-100,0]);
grid on

YPd = EPF2Y(FilterPd, 1); % Y0 is 1S
M_transPd = Y2CMtrans(YPd);
Y = EPF2Y(Filter, 1); % Y0 is 1S
M_trans = Y2CMtrans(Y);

G = zeros(N+2,N+2);
G(1,1) = 1;
G(N+2,N+2) = 1;
G(2:end-1,2:end-1) = eye(N)*(1/FBW/Q);
C = eye(N+2,N+2);
C(1,1) = 0;
C(N+2,N+2) = 0;

CMPd = to_foldedCM(N,M_transPd);
CM = to_foldedCM(N,M_trans);
display(real(CMPd));
% display(imag(CMPd));
[S11Pd, S12Pd, S22Pd] = CMFC_Response(CMPd, C, G, w);
[S11, S12, S22] = CMFC_Response(CM, C, G, w);
figure('name','S-parameters_couplingmatrixPd');
plot(w,db(S11Pd),w,db(S12Pd),w,db(S22Pd),'Linewidth',2);
ylabel('S-parameters (dB)','fontsize',18);
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',18);
set(gca,'linewidth',1.2);
xlabel('\omega (rad/s)','fontsize',18);
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',18);
set(gca,'linewidth',1.2);
% ylim([-100,0]);
grid on

figure('name','S-parameters_couplingmatrix');
plot(w,db(S12Pd) - max(db(S12Pd)),w,db(S12) - max(db(S12)),':','Linewidth',2);
ylabel('Magnitude (dB)','fontsize',18);
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',18);
set(gca,'linewidth',1.2);
xlabel('\omega (rad/s)','fontsize',18);
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',18);
set(gca,'linewidth',1.2);
ylim([-5,0]);
grid on

