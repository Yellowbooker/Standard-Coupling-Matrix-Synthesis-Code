%% Standard Coupling Matrix Synthesis Code
% For General Chebyshev Filter

clear
close all

% N = 4; % order
% Tz = [-3.7431*1i -1.8051*1i 1.5699*1i 6.1910*1i];  % Tz in s-domain
% Tz = Tz/1i; % in w-domain
% RL = 22; % RL level within the passband

% N = 8; % order
% Tz = [1.1185i 1.2973i];  % Tz in s-domain
% Tz = Tz/1i; % in w-domain
% RL = 23; % RL level within the passband

N = 10; % order
Tz = [-1.2i -1.5i -1.1i 1.3i];  % Tz in s-domain
Tz = Tz/1i; % in w-domain
RL = 23; % RL level within the passband


[Fw, Pw] = General_Chebyshev(Tz, N);

w = [-8:0.01:8]; % w(omega) range in the figures
figure('name','General_Chebyshev');
plot(w,db(polyval(Pw,w)./polyval(Fw,w)),'Linewidth',2);
ylabel('Magnitude (dB)','fontsize',18);
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',18);
set(gca,'linewidth',1.2);
xlabel('\omega (rad/s)','fontsize',18);
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',18);
set(gca,'linewidth',1.2);
grid on

Filter = Cheby2EPF(Fw, Pw, RL);
[Maxerror,S_poly] = CheckUnitary(Filter, w);

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
grid on

Y = EPF2Y(Filter, 1); % Y0 is 1S
M_trans = Y2CMtrans(Y);

G = zeros(N+2,N+2);
G(1,1) = 1;
G(N+2,N+2) = 1;
C = eye(N+2,N+2);
C(1,1) = 0;
C(N+2,N+2) = 0;

% CM = to_foldedCM(N,M_trans);
CM = CM2arrow(M_trans);
% CM = BoxExtract(CM,Tz(1)*1i,3);
% CM = BoxExtract(CM,Tz(2)*1i,0);
% CM = TriExtract(CM,Tz(1)*1i,2);
% CM = TriExtract(CM,Tz(2)*1i,1);
CM = CQExtract(CM, [Tz(1),Tz(end)], 7, '13');
% [1,2]
% R = eye(N+2,N+2);
% theta_r = atan((CM(3, 5)-CM(4, 6))/(CM(4, 5) + CM(3, 6)));
% R(3, 3) = cos(theta_r);
% R(4, 4) = cos(theta_r);
% R(3, 4) = -sin(theta_r);
% R(4, 3) = sin(theta_r);
% CM = R*CM*R.';
% [2,3]
% R = eye(N+2,N+2);
% theta_r = atan((CM(4, 6)-CM(3, 5))/(CM(3, 4) + CM(5, 6)));
% R(4, 4) = cos(theta_r);
% R(5, 5) = cos(theta_r);
% R(4, 5) = -sin(theta_r);
% R(5, 4) = sin(theta_r);
% CM = R*CM*R.';
% theta1 = [0:0.01:360];
% for i = 1:length(theta1)
%     CM1(:,:,i) = CQ2EBox6(CM,theta1(i));
% end
% M36 = abs(CM1(3,6,:));
% CM = CQExtract(CM, Tz, 4, '24');
% CM = CQ2EBox6(CM,0.463/180*pi);
% CM = CQ2EBox6(CM,136.53/180*pi);
% CM = CQ2EBox6(CM,58.82/180*pi);
% CM = CQ2EBox6(CM,85.43/180*pi);
% CM = CQ2EBox6(CM,136.6/180*pi);
% CM = CQ2EBox6(CM,180.462/180*pi); %=0.463
% CM = CQ2EBox6(CM,238.8/180*pi); %=58.82
% CM = CQ2EBox6(CM,265.4/180*pi); %=85.43
% CM = CQ2EBox6(CM,316.6/180*pi); %=136.6
Pull = 2;
[~,Bx6CM] = Box6dExtract(CM, [Tz(2),Tz(3)], Pull);
k = 0;
for i = 1:size(Bx6CM,3)
    if abs(Bx6CM(N-Pull-1,N-Pull+2,i)) < 1e-2
        k = k + 1;
        Bx6CMs(:,:,k) = NormalizeCM(Bx6CM(:,:,i));
    end
end
real(Bx6CMs)
CM = Bx6CMs(:,:,1);

[S11, S12, S22] = CMFC_Response(CM, C, G, w);
figure('name','S-parameters_couplingmatrix');
plot(w,db(S11),w,db(S12),'Linewidth',2);
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
