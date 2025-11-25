%% To Box6degree Topology

function [theta_min,M2_min] = Box6dExtract(M1, Tz, Pull)

% TriExtract aims to transform the arrow-type matrix to the Box6degree-type
% coupling matrix.
% Pull is the distance the Box6degree section moves toward the S port.

% By yellowbook, 2025-11-16


    N = length(M1) - 2;

    % to CQ 24
    M1 = CTExtract(M1,Tz(1)*1i,Pull);
    M1 = CTExtract(M1,Tz(2)*1i,Pull-1);
    M1 = Rotate(M1, N-Pull-1, N-Pull+1, N-Pull, 'column');
    dtheta = 1;
    theta1 = [-90:dtheta:90]/180*pi;
   
    
    for i = 1:length(theta1)
        R = eye(N+2,N+2);
        R(N-Pull, N-Pull) = cos(theta1(i));
        R(N-Pull+1, N-Pull+1) = cos(theta1(i));
        R(N-Pull, N-Pull+1) = -sin(theta1(i));
        R(N-Pull+1, N-Pull) = sin(theta1(i));
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
%         M23(i) = M2(2,3);
        M34(i) = M2(N-Pull-1,N-Pull);
        M45(i) = M2(N-Pull,N-Pull+1);
        M56(i) = M2(N-Pull+1,N-Pull+2);
%         M67(i) = M2(6,7);
    end
    
    figure()

    plot(theta1/pi*180,abs(real(M34)),'linewidth',1.5);
    hold on    
    plot(theta1/pi*180,abs(real(M45)),'linewidth',1.5);
    hold on    
    plot(theta1/pi*180,abs(real(M56)),'linewidth',1.5);
    hold on
    legend('M23','M34','M45', 'Location', 'NorthWest')
    ylabel('Value');
    set(gca,'FontName','Times New Roman');
%     set(gca,'FontSize',18);
    set(gca,'linewidth',1.2);
    xlabel('Theta');
    set(gca,'FontName','Times New Roman');
%     set(gca,'FontSize',18);
    set(gca,'linewidth',1.2);
%     grid on

[~,locs34] = findpeaks(1./abs(real(M34)),theta1);
[~,locs45] = findpeaks(1./abs(real(M45)),theta1);
[~,locs56] = findpeaks(1./abs(real(M56)),theta1);

% threshold = 0.01;
% locs34 = find(abs(real(M34))<threshold);
% locs45 = find(abs(real(M45))<threshold);
% locs56 = find(abs(real(M56))<threshold);

% M34
if length(locs34)>=1
    for i = 1:length(locs34)
        [theta_min(i), M2_min(:,:,i)] = golden_section_search34(M1, Pull, ...
            locs34(i)-dtheta/180*pi, locs34(i)+dtheta/180*pi);
    end
end
% M45
if length(locs45)>=1
    for i = 1:length(locs45)
        [theta_min(i+length(locs34)), M2_min(:,:,i+length(locs34))] =...
            golden_section_search45(M1, Pull, ...
            locs45(i)-dtheta/180*pi, locs45(i)+dtheta/180*pi);
    end
end
% M56
if length(locs56)>=1
    for i = 1:length(locs56)
        [theta_min(i+length(locs34)+length(locs45)),...
            M2_min(:,:,i+length(locs34)+length(locs45))] =...
            golden_section_search56(M1, Pull, ...
            locs56(i)-dtheta/180*pi, locs56(i)+dtheta/180*pi);
    end
end
    % 
    % M34 => 0
    % 
    % switch 3 and 5
%     P = eye(N+2,N+2);
%     P(4, 4) = 0;
%     P(6, 6) = 0;
%     P(4, 6) = 1;
%     P(6, 4) = 1;
%     M2 = P*M2*P;


end