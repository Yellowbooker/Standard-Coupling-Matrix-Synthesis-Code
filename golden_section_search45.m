function [x45_min, M2_min45] = golden_section_search45(M1, Pull, a, b)
% 0.618法一维搜索核心代码
% 输入：目标函数f, 搜索区间[a,b], 容差tol
% 输出：极小值点x_min, 极小值f_min
% By Deepseek
tol = 1e-6;
N = length(M1) - 2;
%----M45
    rho = (sqrt(5)-1)/2;  % 黄金分割比0.618
    
    x1 = a + (1-rho)*(b-a);
    x2 = a + rho*(b-a);
    M2 = AfterRotate(M1, Pull, x1);
    f1 = abs(M2(N-Pull,N-Pull+1));
    M2 = AfterRotate(M1, Pull, x2);
    f2 = abs(M2(N-Pull,N-Pull+1));
    
    while (b - a) > tol
        if f1 > f2
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + rho*(b-a);
            M2 = AfterRotate(M1, Pull, x2);
            f2 = abs(M2(N-Pull,N-Pull+1));
        else
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + (1-rho)*(b-a);
            M2 = AfterRotate(M1, Pull, x1);
            f1 = abs(M2(N-Pull,N-Pull+1));
        end
    end
    
    x45_min = (a + b) / 2;
    M2_min45 = AfterRotate(M1, Pull, x45_min);
    P = eye(N+2,N+2);
    P(N-Pull-1, N-Pull-1) = 0;
    P(N-Pull+1, N-Pull+1) = 0;
    P(N-Pull-1, N-Pull+1) = 1;
    P(N-Pull+1, N-Pull-1) = 1;
    P(N-Pull, N-Pull) = 0;
    P(N-Pull+2, N-Pull+2) = 0;
    P(N-Pull, N-Pull+2) = 1;
    P(N-Pull+2, N-Pull) = 1;
    M2_min45 = P*M2_min45*P;
    
end