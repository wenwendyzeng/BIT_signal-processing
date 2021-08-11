clear all;
clc;
%% 参数设置
Tc = 0.1e-6;   %码元宽度
Tp = 1.3e-6;%脉冲宽度
fs  = 20e6; %采样率
bark = [0,0,0,0,0,pi,pi,0,0,pi,0,pi,0];
m = 13;     %码元位数
N = fs * Tp*20;
% N = 13;
%% 巴克码信号模型
% t  = 0:1/fs:(N - 1)/fs;
t  = 0:1/fs:10e-6;
u_t = 0;
for n =0:1:m - 1
    u_t = u_t + rectpuls((t - Tc/2-n*Tc)/Tc)*exp(1i*bark(n + 1));

end

%% 模糊函数定义计算
i = 1;
% for fd = -fs/20:fs / N:fs/20-fs/N
for fd = -10e6/4: 10e6/4 / 800:10e6/4-10e6/4/800
    U = fft(u_t);
    U1 = U.*exp(1i*2*pi*fd*t);
    ka(i,:) = fftshift(ifft(u_t .*conj(ifft(U1))));
%     ka(i,:) = xcorr(u_t,ifft(U1));
    i = i + 1;
end
%% 作图
figure(1);
mesh(abs(ka));%三维网格图形绘制
% ylim([600,1000]);
% xlabel('\tau: us','Fontsize',12,'FontName','Times New Roman');
% ylabel('\xi: MHz','Fontsize',12,'FontName','Times New Roman');
% zlabel('归一化幅度','Fontsize',12,'FontName','宋体');
title('线性调频信号模糊函数三维图','Fontsize',12,'FontName','宋体');
figure(2);
contour(abs(ka))%等高线图形绘制
% ylim([600,1000]);
% xlabel('\tau: us','Fontsize',12,'FontName','Times New Roman');
% ylabel('\xi: MHz','Fontsize',12,'FontName','Times New Roman');
% zlabel('归一化幅度','Fontsize',12,'FontName','宋体');
grid on
title('线性调频信号等高线图','Fontsize',12,'FontName','宋体');