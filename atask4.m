% LFM信号和巴克码信号的模糊函数仿真
clc;
clear all;
%% 参数设置
%初始参数
Tp = 10e-6;     %LFM信号脉冲宽度10us
B = 10e6;       %带宽
fs = 20e6;      %采样率
k = B/Tp;       %调频率
N = fs * Tp;   %脉内采样
% N = 2^(nextpow2(Nr)+1);
%% LFM基带信号模型
t = 0:1/fs:(N- 1)/fs;
u_t = rectpuls((t - Tp)/Tp).*exp(1i*pi*k*(t - Tp/2).^2);

%% 模糊函数定义计算
i = 1;
% t = 0:1/fs:(Nr - 1)/fs;
for fd = -fs/2:fs / N:fs/2-fs/N
    U = fft(u_t);
    U1 = U.*exp(1i*2*pi*fd*t);
    ka(i,:) = fftshift(ifft(u_t .*conj(ifft(U1))));
%     ka(i,:) = xcorr(u_t,ifft(U1));
    i = i + 1;
end

%% 绘图
figure(1)
tao = (-N/2:(N/2-1))*fs/1e6;
fd1 = (-fs/2 : fs/N : fs/2)/1e6;
[t, fd1] = meshgrid(fd,tao);
fd = (-fs/2 : fs/N : fs/2-fs/N)/1e6;
mesh(tao,fd,abs(ka));%三维网格图形绘制
xlabel('\tau: us','Fontsize',12,'FontName','Times New Roman');
ylabel('\xi: MHz','Fontsize',12,'FontName','Times New Roman');
zlabel('归一化幅度','Fontsize',12,'FontName','宋体');
title('线性调频信号模糊函数图','Fontsize',12,'FontName','宋体');
figure(2)
contour(tao,fd1,abs(ka))%等高线图形绘制
xlabel('\tau: us','Fontsize',12,'FontName','Times New Roman');
ylabel('\xi: MHz','Fontsize',12,'FontName','Times New Roman');
zlabel('归一化幅度','Fontsize',12,'FontName','宋体');
grid on
title('线性调频信号等高线图','Fontsize',12,'FontName','宋体');