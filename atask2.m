%% 作业2：PD测速
%%
clc;
clear all;
%%
%参数定义段
f0       = 10e9;                          %载频
tau      = 10e-6;                         %脉冲宽度
B        = 10e6;                          %带宽
C        = 3e8;                           %光速
fs       = 100e6;                         %采样频率
R0       = 3000;                          %初始距离
T        = 100e-6;                        %脉冲重复周期
CPI      = 64;                            %总脉冲个数
v        = 60;                            %目标速度
Gate_echo= 4500;                          %距离波门
%%
%信号构建
Gate_signal = C*tau/2;                       %信号脉宽对应距离
Gate        = Gate_echo+Gate_signal;         %距离波门加脉宽
Gate_Number = round(2*Gate/C*fs);            %波门内采样点个数
t           = 0:1/fs:tau;                    %信号长度
echo_t      = linspace(0,2*Gate/C,Gate_Number);  %波门长度
K           = B/tau;                         %调频系数
% transmit_signal = exp(1i*pi*K*t.^2);             %发射信号
transmit_signal = rectpuls((t-T/2)/T).*exp(1i*pi*K*t.^2); 
%回波信号
for m=1:1:CPI
echo_signal(m,:) = rectpuls((echo_t-2*(R0-(m-1)*v*T)/C-tau/2)/(tau)).*exp(1i*pi*K*(echo_t-2*(R0-(m-1)*v*T)/C).^2-1i*2*pi*f0*round(2*R0/C*fs)+1i*2*pi*(2*f0*v/C)*(m-1)*T)+sqrt(0.1)*(randn(1,Gate_Number)+1i*randn(1,Gate_Number));
end
%% 信号处理
% 脉压
FFT_Number  = 2^nextpow2(length(t)+Gate_Number-1); %FFT点数
Srw  = fft(transmit_signal,FFT_Number);
for m=1:1:CPI
Sw(m,:)    = fft(echo_signal(m,:),FFT_Number);   %距离FFT
Sot(m,:)   = ifft(Sw(m,:).*conj(Srw));
Z1(m,:)    = abs(Sot(m,(1:Gate_Number)));
Z(m,:)     = Z1(m,:)/max(Z1(m,:));
Z(m,:)     = 20*log10(Z(m,:));    %db形式
[maxvalue,maxposition] = max(Z(m,:));         
end
% FFT
for fm=1:1:Gate_Number
Dop(:,fm)   = fft(Sot(:,fm)); %多普勒FFT
A_Dop(:,fm) = fftshift(abs(Dop(:,fm)));
end
%求极大值对应的坐标
[maxvalue,max_V]   = max(A_Dop(:,maxposition));
%PD测速
fd                 = (max_V-33)/CPI/T;
V_pd               = fd*C/2/f0
%% 计算测速范围与测速精度
%测速范围
fd_max             = 1/T/2;
V_max              = fd_max*C/2/f0
%速度分辨率
det_fd             = 1/T/64;
det_v              = det_fd*C/2/f0
%% 绘图
figure(1)
mesh(echo_t*C/2,linspace(-75,75,64),A_Dop);
axis tight;
xlabel('距离：m');
ylabel('速度：m/s');
title('二维距离-多普勒平面');