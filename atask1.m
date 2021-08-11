clc;
clear all;
%% 参数设置
% 初始参数
f0 = 10e9;  %载频10GHz
Tp = 10e-6; %脉冲宽度10us
B  = 10e6;  %信号带宽10MHz
fs = 100e6; %采样率100MHz
R0 = 3000;  %目标初始距离3km
c = 3e8;    %光速
% 后续参数
gate = c * Tp / 2; %信号脉宽内走过的距离
Gate = 4500;       %距离波门
gate_total = gate + Gate;  %总波门 = 距离波门+脉冲内走过的距离
Number  = round(2*gate_total/c*fs);    %采样点个数
echo_t  = linspace(0,2*gate_total/c,Number);  %波门长度
tao = 2 * R0/c;
%% 信号模型构建
% 发射信号
k = B / Tp; %线性调频率
t = 0:1/fs:Tp; %信号长度
st = rectpuls((t-Tp/2)/Tp).*exp(1i*pi*k*t.^2); %参考信号
% st = exp(1i*pi*k*t.^2); %参考信号

% 回波信号
sr = rectpuls((echo_t-tao-Tp/2)/Tp).*exp(1i*pi*k*(echo_t-tao).^2-1i*2*pi*f0*tao);
% sr = rectpuls((echo_t-tao)/Tp).*exp(1i*pi*k*(echo_t-tao).^2-1i*2*pi*f0*round(2*R0/c*fs));
%% 作图
figure(1);
subplot(211)
plot(t * c / 2,real(st))
% axis([0,3000,-1,1]);
xlabel('距离/m');
ylabel('实部');
title('参考信号I路');

subplot(212)
plot(t * c / 2,imag(st))
axis tight;
xlabel('距离/m');
ylabel('虚部');
title('参考信号Q路');

figure (2)
subplot(211)
plot(echo_t * c/2,real(sr))
axis([3000,5000,-2,2]);
xlabel('距离/m');
ylabel('实部');
title('回波信号I路');

subplot(212)
plot(echo_t * c/2,imag(sr))
axis([3000,5000,-2,2]);
xlabel('距离/m');
ylabel('虚部');
title('回波信号Q路');
%% 脉冲压缩（频域匹配）
FFT_Number  = 2^nextpow2(length(t)+Number-1); %FFT点数
S1  = fft(st,FFT_Number);   %参考信号FFT变换
S2  = fft(sr,FFT_Number);   %回波信号FFT变换
s3  = ifft(S2.*conj(S1));   %脉冲压缩
Z1  = abs(s3(1:Number)); 
Z   = Z1/max(Z1);
Z2  = 20*log10(Z);          %db形式

%% 频谱图
figure(3);
plot(linspace(-0.1,0.1,length(t)),abs(fftshift(fft(st))));
% plot(linspace(t * fs*(-0.5:0.5),length(t)),abs(fftshift(fft(st))));
xlabel('频率');
ylabel('幅度/V');
title('参考信号频谱图');

figure(4);
plot(linspace(-0.1,0.1,length(echo_t)),abs(fftshift(fft(sr))));
xlabel('频率');
ylabel('幅度/V');
title('回波信号频谱图');

figure(5)
plot(echo_t * c /2,Z2);
xlabel('距离/m');
axis([2600,3400,-40,0]);
ylabel('幅度/dB');
title('脉压后的信号图');