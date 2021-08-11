clc;
clear all;
%% 参数设置
% 初始参数
f0 = 10e9;      %载频10GHz
Tp = 10e-6;     %脉冲宽度10us
Tr = 100e-6;    %脉冲重复周期100us
B  = 10e6;      %信号带宽10MHz
fs = 100e6;     %采样率100MHz
R0 = 3000;      %目标初始距离3km
c = 3e8;        %光速
v = 60;         %径向速度
m = 64;         %64为总脉冲个数
SNR = 10;       %信噪比(dB)
A = 1;          %信号幅度
% 后续参数
gate = c * Tp / 2;          %信号脉宽内走过的距离
Gate = 4500;                %距离波门
gate_total = gate + Gate;   %总波门 = 距离波门+脉冲内走过的距离
Number  = round(2*gate_total/c*fs);    %采样点个数
echo_t  = linspace(0,2*gate_total/c,Number);  %波门长度
% tao = 2 * (R0 - m * Tr * v)/c;
% fd = 2 * f0 * v / c; 

%% 噪声构建
%噪声功率，SNR=10log10(A^2/sigmaf),A^2为信号功率
sigmaf = A^2/(10^(SNR/10));  
noise = sqrt(sigmaf/2)*(randn(1,Number)+1i*randn(1,Number));
%% 信号模型构建
% 发射信号
k = B / Tp; %线性调频率
t = 0:1/fs:Tp; %信号长度
st = rectpuls((t-Tp/2)/Tp).*exp(1i*pi*k*t.^2); %参考信号

% 回波信号
for i = 1:1:m
    tao = 2 * (R0 - i * Tr * v)/c;
    %无噪声
    sr(i,:) = rectpuls((echo_t-2*(R0-(i-1)*v*Tr)/c-Tp/2)/Tp).*exp(1i*pi*k*(echo_t-2*(R0-(i-1)*v*Tr)/c).^2-1i*2*pi*f0*round(2*R0/c*fs)+1i*2*pi*(2*f0*v/c)*(i-1)*Tr);
    %有噪声
    sr_noise(i,:) = sr(i,:) + noise;
end
%% 脉冲压缩（频域匹配）
FFT_Number  = 2^nextpow2(length(t)+Number-1); %FFT点数
S1  = fft(st,FFT_Number);   %参考信号FFT变换
%快时间域
for i = 1:1:m
    S2_noise(i,:)  = fft(sr_noise(i,:),FFT_Number);   %回波信号FFT变换
    s3_noise(i,:)  = ifft(S2_noise(i,:).*conj(S1));   %脉冲压缩
    Z1_noise(i,:)  = abs(s3_noise(1:Number)); 
    Z_noise(i,:)   = Z1_noise(i,:)/max(Z1_noise(i,:));
    Z2_noise(i,:)  = 20*log10(Z_noise(i,:));          %db形式
    
    S2(i,:)  = fft(sr(i,:),FFT_Number);   %回波信号FFT变换
    s3(i,:)  = ifft(S2(i,:).*conj(S1));   %脉冲压缩
    Z1(i,:)  = abs(s3(1:Number)); 
    Z(i,:)   = Z1(i,:)/max(Z1_noise(i,:));
    Z2(i,:)  = 20*log10(Z(i,:));          %db形式
end

%% PD测速
%多普勒变换（慢时间域）
for ii = 1:1:Number
    Sb_noise(:,ii) = fft(s3_noise(:,ii));
    Sb_a_noise(:,ii) = fftshift(abs(Sb_noise(:,ii)));
    
    Sb(:,ii) = fft(s3(:,ii));
    Sb_a(:,ii) = fftshift(abs(Sb(:,ii)));
end

%% 形心法
[maxvalue_R_noise,maxposition1]  = max(Sb_a_noise,[],2);  
[maxvalue_V_noise,max_V_noise]   = max(maxvalue_R_noise);
max_R_noise                     = maxposition1(max_V_noise);

[maxvalue_R,maxposition2]  = max(Sb_a,[],2);  
[maxvalue_V,max_V]        = max(maxvalue_R);
max_R                     = maxposition2(max_V);
%有噪声
R_noise = c/2/fs*(Sb_a_noise(max_V_noise,max_R_noise-3)*(max_R_noise-3)+...
    Sb_a_noise(max_V_noise,max_R_noise-2)*(max_R_noise-2)+...
    Sb_a_noise(max_V_noise,max_R_noise-1)*(max_R_noise-1)+...
    Sb_a_noise(max_V_noise,max_R_noise)*(max_R_noise)+...
    Sb_a_noise(max_V_noise,max_R_noise+1)*(max_R_noise+1)+...
    Sb_a_noise(max_V_noise,max_R_noise+2)*(max_R_noise+2)+...
    Sb_a_noise(max_V_noise,max_R_noise+3)*(max_R_noise+3))...
    /sum(Sb_a_noise(max_V,max_R_noise-3:max_R_noise+3));
V_noise           = ...
   (((Sb_a_noise(max_V-3,max_R_noise)*(max_V_noise-3)+...
   Sb_a_noise(max_V_noise-2,max_R_noise)*(max_V_noise-2)+...
   Sb_a_noise(max_V_noise-1,max_R_noise)*(max_V_noise-1)+...
   Sb_a_noise(max_V_noise,max_R_noise)*(max_V_noise)+...
   Sb_a_noise(max_V_noise+1,max_R_noise)*(max_V_noise+1)+...
   Sb_a_noise(max_V_noise+2,max_R_noise)*(max_V_noise+2)+...
   Sb_a_noise(max_V_noise+3,max_R_noise)*(max_V_noise+3))...
   /sum(Sb_a_noise(max_V_noise-3:max_V_noise+3,max_R_noise)))-33)/m/Tr*c/2/f0;
% 无噪声
R = c/2/fs*(Sb_a(max_V,max_R-3)*(max_R-3)+...
    Sb_a(max_V,max_R-2)*(max_R-2)+...
    Sb_a(max_V,max_R-1)*(max_R-1)+...
    Sb_a(max_V,max_R)*(max_R)+...
    Sb_a(max_V,max_R+1)*(max_R+1)+...
    Sb_a(max_V,max_R+2)*(max_R+2)+...
    Sb_a(max_V,max_R+3)*(max_R+3))...
    /sum(Sb_a(max_V,max_R-3:max_R+3));
V           = ...
   (((Sb_a(max_V-3,max_R)*(max_V-3)+...
   Sb_a(max_V-2,max_R)*(max_V-2)+...
   Sb_a(max_V-1,max_R)*(max_V-1)+...
   Sb_a(max_V,max_R)*(max_V)+...
   Sb_a(max_V+1,max_R)*(max_V+1)+...
   Sb_a(max_V+2,max_R)*(max_V+2)+...
   Sb_a(max_V+3,max_R)*(max_V+3))...
   /sum(Sb_a(max_V-3:max_V+3,max_R)))-33)/m/Tr*c/2/f0;
disp(['目标回波无噪声时的距离为: ', num2str(R),'m']);
disp(['目标回波无噪声时的速度为: ', num2str(V),'m/s']);
disp(['目标回波有噪声时的距离为: ', num2str(R_noise),'m']);
disp(['目标回波有噪声时的速度为: ', num2str(V_noise),'m/s']);