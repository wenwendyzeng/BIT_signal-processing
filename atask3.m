clc;
clear all;
%% ��������
% ��ʼ����
f0 = 10e9;      %��Ƶ10GHz
Tp = 10e-6;     %������10us
Tr = 100e-6;    %�����ظ�����100us
B  = 10e6;      %�źŴ���10MHz
fs = 100e6;     %������100MHz
R0 = 3000;      %Ŀ���ʼ����3km
c = 3e8;        %����
v = 60;         %�����ٶ�
m = 64;         %64Ϊ���������
SNR = 10;       %�����(dB)
A = 1;          %�źŷ���
% ��������
gate = c * Tp / 2;          %�ź��������߹��ľ���
Gate = 4500;                %���벨��
gate_total = gate + Gate;   %�ܲ��� = ���벨��+�������߹��ľ���
Number  = round(2*gate_total/c*fs);    %���������
echo_t  = linspace(0,2*gate_total/c,Number);  %���ų���
% tao = 2 * (R0 - m * Tr * v)/c;
% fd = 2 * f0 * v / c; 

%% ��������
%�������ʣ�SNR=10log10(A^2/sigmaf),A^2Ϊ�źŹ���
sigmaf = A^2/(10^(SNR/10));  
noise = sqrt(sigmaf/2)*(randn(1,Number)+1i*randn(1,Number));
%% �ź�ģ�͹���
% �����ź�
k = B / Tp; %���Ե�Ƶ��
t = 0:1/fs:Tp; %�źų���
st = rectpuls((t-Tp/2)/Tp).*exp(1i*pi*k*t.^2); %�ο��ź�

% �ز��ź�
for i = 1:1:m
    tao = 2 * (R0 - i * Tr * v)/c;
    %������
    sr(i,:) = rectpuls((echo_t-2*(R0-(i-1)*v*Tr)/c-Tp/2)/Tp).*exp(1i*pi*k*(echo_t-2*(R0-(i-1)*v*Tr)/c).^2-1i*2*pi*f0*round(2*R0/c*fs)+1i*2*pi*(2*f0*v/c)*(i-1)*Tr);
    %������
    sr_noise(i,:) = sr(i,:) + noise;
end
%% ����ѹ����Ƶ��ƥ�䣩
FFT_Number  = 2^nextpow2(length(t)+Number-1); %FFT����
S1  = fft(st,FFT_Number);   %�ο��ź�FFT�任
%��ʱ����
for i = 1:1:m
    S2_noise(i,:)  = fft(sr_noise(i,:),FFT_Number);   %�ز��ź�FFT�任
    s3_noise(i,:)  = ifft(S2_noise(i,:).*conj(S1));   %����ѹ��
    Z1_noise(i,:)  = abs(s3_noise(1:Number)); 
    Z_noise(i,:)   = Z1_noise(i,:)/max(Z1_noise(i,:));
    Z2_noise(i,:)  = 20*log10(Z_noise(i,:));          %db��ʽ
    
    S2(i,:)  = fft(sr(i,:),FFT_Number);   %�ز��ź�FFT�任
    s3(i,:)  = ifft(S2(i,:).*conj(S1));   %����ѹ��
    Z1(i,:)  = abs(s3(1:Number)); 
    Z(i,:)   = Z1(i,:)/max(Z1_noise(i,:));
    Z2(i,:)  = 20*log10(Z(i,:));          %db��ʽ
end

%% PD����
%�����ձ任����ʱ����
for ii = 1:1:Number
    Sb_noise(:,ii) = fft(s3_noise(:,ii));
    Sb_a_noise(:,ii) = fftshift(abs(Sb_noise(:,ii)));
    
    Sb(:,ii) = fft(s3(:,ii));
    Sb_a(:,ii) = fftshift(abs(Sb(:,ii)));
end

%% ���ķ�
[maxvalue_R_noise,maxposition1]  = max(Sb_a_noise,[],2);  
[maxvalue_V_noise,max_V_noise]   = max(maxvalue_R_noise);
max_R_noise                     = maxposition1(max_V_noise);

[maxvalue_R,maxposition2]  = max(Sb_a,[],2);  
[maxvalue_V,max_V]        = max(maxvalue_R);
max_R                     = maxposition2(max_V);
%������
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
% ������
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
disp(['Ŀ��ز�������ʱ�ľ���Ϊ: ', num2str(R),'m']);
disp(['Ŀ��ز�������ʱ���ٶ�Ϊ: ', num2str(V),'m/s']);
disp(['Ŀ��ز�������ʱ�ľ���Ϊ: ', num2str(R_noise),'m']);
disp(['Ŀ��ز�������ʱ���ٶ�Ϊ: ', num2str(V_noise),'m/s']);