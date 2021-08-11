%% ��ҵ2��PD����
%%
clc;
clear all;
%%
%���������
f0       = 10e9;                          %��Ƶ
tau      = 10e-6;                         %������
B        = 10e6;                          %����
C        = 3e8;                           %����
fs       = 100e6;                         %����Ƶ��
R0       = 3000;                          %��ʼ����
T        = 100e-6;                        %�����ظ�����
CPI      = 64;                            %���������
v        = 60;                            %Ŀ���ٶ�
Gate_echo= 4500;                          %���벨��
%%
%�źŹ���
Gate_signal = C*tau/2;                       %�ź������Ӧ����
Gate        = Gate_echo+Gate_signal;         %���벨�ż�����
Gate_Number = round(2*Gate/C*fs);            %�����ڲ��������
t           = 0:1/fs:tau;                    %�źų���
echo_t      = linspace(0,2*Gate/C,Gate_Number);  %���ų���
K           = B/tau;                         %��Ƶϵ��
% transmit_signal = exp(1i*pi*K*t.^2);             %�����ź�
transmit_signal = rectpuls((t-T/2)/T).*exp(1i*pi*K*t.^2); 
%�ز��ź�
for m=1:1:CPI
echo_signal(m,:) = rectpuls((echo_t-2*(R0-(m-1)*v*T)/C-tau/2)/(tau)).*exp(1i*pi*K*(echo_t-2*(R0-(m-1)*v*T)/C).^2-1i*2*pi*f0*round(2*R0/C*fs)+1i*2*pi*(2*f0*v/C)*(m-1)*T)+sqrt(0.1)*(randn(1,Gate_Number)+1i*randn(1,Gate_Number));
end
%% �źŴ���
% ��ѹ
FFT_Number  = 2^nextpow2(length(t)+Gate_Number-1); %FFT����
Srw  = fft(transmit_signal,FFT_Number);
for m=1:1:CPI
Sw(m,:)    = fft(echo_signal(m,:),FFT_Number);   %����FFT
Sot(m,:)   = ifft(Sw(m,:).*conj(Srw));
Z1(m,:)    = abs(Sot(m,(1:Gate_Number)));
Z(m,:)     = Z1(m,:)/max(Z1(m,:));
Z(m,:)     = 20*log10(Z(m,:));    %db��ʽ
[maxvalue,maxposition] = max(Z(m,:));         
end
% FFT
for fm=1:1:Gate_Number
Dop(:,fm)   = fft(Sot(:,fm)); %������FFT
A_Dop(:,fm) = fftshift(abs(Dop(:,fm)));
end
%�󼫴�ֵ��Ӧ������
[maxvalue,max_V]   = max(A_Dop(:,maxposition));
%PD����
fd                 = (max_V-33)/CPI/T;
V_pd               = fd*C/2/f0
%% ������ٷ�Χ����پ���
%���ٷ�Χ
fd_max             = 1/T/2;
V_max              = fd_max*C/2/f0
%�ٶȷֱ���
det_fd             = 1/T/64;
det_v              = det_fd*C/2/f0
%% ��ͼ
figure(1)
mesh(echo_t*C/2,linspace(-75,75,64),A_Dop);
axis tight;
xlabel('���룺m');
ylabel('�ٶȣ�m/s');
title('��ά����-������ƽ��');