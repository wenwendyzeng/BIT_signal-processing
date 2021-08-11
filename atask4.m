% LFM�źźͰͿ����źŵ�ģ����������
clc;
clear all;
%% ��������
%��ʼ����
Tp = 10e-6;     %LFM�ź�������10us
B = 10e6;       %����
fs = 20e6;      %������
k = B/Tp;       %��Ƶ��
N = fs * Tp;   %���ڲ���
% N = 2^(nextpow2(Nr)+1);
%% LFM�����ź�ģ��
t = 0:1/fs:(N- 1)/fs;
u_t = rectpuls((t - Tp)/Tp).*exp(1i*pi*k*(t - Tp/2).^2);

%% ģ�������������
i = 1;
% t = 0:1/fs:(Nr - 1)/fs;
for fd = -fs/2:fs / N:fs/2-fs/N
    U = fft(u_t);
    U1 = U.*exp(1i*2*pi*fd*t);
    ka(i,:) = fftshift(ifft(u_t .*conj(ifft(U1))));
%     ka(i,:) = xcorr(u_t,ifft(U1));
    i = i + 1;
end

%% ��ͼ
figure(1)
tao = (-N/2:(N/2-1))*fs/1e6;
fd1 = (-fs/2 : fs/N : fs/2)/1e6;
[t, fd1] = meshgrid(fd,tao);
fd = (-fs/2 : fs/N : fs/2-fs/N)/1e6;
mesh(tao,fd,abs(ka));%��ά����ͼ�λ���
xlabel('\tau: us','Fontsize',12,'FontName','Times New Roman');
ylabel('\xi: MHz','Fontsize',12,'FontName','Times New Roman');
zlabel('��һ������','Fontsize',12,'FontName','����');
title('���Ե�Ƶ�ź�ģ������ͼ','Fontsize',12,'FontName','����');
figure(2)
contour(tao,fd1,abs(ka))%�ȸ���ͼ�λ���
xlabel('\tau: us','Fontsize',12,'FontName','Times New Roman');
ylabel('\xi: MHz','Fontsize',12,'FontName','Times New Roman');
zlabel('��һ������','Fontsize',12,'FontName','����');
grid on
title('���Ե�Ƶ�źŵȸ���ͼ','Fontsize',12,'FontName','����');