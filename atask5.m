clear all;
clc;
%% ��������
Tc = 0.1e-6;   %��Ԫ���
Tp = 1.3e-6;%������
fs  = 20e6; %������
bark = [0,0,0,0,0,pi,pi,0,0,pi,0,pi,0];
m = 13;     %��Ԫλ��
N = fs * Tp*20;
% N = 13;
%% �Ϳ����ź�ģ��
% t  = 0:1/fs:(N - 1)/fs;
t  = 0:1/fs:10e-6;
u_t = 0;
for n =0:1:m - 1
    u_t = u_t + rectpuls((t - Tc/2-n*Tc)/Tc)*exp(1i*bark(n + 1));

end

%% ģ�������������
i = 1;
% for fd = -fs/20:fs / N:fs/20-fs/N
for fd = -10e6/4: 10e6/4 / 800:10e6/4-10e6/4/800
    U = fft(u_t);
    U1 = U.*exp(1i*2*pi*fd*t);
    ka(i,:) = fftshift(ifft(u_t .*conj(ifft(U1))));
%     ka(i,:) = xcorr(u_t,ifft(U1));
    i = i + 1;
end
%% ��ͼ
figure(1);
mesh(abs(ka));%��ά����ͼ�λ���
% ylim([600,1000]);
% xlabel('\tau: us','Fontsize',12,'FontName','Times New Roman');
% ylabel('\xi: MHz','Fontsize',12,'FontName','Times New Roman');
% zlabel('��һ������','Fontsize',12,'FontName','����');
title('���Ե�Ƶ�ź�ģ��������άͼ','Fontsize',12,'FontName','����');
figure(2);
contour(abs(ka))%�ȸ���ͼ�λ���
% ylim([600,1000]);
% xlabel('\tau: us','Fontsize',12,'FontName','Times New Roman');
% ylabel('\xi: MHz','Fontsize',12,'FontName','Times New Roman');
% zlabel('��һ������','Fontsize',12,'FontName','����');
grid on
title('���Ե�Ƶ�źŵȸ���ͼ','Fontsize',12,'FontName','����');