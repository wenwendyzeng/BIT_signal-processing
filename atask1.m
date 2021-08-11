clc;
clear all;
%% ��������
% ��ʼ����
f0 = 10e9;  %��Ƶ10GHz
Tp = 10e-6; %������10us
B  = 10e6;  %�źŴ���10MHz
fs = 100e6; %������100MHz
R0 = 3000;  %Ŀ���ʼ����3km
c = 3e8;    %����
% ��������
gate = c * Tp / 2; %�ź��������߹��ľ���
Gate = 4500;       %���벨��
gate_total = gate + Gate;  %�ܲ��� = ���벨��+�������߹��ľ���
Number  = round(2*gate_total/c*fs);    %���������
echo_t  = linspace(0,2*gate_total/c,Number);  %���ų���
tao = 2 * R0/c;
%% �ź�ģ�͹���
% �����ź�
k = B / Tp; %���Ե�Ƶ��
t = 0:1/fs:Tp; %�źų���
st = rectpuls((t-Tp/2)/Tp).*exp(1i*pi*k*t.^2); %�ο��ź�
% st = exp(1i*pi*k*t.^2); %�ο��ź�

% �ز��ź�
sr = rectpuls((echo_t-tao-Tp/2)/Tp).*exp(1i*pi*k*(echo_t-tao).^2-1i*2*pi*f0*tao);
% sr = rectpuls((echo_t-tao)/Tp).*exp(1i*pi*k*(echo_t-tao).^2-1i*2*pi*f0*round(2*R0/c*fs));
%% ��ͼ
figure(1);
subplot(211)
plot(t * c / 2,real(st))
% axis([0,3000,-1,1]);
xlabel('����/m');
ylabel('ʵ��');
title('�ο��ź�I·');

subplot(212)
plot(t * c / 2,imag(st))
axis tight;
xlabel('����/m');
ylabel('�鲿');
title('�ο��ź�Q·');

figure (2)
subplot(211)
plot(echo_t * c/2,real(sr))
axis([3000,5000,-2,2]);
xlabel('����/m');
ylabel('ʵ��');
title('�ز��ź�I·');

subplot(212)
plot(echo_t * c/2,imag(sr))
axis([3000,5000,-2,2]);
xlabel('����/m');
ylabel('�鲿');
title('�ز��ź�Q·');
%% ����ѹ����Ƶ��ƥ�䣩
FFT_Number  = 2^nextpow2(length(t)+Number-1); %FFT����
S1  = fft(st,FFT_Number);   %�ο��ź�FFT�任
S2  = fft(sr,FFT_Number);   %�ز��ź�FFT�任
s3  = ifft(S2.*conj(S1));   %����ѹ��
Z1  = abs(s3(1:Number)); 
Z   = Z1/max(Z1);
Z2  = 20*log10(Z);          %db��ʽ

%% Ƶ��ͼ
figure(3);
plot(linspace(-0.1,0.1,length(t)),abs(fftshift(fft(st))));
% plot(linspace(t * fs*(-0.5:0.5),length(t)),abs(fftshift(fft(st))));
xlabel('Ƶ��');
ylabel('����/V');
title('�ο��ź�Ƶ��ͼ');

figure(4);
plot(linspace(-0.1,0.1,length(echo_t)),abs(fftshift(fft(sr))));
xlabel('Ƶ��');
ylabel('����/V');
title('�ز��ź�Ƶ��ͼ');

figure(5)
plot(echo_t * c /2,Z2);
xlabel('����/m');
axis([2600,3400,-40,0]);
ylabel('����/dB');
title('��ѹ����ź�ͼ');