%�������3�����Ͼ����㶯�ز��źţ�����RD����
clc;
clear all;

%% ��������
c = 3*10^8;     %����
fc = 5.3*10^9;  %��Ƶ5.3GHz
lamda = c/fc;   %����
v = 150;        %�״��ٶ�150m/s
B = 150*10^6;   %�����źŴ���150MHz
Daz = 2;        %��λ���߳ߴ�2m
R0 = 20*10^3;   %�ο���б��20km
Tr = 1.5*10^-6;  %�������ʱ��
Kr = B/Tr;      %���Ե�Ƶ�źŵ�Ƶб��
H = 1000;       %�״�߶�
Yc = sqrt(R0^2-H^2);    %��������  
center = [0,Yc,0];      %������������
CR = Daz/2;             %����ֱ���
thetaaz = lamda/Daz;    %�������
Fs = 2.5*B;             %����Ƶ��
Dsar = lamda*R0/(2*CR);  %�ϳɿ׾�����
Tsar = Dsar/v;          %һ���ϳɿ׾�ʱ��
Ka = -2*v^2/(lamda*R0);  %�����յ�Ƶб��
Ba = abs(Ka*Tsar);      %�����մ���
PRF = 1.2*Ba;           %�����ظ�Ƶ��
X = 150;Y = 150;        %������С
Rmin = sqrt((Yc-Y/2)^2+H^2);  %������С����
Rmax = sqrt((Yc+Y/2)^2+H^2+(Dsar/2)^2);  %����������
Nfast = ceil(((2*(Rmax-Rmin)/c+Tr)*Fs));  %��ʱ��ά��������
tf = 2*R0/c+(-(Nfast/2):(Nfast/2)-1)/Fs;  %��ʱ���������
Nslow = ceil((X+Dsar)/v*PRF);  %��ʱ��ά����λ�򣩲�������
ts = (-(Nslow/2):(Nslow/2)-1)/PRF;  %��ʱ���������
pos = [-50,50,0,3];  %Ŀ��������ĵ�λ��,[x,y,z,rcs],������Ϊ������ϵ��
disp('Ŀ��λ�ã���λ��б�࣬�߶ȣ���');
Rpos(1:3) = pos(1:3)+center  %Ŀ�����λ��
Rpos(:,4) = pos(:,4);

%% �ز��ź�����
signal = zeros(Nfast,Nslow);
Xs = ts.*v-Rpos(1);  
Ys = 0-Rpos(2);
Zs = H-Rpos(3);
sigma = Rpos(4);  %����ϵ��
R = sqrt(Xs.^2+Ys^2+Zs^2);  %б��
tau = 2*R/c;  %ʱ��
Tfast = tf'*ones(1,Nslow)-ones(Nfast,1)*tau;   %���ǿ�ʱ��  ������ʱ��
Phase = pi*Kr.*Tfast.^2-(4*pi/lamda)*ones(Nfast,1)*R;   %��λ�ӳ�
signal = signal+sigma*exp(1i*Phase).*(abs(Tfast)<=Tr/2).*(ones(Nfast,1)*(abs(Xs)<=Dsar/2));  %�ز�
S = fftshift(fft(fftshift(signal)));

%% ��ά��ѹ
hf = exp(1i*pi*Kr*(tf-2*R0/c).^2).*(abs(tf-2*R0/c)<=Tr/2);  %������ο�����
Hf = (fftshift(fft(fftshift(hf))).')*ones(1,Nslow);
ComF = S.*conj(Hf);  %������ƥ���˲�
Sr = fftshift(ifft(fftshift(ComF)));  %������IFFT
Coms = fftshift(fft(fftshift(Sr.'))).';  %��λ��FFT
hs = exp(1i*pi*Ka*ts.^2).*(abs(ts)<Tsar/2);  %��λ��ο�����
Hs = ones(Nfast,1)*fftshift(fft(fftshift(hs)));
ComS = Coms.*conj(Hs);  %��λ��ƥ���˲�
Saz = fftshift(ifft(fftshift(ComS.'))).';

%% RD�㷨,sinc��ֵ
Coms_rcmc = zeros(Nfast,Nslow);
N = 10;  %��ֵ����
Rp = sqrt(sum((Rpos(2:3)-[0,H]).^2));  %Ŀ�굽�״���������
h = waitbar(0,'��ֵ��');  %����һ��������
for m = 1:Nslow  %��ʱ��
    for n = N/2+1:Nfast  %��ʱ��
      %����ƫ����
      deltaR = (lamda/v)^2*(Rp+(n-Nfast/2)*c/2/Fs)*((m-Nslow/2)/Nslow*PRF)^2/8;
      DU = deltaR/(c/2/Fs);  %ƫ�ƾ��뵥Ԫ
      deltaDU = DU-floor(DU);  %ƫ�ƾ��뵥ԪС������
      for k = -N/2:N/2-1
          if (n+floor(DU)+k)>Nfast %�����߽�
              Coms_rcmc(n,m) = Coms_rcmc(n,m)+Coms(Nfast,m)*sinc(DU-k);
          else
              Coms_rcmc(n,m) = Coms_rcmc(n,m)+Coms(n+floor(DU)+k,m)*sinc(deltaDU-k);
          end
      end
  end
  waitbar(m/Nslow);
end
close(h);  %�رս�����
ComS_rcmc = Coms_rcmc.*conj(Hs);  %��λ��ѹ��
Saz_rcmc = fftshift(ifft(fftshift(ComS_rcmc.'))).';

%% ��ͼ
rf = c*tf/2;  %����
az = v*ts;  %��λ
faz = (-Nslow/2:Nslow/2-1)/Nslow*PRF;  %������Ƶ��
figure(1);
[f,Rf] = meshgrid(faz,rf);
subplot(121);
mesh(f,Rf,abs(Coms));view(0,90);
title('(a) δRCMC');
xlabel('������/Hz');ylabel('б��R/m');
subplot(122);
mesh(f,Rf,abs(Coms_rcmc));view(0,90);
title('(b) RCMC');
xlabel('������/Hz');ylabel('б��R/m');

figure(2);
[Az,Rf] = meshgrid(az,rf);
mesh(Az,Rf,abs(Saz));
% view(0,90);
title('��ά��ѹ������');
xlabel('��λx/m');ylabel('б��R/m');

figure(3);
mesh(Az,Rf,abs(Saz_rcmc));
% view(0,90);
title('RD�㷨������');
xlabel('��λx/m');ylabel('б��R/m');

figure(4);
contour(abs(Saz_rcmc))%�ȸ���ͼ�λ���