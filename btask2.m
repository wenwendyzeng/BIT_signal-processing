%仿真具有3个以上距离徙动回波信号，进行RD成像
clc;
clear all;

%% 基本参数
c = 3*10^8;     %波速
fc = 5.3*10^9;  %载频5.3GHz
lamda = c/fc;   %波长
v = 150;        %雷达速度150m/s
B = 150*10^6;   %发射信号带宽150MHz
Daz = 2;        %方位天线尺寸2m
R0 = 20*10^3;   %参考点斜距20km
Tr = 1.5*10^-6;  %脉冲持续时间
Kr = B/Tr;      %线性调频信号调频斜率
H = 1000;       %雷达高度
Yc = sqrt(R0^2-H^2);    %成像中线  
center = [0,Yc,0];      %场景中心坐标
CR = Daz/2;             %横向分辨率
thetaaz = lamda/Daz;    %波束宽度
Fs = 2.5*B;             %采样频率
Dsar = lamda*R0/(2*CR);  %合成孔径长度
Tsar = Dsar/v;          %一个合成孔径时间
Ka = -2*v^2/(lamda*R0);  %多普勒调频斜率
Ba = abs(Ka*Tsar);      %多普勒带宽
PRF = 1.2*Ba;           %脉冲重复频率
X = 150;Y = 150;        %场景大小
Rmin = sqrt((Yc-Y/2)^2+H^2);  %场景最小距离
Rmax = sqrt((Yc+Y/2)^2+H^2+(Dsar/2)^2);  %场景最大距离
Nfast = ceil(((2*(Rmax-Rmin)/c+Tr)*Fs));  %快时间维采样点数
tf = 2*R0/c+(-(Nfast/2):(Nfast/2)-1)/Fs;  %快时间采样序列
Nslow = ceil((X+Dsar)/v*PRF);  %慢时间维（方位向）采样点数
ts = (-(Nslow/2):(Nslow/2)-1)/PRF;  %慢时间采样序列
pos = [-50,50,0,3];  %目标相对中心点位置,[x,y,z,rcs],第四列为后向反射系数
disp('目标位置（方位，斜距，高度）：');
Rpos(1:3) = pos(1:3)+center  %目标绝对位置
Rpos(:,4) = pos(:,4);

%% 回波信号生成
signal = zeros(Nfast,Nslow);
Xs = ts.*v-Rpos(1);  
Ys = 0-Rpos(2);
Zs = H-Rpos(3);
sigma = Rpos(4);  %反射系数
R = sqrt(Xs.^2+Ys^2+Zs^2);  %斜距
tau = 2*R/c;  %时延
Tfast = tf'*ones(1,Nslow)-ones(Nfast,1)*tau;   %列是快时间  行是慢时间
Phase = pi*Kr.*Tfast.^2-(4*pi/lamda)*ones(Nfast,1)*R;   %相位延迟
signal = signal+sigma*exp(1i*Phase).*(abs(Tfast)<=Tr/2).*(ones(Nfast,1)*(abs(Xs)<=Dsar/2));  %回波
S = fftshift(fft(fftshift(signal)));

%% 二维脉压
hf = exp(1i*pi*Kr*(tf-2*R0/c).^2).*(abs(tf-2*R0/c)<=Tr/2);  %距离向参考函数
Hf = (fftshift(fft(fftshift(hf))).')*ones(1,Nslow);
ComF = S.*conj(Hf);  %距离向匹配滤波
Sr = fftshift(ifft(fftshift(ComF)));  %距离向IFFT
Coms = fftshift(fft(fftshift(Sr.'))).';  %方位向FFT
hs = exp(1i*pi*Ka*ts.^2).*(abs(ts)<Tsar/2);  %方位向参考函数
Hs = ones(Nfast,1)*fftshift(fft(fftshift(hs)));
ComS = Coms.*conj(Hs);  %方位向匹配滤波
Saz = fftshift(ifft(fftshift(ComS.'))).';

%% RD算法,sinc插值
Coms_rcmc = zeros(Nfast,Nslow);
N = 10;  %插值点数
Rp = sqrt(sum((Rpos(2:3)-[0,H]).^2));  %目标到雷达的最近距离
h = waitbar(0,'插值中');  %生成一个进度条
for m = 1:Nslow  %慢时间
    for n = N/2+1:Nfast  %快时间
      %距离偏移量
      deltaR = (lamda/v)^2*(Rp+(n-Nfast/2)*c/2/Fs)*((m-Nslow/2)/Nslow*PRF)^2/8;
      DU = deltaR/(c/2/Fs);  %偏移距离单元
      deltaDU = DU-floor(DU);  %偏移距离单元小数部分
      for k = -N/2:N/2-1
          if (n+floor(DU)+k)>Nfast %超出边界
              Coms_rcmc(n,m) = Coms_rcmc(n,m)+Coms(Nfast,m)*sinc(DU-k);
          else
              Coms_rcmc(n,m) = Coms_rcmc(n,m)+Coms(n+floor(DU)+k,m)*sinc(deltaDU-k);
          end
      end
  end
  waitbar(m/Nslow);
end
close(h);  %关闭进度条
ComS_rcmc = Coms_rcmc.*conj(Hs);  %方位向压缩
Saz_rcmc = fftshift(ifft(fftshift(ComS_rcmc.'))).';

%% 画图
rf = c*tf/2;  %距离
az = v*ts;  %方位
faz = (-Nslow/2:Nslow/2-1)/Nslow*PRF;  %多普勒频率
figure(1);
[f,Rf] = meshgrid(faz,rf);
subplot(121);
mesh(f,Rf,abs(Coms));view(0,90);
title('(a) 未RCMC');
xlabel('多普勒/Hz');ylabel('斜距R/m');
subplot(122);
mesh(f,Rf,abs(Coms_rcmc));view(0,90);
title('(b) RCMC');
xlabel('多普勒/Hz');ylabel('斜距R/m');

figure(2);
[Az,Rf] = meshgrid(az,rf);
mesh(Az,Rf,abs(Saz));
% view(0,90);
title('二维脉压成像结果');
xlabel('方位x/m');ylabel('斜距R/m');

figure(3);
mesh(Az,Rf,abs(Saz_rcmc));
% view(0,90);
title('RD算法成像结果');
xlabel('方位x/m');ylabel('斜距R/m');

figure(4);
contour(abs(Saz_rcmc))%等高线图形绘制