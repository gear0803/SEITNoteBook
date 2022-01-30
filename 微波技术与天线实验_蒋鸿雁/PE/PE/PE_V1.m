clear;
close all;

%% 初始化，
cLight = 3e8;
frq = 1e+09;
K0 = 2*pi*frq/cLight;
lamda = cLight/frq;
dx = 0.5*lamda;
dz = 0.5*lamda;
X= 1000; %水平距离
Z= 100;  %高度
Xmax = floor(X/dx + 1);
Zmax = floor(Z/dz + 1);
dxt = 0:dx:X;
dzt = 0:dz:Z;

%% 设置天线
h_ant = 15;                   %天线高度
H_ant = floor(h_ant/dz + 1); 
pol = 0;                      %pol=0代表天线水平极化 1为垂直极化
beta = deg2rad(3);            %束宽
ceta = deg2rad(0);            %天线倾角
U = zeros(Xmax,Zmax);         %电场值
for j = 1:Zmax
    U(1,j) = (K0*beta/(2*((2*pi*log(2))^0.5)))*exp(-1i*K0*(j-1)*dz*ceta)...
        *exp(-(K0*beta*((j-1)*dz-(H_ant-1)*dz))^2/(8*log(2)));
end
U(1,:)=U(1,:)/max(U(1,:));    %归一化

%% 设置吸收窗函数
f1 = floor(Zmax/3);
filt = zeros(1,f1);
for n = 1:f1
    filt(1,n) = 0.5+0.5*cos(pi*n/f1);
end

%% 设置阻抗边界条件
%中等干燥地面
eps_g = 30;
sigma_g = 0.01;
eps_e = eps_g +1i*60*lamda*sigma_g;
if(pol==0)
    eta = (eps_e-1)^0.5;
else 
    eta = (eps_e-1)^0.5/eps_e;
end

%% 设置三对角矩阵元素
A = -2+4*1i*K0*dz^2/dx;
B = 2+4*1i*K0*dz^2/dx;
main = A*ones(1,Zmax);   %主对角线
sup = ones(1,Zmax-1);    %上对角
sub = ones(1,Zmax-1);    %下对角
rhs = zeros(1,Zmax);     %方程右边

%% 主计算
for n = 2:Xmax
    rhs(1) = B*U(n-1,1)-U(n-1,2);
    for m = 2:Zmax -1 
        rhs(m) = -U(n-1,m-1)+B*U(n-1,m)-U(n-1,m+1);
    end
    rhs(Zmax) = -U(n-1,Zmax-1)+B*U(n-1,Zmax);
    
    U(n,:) = Solve_Tridiagonal(main,sup,sub,rhs);
    U(n,Zmax-f1+1:Zmax)= U(n,Zmax-f1+1:Zmax).*filt;
    
end

% save('data.mat','-v7.3')

figure;
imagesc(dxt,dzt,abs(U'))
axis xy;
set(gca,'Fontsize',16);
xlim([0,X]);
ylim([0,Z]);
colormap(jet);
colorbar('EastOutside','Fontsize',16)
title('FIeld intensity')
xlabel('Propagatin Distance /km')
ylabel('Propagation Height /m');
saveas(gcf,'U','tif');

%% 解释：
%{
1、Ex的最大值在左边Z=15处，是因为设置了天线高度就为15
2、天线从Z=15,X=0处发射信号，向右边传播，在水平传播中有损耗，所以幅值不断减小
3、随着高度的增大，电场的幅值由于扩散而逐渐减小。
4、在X=800,Z=70处的电场幅值消失，这是由于顶部使用了吸收窗，让电磁波在吸收层中缓慢减至0
5、电场有干涉波束现象，是不同电磁波在空间重叠产生的结果。
5、X=300以后开始出现波前倾斜现象，这是由于地面损耗造成电场向传播方向倾斜的一种现象

%}




