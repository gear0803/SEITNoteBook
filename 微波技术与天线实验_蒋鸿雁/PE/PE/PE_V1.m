clear;
close all;

%% ��ʼ����
cLight = 3e8;
frq = 1e+09;
K0 = 2*pi*frq/cLight;
lamda = cLight/frq;
dx = 0.5*lamda;
dz = 0.5*lamda;
X= 1000; %ˮƽ����
Z= 100;  %�߶�
Xmax = floor(X/dx + 1);
Zmax = floor(Z/dz + 1);
dxt = 0:dx:X;
dzt = 0:dz:Z;

%% ��������
h_ant = 15;                   %���߸߶�
H_ant = floor(h_ant/dz + 1); 
pol = 0;                      %pol=0��������ˮƽ���� 1Ϊ��ֱ����
beta = deg2rad(3);            %����
ceta = deg2rad(0);            %�������
U = zeros(Xmax,Zmax);         %�糡ֵ
for j = 1:Zmax
    U(1,j) = (K0*beta/(2*((2*pi*log(2))^0.5)))*exp(-1i*K0*(j-1)*dz*ceta)...
        *exp(-(K0*beta*((j-1)*dz-(H_ant-1)*dz))^2/(8*log(2)));
end
U(1,:)=U(1,:)/max(U(1,:));    %��һ��

%% �������մ�����
f1 = floor(Zmax/3);
filt = zeros(1,f1);
for n = 1:f1
    filt(1,n) = 0.5+0.5*cos(pi*n/f1);
end

%% �����迹�߽�����
%�еȸ������
eps_g = 30;
sigma_g = 0.01;
eps_e = eps_g +1i*60*lamda*sigma_g;
if(pol==0)
    eta = (eps_e-1)^0.5;
else 
    eta = (eps_e-1)^0.5/eps_e;
end

%% �������ԽǾ���Ԫ��
A = -2+4*1i*K0*dz^2/dx;
B = 2+4*1i*K0*dz^2/dx;
main = A*ones(1,Zmax);   %���Խ���
sup = ones(1,Zmax-1);    %�϶Խ�
sub = ones(1,Zmax-1);    %�¶Խ�
rhs = zeros(1,Zmax);     %�����ұ�

%% ������
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

%% ���ͣ�
%{
1��Ex�����ֵ�����Z=15��������Ϊ���������߸߶Ⱦ�Ϊ15
2�����ߴ�Z=15,X=0�������źţ����ұߴ�������ˮƽ����������ģ����Է�ֵ���ϼ�С
3�����Ÿ߶ȵ����󣬵糡�ķ�ֵ������ɢ���𽥼�С��
4����X=800,Z=70���ĵ糡��ֵ��ʧ���������ڶ���ʹ�������մ����õ�Ų������ղ��л�������0
5���糡�и��沨�������ǲ�ͬ��Ų��ڿռ��ص������Ľ����
5��X=300�Ժ�ʼ���ֲ�ǰ��б�����������ڵ��������ɵ糡�򴫲�������б��һ������

%}




