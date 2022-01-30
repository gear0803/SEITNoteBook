%***********************************************************************
%     3-D FDTD code with PEC boundaries
%***********************************************************************
%
%     Program author: Susan C. Hagness
%                     Department of Electrical and Computer Engineering
%                     University of Wisconsin-Madison
%                     1415 Engineering Drive
%                     Madison, WI 53706-1691
%                     608-265-5739
%                     hagness@engr.wisc.edu
%
%     Date of this version:  February 2000
%
%     This MATLAB M-file implements the finite-difference time-domain
%     solution of Maxwell's curl equations over a three-dimensional
%     Cartesian space lattice comprised of uniform cubic grid cells.
%     
%     To illustrate the algorithm, an air-filled rectangular cavity 
%     resonator is modeled.  The length, width, and height of the 
%     cavity are 10.0 cm (x-direction), 4.8 cm (y-direction), and 
%     2.0 cm (z-direction), respectively.
%
%     The computational domain is truncated using PEC boundary 
%     conditions:
%          ex(i,j,k)=0 on the j=1, j=jb, k=1, and k=kb planes
%          ey(i,j,k)=0 on the i=1, i=ib, k=1, and k=kb planes
%          ez(i,j,k)=0 on the i=1, i=ib, j=1, and j=jb planes
%     These PEC boundaries form the outer lossless walls of the cavity.
%
%     The cavity is excited by an additive current source oriented
%     along the z-direction.  The source waveform is a differentiated 
%     Gaussian pulse given by 
%          J(t)=-J0*(t-t0)*exp(-(t-t0)^2/tau^2), 
%     where tau=50 ps.  The FWHM spectral bandwidth of this zero-dc-
%     content pulse is approximately 7 GHz. The grid resolution 
%     (dx = 2 mm) was chosen to provide at least 10 samples per 
%     wavelength up through 15 GHz.
%
%     To execute this M-file, type "fdtd3D" at the MATLAB prompt.
%     This M-file displays the FDTD-computed Ez fields at every other
%     time step, and records those frames in a movie matrix, M, which 
%     is played at the end of the simulation using the "movie" command.
%
%***********************************************************************

clear all
clc
close all

%***********************************************************************
%     Fundamental constants
%***********************************************************************

cc=2.99792458e8;                 %真空光速
muz=4.0*pi*1.0e-7;               %真空磁导率
epsz=1.0/(cc*cc*muz);            %真空介电常数
frequency=3e6;                   %电波频率

lambda = 3e8/frequency;          %波长 100
%***********************************************************************
%     Grid parameters
%***********************************************************************

ie=100;       %总计算范围（设置为正方体） 一半：50个dx对应5个波长
je=100;
ke=100;

ib=ie+1;     
jb=je+1;   
kb=ke+1;   

is=ie/2;       %发射源位置
js=je/2;       
ks=ke/2;

kobs=5;

dx=10;          %空间差分步长，单位m   100/10 = 10
dt=dx/(2.0*cc);    %时间差分步长，此处为走过半个空间差分步长的时间（走过5m的时间），可修改，但不得大于CFL条件
                    %一个dt为走过1/20波长的时间，所以20dt为走过一个波长的时间
nmax=100;          %总时间迭代次数5T = 5*20dt = 100dt


%***********************************************************************
%     Material parameters
%***********************************************************************

eps=1.0;
sig=0.0;        

%***********************************************************************
%     Updating coefficients
%***********************************************************************

ca=(1.0-(dt*sig)/(2.0*epsz*eps))/(1.0+(dt*sig)/(2.0*epsz*eps));
cb=(dt/epsz/eps/dx)/(1.0+(dt*sig)/(2.0*epsz*eps));
da=1.0;
db=dt/muz/dx;

%***********************************************************************
%     Field arrays
%***********************************************************************

ex=zeros(ie,jb,kb);
ey=zeros(ib,je,kb);
ez=zeros(ib,jb,ke);
hx=zeros(ib,je,ke);
hy=zeros(ie,jb,ke);
hz=zeros(ie,je,kb);

%***********************************************************************
%     BEGIN TIME-STEPPING LOOP
%***********************************************************************
% 申明source  wavnum1  wavnum2
source = zeros(1,nmax);
wavnum1 = zeros(1,nmax);
wavnum2 = zeros(1,nmax);

for n=1:nmax
   
%***********************************************************************
%     Update electric fields
%***********************************************************************

ex(1:ie,2:je,2:ke)=ca*ex(1:ie,2:je,2:ke)+...
                   cb*(hz(1:ie,2:je,2:ke)-hz(1:ie,1:je-1,2:ke)+...
                       hy(1:ie,2:je,1:ke-1)-hy(1:ie,2:je,2:ke));

ey(2:ie,1:je,2:ke)=ca*ey(2:ie,1:je,2:ke)+...
                   cb*(hx(2:ie,1:je,2:ke)-hx(2:ie,1:je,1:ke-1)+...
                       hz(1:ie-1,1:je,2:ke)-hz(2:ie,1:je,2:ke));
                    
ez(2:ie,2:je,1:ke)=ca*ez(2:ie,2:je,1:ke)+...
                   cb*(hx(2:ie,1:je-1,1:ke)-hx(2:ie,2:je,1:ke)+...
                       hy(2:ie,2:je,1:ke)-hy(1:ie-1,2:je,1:ke));
                    


ex(is,js,ks)=sin(2*pi*frequency*dt*n);
source(n)=ex(is,js,ks);                %发射源（z=0处），x极化波 
wavnum1(n)=ex(is,js,ks+10);     %沿z方向传播，距发射源一个波长位置的观察点（10*dx=100m=波长）
wavnum2(n)=ex(is,js,ks+20);     %沿z方向传播，距发射源一个波长位置的观察点（20*dx=200m=2波长）


%***********************************************************************
%     Update magnetic fields
%***********************************************************************

hx(2:ie,1:je,1:ke)=hx(2:ie,1:je,1:ke)+...
                   db*(ey(2:ie,1:je,2:kb)-ey(2:ie,1:je,1:ke)+...
                       ez(2:ie,1:je,1:ke)-ez(2:ie,2:jb,1:ke));
                
hy(1:ie,2:je,1:ke)=hy(1:ie,2:je,1:ke)+...
                   db*(ex(1:ie,2:je,1:ke)-ex(1:ie,2:je,2:kb)+...
                       ez(2:ib,2:je,1:ke)-ez(1:ie,2:je,1:ke));
                
hz(1:ie,1:je,2:ke)=hz(1:ie,1:je,2:ke)+...
                   db*(ex(1:ie,2:jb,2:ke)-ex(1:ie,1:je,2:ke)+...
                       ey(1:ie,1:je,2:ke)-ey(2:ib,1:je,2:ke));
                    

%***********************************************************************
%     END TIME-STEPPING LOOP
%***********************************************************************

end

figure;
plot(source)
xlabel('Time step');
ylabel('Amplitude of Ex');
title('dt=d/(2*c)  z = 0')

figure;
plot(wavnum1)
xlabel('Time step');
ylabel('Amplitude of Ex');
title('dt=d/(2*c) z = \lambda')

figure;
plot(wavnum2)
xlabel('Time step');
ylabel('Amplitude of Ex');
title('dt=d/(2*c) z = 2\lambda')


%% 解释 时间分辨率为：dt=d/(2*c)  
%{
Z=0：
波谷的周期为20个timestep，因为一个dt为走过1/20波长的时间，
所以20dt为走过一个波长的时间
波的幅值为1，因为实在发射源处

Z=\lambda 一个波长:
1、在波长为1时，理论要20个timestep才能有信号，因为要经过一个波长的时间才能传到
2、但其实11步就有数据（看wavnum1变量），是因为时间步长为1到9时源信号还没有传过去 
    而11开始有很小幅值的信号，该超光速存在是因为计算机在数值计算中存在数值离散，
    中心差分导致的，而可以忽略不计是因为算法本身设置就可以大致满足物理规律，
    因此可以当做光速是可以传到20才开始有信号
3、在20~40个timestep时波形不稳定，这是因为发射源不是由电流源激励，
    所以波形会有大的跳变，跳完就稳定了 所以后面周期差为20timestep
4、幅值相比于Z=0处变小，是因为信号在空间传输的损耗引起的，满足物理规律

Z = 2*\lambda 两个波长：
分析与Z=1个波长类似：
1、在波长为2时，理论要40个timestep才能有信号，因为要经过2个波长的时间才能传到
2、但其实21步就有数据，是因为时间步长为1到20时源信号还没有传过去 
    而21timestep时存在很小幅值的信号，该超光速存在是因为计算机在数值计算中存在数值离散，
    中心差分导致的，而可以忽略不计是因为算法本身设置就可以大致满足物理规律，
    因此可以当做光速是可以传到40才开始有信号
3、在40~60个timestep时波形不稳定，这是因为发射源不是由电流源激励，
    所以波形会有大的跳变，跳完就稳定了 所以后面周期差为20timestep
4、幅值相比于Z=1处变小，是因为信号在空间传输的损耗引起的，满足物理规律
%}


