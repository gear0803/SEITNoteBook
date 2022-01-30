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

cc=2.99792458e8;                 %��չ���
muz=4.0*pi*1.0e-7;               %��մŵ���
epsz=1.0/(cc*cc*muz);            %��ս�糣��
frequency=3e6;                   %�粨Ƶ��

lambda = 3e8/frequency;          %���� 100
%***********************************************************************
%     Grid parameters
%***********************************************************************

ie=100;       %�ܼ��㷶Χ������Ϊ�����壩 һ�룺50��dx��Ӧ5������
je=100;
ke=100;

ib=ie+1;     
jb=je+1;   
kb=ke+1;   

is=ie/2;       %����Դλ��
js=je/2;       
ks=ke/2;

kobs=5;

dx=10;          %�ռ��ֲ�������λm   100/10 = 10
dt=dx/(2.0*cc);    %ʱ���ֲ������˴�Ϊ�߹�����ռ��ֲ�����ʱ�䣨�߹�5m��ʱ�䣩�����޸ģ������ô���CFL����
                    %һ��dtΪ�߹�1/20������ʱ�䣬����20dtΪ�߹�һ��������ʱ��
nmax=100;          %��ʱ���������5T = 5*20dt = 100dt


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
% ����source  wavnum1  wavnum2
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
source(n)=ex(is,js,ks);                %����Դ��z=0������x������ 
wavnum1(n)=ex(is,js,ks+10);     %��z���򴫲����෢��Դһ������λ�õĹ۲�㣨10*dx=100m=������
wavnum2(n)=ex(is,js,ks+20);     %��z���򴫲����෢��Դһ������λ�õĹ۲�㣨20*dx=200m=2������


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


%% ���� ʱ��ֱ���Ϊ��dt=d/(2*c)  
%{
Z=0��
���ȵ�����Ϊ20��timestep����Ϊһ��dtΪ�߹�1/20������ʱ�䣬
����20dtΪ�߹�һ��������ʱ��
���ķ�ֵΪ1����Ϊʵ�ڷ���Դ��

Z=\lambda һ������:
1���ڲ���Ϊ1ʱ������Ҫ20��timestep�������źţ���ΪҪ����һ��������ʱ����ܴ���
2������ʵ11���������ݣ���wavnum1������������Ϊʱ�䲽��Ϊ1��9ʱԴ�źŻ�û�д���ȥ 
    ��11��ʼ�к�С��ֵ���źţ��ó����ٴ�������Ϊ���������ֵ�����д�����ֵ��ɢ��
    ���Ĳ�ֵ��µģ������Ժ��Բ�������Ϊ�㷨�������þͿ��Դ�������������ɣ�
    ��˿��Ե��������ǿ��Դ���20�ſ�ʼ���ź�
3����20~40��timestepʱ���β��ȶ���������Ϊ����Դ�����ɵ���Դ������
    ���Բ��λ��д�����䣬������ȶ��� ���Ժ������ڲ�Ϊ20timestep
4����ֵ�����Z=0����С������Ϊ�ź��ڿռ䴫����������ģ������������

Z = 2*\lambda ����������
������Z=1���������ƣ�
1���ڲ���Ϊ2ʱ������Ҫ40��timestep�������źţ���ΪҪ����2��������ʱ����ܴ���
2������ʵ21���������ݣ�����Ϊʱ�䲽��Ϊ1��20ʱԴ�źŻ�û�д���ȥ 
    ��21timestepʱ���ں�С��ֵ���źţ��ó����ٴ�������Ϊ���������ֵ�����д�����ֵ��ɢ��
    ���Ĳ�ֵ��µģ������Ժ��Բ�������Ϊ�㷨�������þͿ��Դ�������������ɣ�
    ��˿��Ե��������ǿ��Դ���40�ſ�ʼ���ź�
3����40~60��timestepʱ���β��ȶ���������Ϊ����Դ�����ɵ���Դ������
    ���Բ��λ��д�����䣬������ȶ��� ���Ժ������ڲ�Ϊ20timestep
4����ֵ�����Z=1����С������Ϊ�ź��ڿռ䴫����������ģ������������
%}


