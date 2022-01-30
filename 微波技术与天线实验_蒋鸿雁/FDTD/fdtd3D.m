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

cc=2.99792458e8;            %真空光速
muz=4.0*pi*1.0e-7;          %真空磁导率
epsz=1.0/(cc*cc*muz);       %真空介电常数
frequency=3e8;                %电波频率

%***********************************************************************
%     Grid parameters
%***********************************************************************

ie=80;       %总计算范围（设置为正方体）
je=80;
ke=80;

ib=ie+1;     
jb=je+1;   
kb=ke+1;   

is=ie/2;       %发射源位置
js=je/2;       
ks=ke/2;

kobs=5;

dx=0.1;          %空间差分步长，单位m
dt=dx/(2.0*cc);    %时间差分步长，此处为走过半个空间差分步长的时间（走过0.05m的时间），可修改，但不得大于CFL条件

nmax=120;          %总时间迭代次数  0~5T  2*10dt 为1个周期


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
wavnum1(n)=ex(is,js,ks+10);     %沿z方向传播，距发射源一个波长位置的观察点（10*dx=1m=波长）
wavnum2(n)=ex(is,js,ks+20);


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

plot(source)
xlabel('Time step');
ylabel('Amplitude of Ex');

figure
plot(wavnum1)
xlabel('Time step');
ylabel('Amplitude of Ex');

figure
plot(wavnum2)
xlabel('Time step');
ylabel('Amplitude of Ex');

% movie(gcf,M,0,10,rect);
