%2D FFT beamforming
%clc;
clear all;
%close all;
c=3e8;
%frequency and wavelength of observing source
f=1420e6;
lambda=c/f;
ls = 2;
%distance between antennas in y direction
dx = lambda/ls;
%distance between antennas in x direction
dy = lambda/ls;
%propagation constant
k = 2*pi/lambda;
%Number of antennas in y direction
M = 6;
%Number of antennas in x direction
N = 6;
%Progressive Phase shift in y and x direction 
Bx = 0;
By = 0;
%Find maximun of M and N and create a square matrix of max of M and N if M and N are not equal
if (M == N)
  A =  ones(M,N);
else 
  r = max(M,N);
B = ones(r,r);
A = ones(M,N);
%Amplitudes of each antenna in 8X8 array: 64 elements
A(size(B,1), size(B,2)) = 0;
endif


%Im(:,1) = [1+j 0.5+0.5j 0.6+0.6j 0.8+0.8j 1+j 0.5+0.5j 0.6+0.6j 0.8+0.8j];
%In(1,:) = [1+j 0.5+0.5j 0.6+0.6j 0.8+0.8j 1+j 0.5+0.5j 0.6+0.6j 0.8+0.8j];
%Zenith angle
thetad900 = -90:1:0;
thetar900 = deg2rad(thetad900);
%Azimuth angle
phi = 45;
phir = deg2rad(phi);
%Initial Array Factor 
AFx = 0;
AFy = 0;
w =  zeros(1,M);
z = zeros(1,N);
%AF for all the columns
m = 1:1:M;
for y = 1:N
 for p =  1:length(thetar900)
 wz9 = A(y,m).*exp(1i.*(m-1)*(k.*dx.*sin(thetar900(p)).*cos(phir))); 
 AFxz9(p) = sum(wz9);
endfor
endfor

%AF for all the rows
n = 1:1:M;
for x = 1:M
 for p =  1:length(thetar900)
 zz9 = A(x,n).*exp(1i.*(n-1)*(k.*dy.*sin(thetar900(p)).*sin(phir)));
 AFyz9(p) = sum(zz9);
endfor
endfor

%multiply array factor of rows and columns
%AFxl = 20*log10(AFx);
subplot(3,3,1)
plot(thetad900,AFxz9);
ylabel('Array factor of x');
xlabel('zenith-90-0,azi-0');
title('AFx')
%AFyl = 20*log10(AFy);
subplot(3,3,2)
plot(thetad900,AFyz9);
ylabel('Array factor of y');
xlabel('zenith-90-0,azi-0');
title('AFy');

  AFz9 = AFxz9.*AFyz9;
  %y = 20*log10(AF);
  subplot(3,3,3)
plot(thetad900,AFz9);
ylabel('Array factor of x and y');
xlabel('Czenith-90-0,azi-0');
title('Product of AFx and AFy')



%Zenith angle
thetad090 = 0:1:90;
thetar090 = deg2rad(thetad090);
%Azimuth angle
phi = 225;
phir = deg2rad(phi);
%Initial Array Factor 
AFx = 0;
AFy = 0;
w =  zeros(1,M);
z = zeros(1,N);
%AF for all the columns
m = 1:1:M;
for y = 1:N
 for p =  1:length(thetar090)
 wz0 = A(y,m).*exp(1i.*(m-1)*(k.*dx.*sin(thetar090(p)).*cos(phir))); 
 AFxz0(p) = sum(wz0);
endfor
endfor

%AF for all the rows
n = 1:1:M;
for x = 1:M
 for p =  1:length(thetar090)
 zz0 = A(x,n).*exp(1i.*(n-1)*(k.*dy.*sin(thetar090(p)).*sin(phir)));
 AFyz0(p) = sum(zz0);
endfor
endfor

%multiply array factor of rows and columns
%AFxl = 20*log10(AFx);
subplot(3,3,4)
plot(thetad090,AFxz0);
ylabel('Array factor of x');
xlabel('zenith-0-90,azi-180');
title('AFx')
%AFyl = 20*log10(AFy);
subplot(3,3,5)
plot(thetad090,AFyz0);
ylabel('Array factor of y');
xlabel('zenith-0-90,azi-180');
title('AFy');

  AFz0 = AFxz0.*AFyz0;
  %y = 20*log10(AF);
  subplot(3,3,6)
plot(thetad090,AFz0);
ylabel('Array factor of x and y');
xlabel('zenith-0-90,azi-180');
title('Product of AFx and AFy')

AF = [AFz9 AFz0];
theta = [thetad900 thetad090];
subplot(3,3,8)
plot(theta,AF);
legend(num2str(1/ls));
ylabel('Array factor');
xlabel('zenith--90 to 90,azi-0 and 180');
title('beam response')
