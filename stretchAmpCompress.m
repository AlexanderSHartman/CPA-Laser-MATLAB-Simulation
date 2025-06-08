%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% SA01 Program x
% Program for pulse Streching, Amplification in CPA System x
% Stretching fiber : Lst, beta_st, alfa_st x
% Amplification fiber: Lam, beta_am, g_am, alfa_am x
% Using: RK4IP algorithm (Hunt JLT 2007) x
% fucntions: x
% Propagation_EDF.m for propagation in EDF x
% Propagation_SMF.m for propagation in SMF x
% x
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
close all
clc
if ~exist('SA01','dir')
mkdir('SA01');
end
%if ~exist('SA01b','dir')
% mkdir('SA01b');
%end
%if ~exist('SA01c','dir')
% mkdir('SA01c');
%end
%% ---- parameters for the simulation and a part of the pulse ---------
c0 = 3e8; % (m/s)
lambdaP = 1550e-9; % (m)% pump wavelength (3000 , 3500, 6500, 7250 nm)
wP = 2*pi*c0/lambdaP; % Angular frequency of pump @ 1560nm
%TFWHM = 500e-15; % (s) full width at half maximum of pump pulse

TFWHM = 175e-15; %pulse width given in paper 175fs

%Tmax = 600e-12; % (s), half width of time window
Tmax = 600e-12;
Tmin = -Tmax; % Tmin < t < Tmax
tol = 1e-3; % tollerence
global t % vector time, globally used (even in functions)
N = 2^14; % number of points used to present in time domain
deltat = 2*Tmax/(N-1); % (s), resolution of time
t = linspace(Tmin,Tmax-deltat,N);
global f % frequency, globally used
deltaf = 1/(2*Tmax); % resolution of frequency
f = [0:deltaf:(N/2-1)*deltaf,-N*deltaf/2:deltaf:-deltaf]; % vertor freq
w = 2*pi*f; % vector angular freq
wshift = -w + wP; % plus the freq with its central freq.
% the minus sign is to compensate the fact that
% defenition of FFt in matlab and in the GNLSE are
% in opposited sign
wl = 1e9*2*pi*c0./wshift;%(nm), shift to nm to present
tm = 1e12*t; % (ps), shift to ps to present
limit_t = [-50 50]; % (ps)
limit_f = [1520 1580]; % (nm)
%% -----------This part constructs the pump @1550nm ----------
C = 0; % chirp of the pulse
form_pulse = 'gaus' % 'gaus'
if form_pulse == 'sech'
T0 = TFWHM/(2*log(1+sqrt(2))); % if pulses are secant-hyperbol (soliton)
U0 = sech(t/T0); %
else
T0 = TFWHM/(2*sqrt(log(2))); % if pulses are gaussian (bellshape)
U0 = exp(-0.5*(t/T0).^2*(1+i*C));
end
% ------------- from mode-locked laser before EDF ------------------------
R = 204e6; % repetation rate 204 MHz
Tau = 1/R;

%Pin = 1*1e-3 % output of seed laser is about 1mW
Pavg = 140e-3 %average power is 140mW from paper

%P0 = (0.93*Pavg*Tau)/TFWHM % (W), peak power of pulses (formula used to
%calculate P0
P0 = 3.6471e3;


Uint = sqrt(P0)*U0; % (W^0.5), Input pulse E-field, in time domain
Uinf = fft(Uint); % Input spectrum, freq domain
Iint = abs(Uint).^2;
Iinf = abs(Uinf).^2;
check0 = sum(Iint); % Input power = integral of (Uint)^2
%---------------
t_shock = 1/wP; % (s), paramet of self steepening & shock formation
fR = 0.18; % factor of contribution of Raman effect
t1 = 12.2e-15; % (s), paramet in approx raman response function
t2 = 32e-15; % (s), paramet in approx raman response function,
% 1/t2 is the bandwidth of Raman gain
tr = t - t(1); % shift time vector to the root. starting point is 0
% Raman response function approx by Blow and Wood, 1989.
hR = (t1^2+t2^2)/t1/t2^2*exp(-tr/t2).*sin(tr/t1);
hR = hR./trapz(tr,hR); % (??? why divided, not affect to results ???)
hR_f = fft(hR); % raman response function in the freq domain
tt = 50e-15; % effect of gain in dispers value [Agrawal OL1991]
% --- This part for streching fiber ---------------------------------------
% Assume fiber is SMF-28
%-----------------------------------------------------------------------
Lst = 300; % Length of Strecher (m)
La1 = 10; % Length of Amplifier (m)
Lcm = 110; % Length of Compressor (m)
n2 = 2.6e-20; % nonlinear refractive index
MFD1 = 10e-6 ; % MFD of the fiber, assume 10um (check)
AdB_st = 7.5; % (dB/km) loss of strech fiber
alfa_st= AdB_st/4343; % (1/m) loss
Aeff1 = pi*(MFD1/2)^2; % effective core area
gamma_st = 2*pi*n2/lambdaP/Aeff1 % 0.0013/((2)^2);
% [beta1, beta2, beta3...] in ps^k/m
beta_st = [0 -6.2 2.6 -10.40].*[0 1e-2 1e-4 1e-7];
% change unit of beta to (s^k/m)
for ii = 1:length(beta_st)
beta_st(ii) = (1e-12)^(ii)*beta_st(ii);
end

%% -- this part is for gain fiber
AdB_am1 = 7.5; % dB/km loss of silica fiber
alpha_am1=AdB_am1/4343; % loss in unit 1/m
MFD2a = 10e-6 ; % MFD of the fiber, assume 10um (check)
Aeff2a = pi*(MFD2a/2)^2; % effective core area
n2a1 = 2.6e-20; % nonlinear refractive index
gamma_am = 2*pi*n2a1/lambdaP/Aeff2a % 0.0013/((2)^2);
% [beta1, beta2, beta3...] in ps^k/m
beta_am1 = [0 -6.2 2.6 -10.40].*[0 1e-2 1e-4 1e-7];
% change unit of beta to (s^k/m)
for ii = 1:length(beta_am1)
beta_am1(ii) = (1e-12)^(ii)*beta_am1(ii);
end

%% -- this part is for compressing fiber
%% -- PM 1550 fiber (used in the paper, datasheet from thorlabs)
%% -- The PM is a SMF fiber and can be simulated using Propagation_SMF
AdB_cm1 = 1; % dB/km loss of silica fiber (datasheet says <1 @ 1.55um)
alpha_cm1=AdB_cm1/4343; % loss in unit 1/m
MFD2cm = 10.1e-6 ; % MFD of the fiber
Aeff2cm = pi*(MFD2cm/2)^2; % effective core area
n2cm1 = 2.6e-20; % nonlinear refractive index (for silica fibers)
gamma_cm = 2*pi*n2cm1/lambdaP/Aeff2cm % 0.0013/((2)^2);
% [beta1, beta2, beta3...] in ps^k/m
% Using "virtual fiber" techniques as discussed in dan nguyn book chapter 5
%for the compressor
beta_cm1 = 5*[0 6.2 -5*2.6 10000*10.40].*[0 1e-2 1e-4 1e-7];
% change unit of beta to (s^k/m)
for ii = 1:length(beta_cm1)
beta_cm1(ii) = (1e-12)^(ii)*beta_cm1(ii);
end

%-------Progagation in strecher fiber -------------------
% signal 1920nm propagates throu SMF w/o gain
CL1 = 0.5;
TotalLoss = (AdB_st*Lst/1000 + AdB_am1*La1/1000) %+ AdB_cm1*Lcm/1000)
TotalGain = TotalLoss + CL1 + 17
GdB = TotalGain/La1 % gain dB/m

%This is the unsaturated gain coefficient g0 of the system
g1 = GdB/4.343 % linear gain in per m
Esat = 20e-9 ; % 20nJ saturation energy

Nst = 5;
dLst = Lst/Nst;
Nam = 3;
dLam = La1/Nam;
%Ncm = 5;
Ncm = 5;
dLcm = Lcm/Ncm;
Nz = Nst + Nam + Ncm;
Ltot = Lst + La1 + Lcm;
Lz = 0;
figure(1)
%plot3(tm,0*ones(size(tm)),Iint/max(Iint),'r','LineWidth',2);
%setup the 3D plot
plot3(tm,0*ones(size(tm)),zeros(size(tm)));
hold on;


%3 Stages Loop
for nz=1:Nz+1

% Stretcher Stage
if nz < Nst+1 
Lz(nz) = nz*dLst;
Uout1 = Propagation_SMF(Uint,tol,Lz(nz),alfa_st,beta_st,gamma_st,t_shock,fR,hR_f);%
Uouf1 = fft(Uout1);
Iout1 = zeros(16384, nz);
Iout1(:,nz) = abs(Uout1.').^2; % intensity of output
Inormt1(:,nz) = abs(Uout1./max(Uout1)).^2; % intensity of output
Iouf1 = zeros(16384, nz);
Iouf1(:,nz) = abs(Uouf1).^2; % intensity in freq domain
Inormf1(:,nz)= abs(Uouf1./max(Uouf1)).^2; % intensity in freq domain

%plot the stretch
plot3(tm,Lz(nz)*ones(size(tm)),Iout1(:,nz),'b','LineWidth',0.5);
% plot3(tm,Lz(nz)*ones(size(tm)),Inormt1(:,nz),'b','LineWidth',0.5);
hold on;
if nz==Nst
check1 = sum(abs(Uout1).^2);
P1 = check1/check0*Pavg
end

% Amplification Stage
elseif nz<=Nst+Nam+1 
L2(nz) = (nz-Nst-1)*dLam;
Lz(nz) = Lst + L2(nz);
Uin2 = sqrt(1)*Uout1;
Uout2 = Propagation_EDF(Uout1,tol,L2(nz),g1,Esat,alpha_am1,beta_am1,gamma_am,t_shock,fR,hR_f);%
Uouf2 = fft(Uout2);
Iout2(:,nz) = abs(Uout2).^2;
Inormt2(:,nz) = abs(Uout2./max(Uout2)).^2;
Iouf2(:,nz) = abs(Uouf2).^2;
Inormf2(:,nz) = abs(Uouf2./max(Uouf2)).^2;
plot3(tm,Lz(nz)*ones(size(tm)),Iout2(:,nz),'g','LineWidth',0.5);
% plot3(tm,Lz(nz)*ones(size(tm)),Inormt2(:,nz),'g','LineWidth',0.5);
if nz==Nst + Nam
check2 = sum(abs(Uout2).^2);
P2 = check2/check0*Pavg
end

%Compression Stage
else nz<=Nst+Nam+1 

L3(nz) = (nz-Nst-Ncm-1)*dLcm;
Lz(nz) = Lst + La1 + L3(nz);
Uin3 = sqrt(1)*Uout2;
%Uout3 = Propagation_EDF(Uout2,tol,L3(nz),0.01,Esat,alpha_cm1,beta_cm1,gamma_cm,t_shock,fR,hR_f);%
Uout3 = Propagation_SMF(Uout2,tol,L3(nz),alpha_cm1,beta_cm1,gamma_cm,t_shock,fR,hR_f);%
Uouf3 = fft(Uout3);
Iout3(:,nz) = abs(Uout3).^2;
Inormt3(:,nz) = abs(Uout3./max(Uout3)).^2;
Iouf3(:,nz) = abs(Uouf3).^2;
Inormf3(:,nz) = abs(Uouf3./max(Uouf3)).^2;

%Some of the data is noisy/looks really time-jittered, Ignore it
if (nz~=10 && nz~=12)
plot3(tm,Lz(nz)*ones(size(tm)),Iout3(:,nz),'m','LineWidth',0.5);
nz
end
%plot3(tm,Lz(nz)*ones(size(tm)),Iout3(:,nz),'m','LineWidth',0.5);
% plot3(tm,Lz(nz)*ones(size(tm)),Inormt2(:,nz),'g','LineWidth',0.5);
if nz==Nz
check3 = sum(abs(Uout3).^2);
P3 = check3/check0*Pavg
end
end


%-----This part presents the DATA -----------------------------------

axis([-50 50 0 Ltot 0 20])
xlabel('Time (ps)', 'FontSize',16);
ylabel('L(m)', 'FontSize',16);
zlabel('Power (W)', 'FontSize',16);


grid on;
cd SA01
saveas(gcf, strcat('SAC01R','.bmp'));
cd ../
drawnow
end





figure(2); % present time signal
%plot(tm,Iint/max(Iint),'LineWidth',2,'Color','r');
%semilogy(tm,Iint/max(Iint),'LineWidth',2,'Color','k');
hold on;
plot(tm,Iout1(:,Nst),'LineWidth',2,'Color','b');
plot(tm,Iout2(:,Nst+Nam),'LineWidth',2,'Color','g');
plot(tm, Iout3(:,Nst+Nam+Ncm), 'Linewidth', 2, 'Color', 'm');
%semilogy(tm,Iout1(:,Nst),'LineWidth',2,'Color','b');
%semilogy(tm,Iout2(:,Nst+Nam),'LineWidth',2,'Color','g');
hold off;
legend('Stretcher', 'Amplifier', 'Compressor');
xlabel ('Time (ps)','FontSize',14);
ylabel ('Power(W)','FontSize',14);
grid on;
set(gca,'FontSize',12);
axis([-200 200 0 10]);
%set(gca,'xtick',[-600 -500 -400 -300 -200 -100 0 100 200 300 400 500 600])
title('Peak Power','FontSize',14)
%cd SA01b
% saveas(gcf, strcat('SA01Rb','.bmp'));
%cd ../
figure(3); % present time signal
hold on;
plot(wl,Iinf/max(Iinf),'LineWidth',2,'Color','r');
%semilogy(wl,Iinf/max(Iinf),'LineWidth',2,'Color','k');
plot(wl,Iouf1(:,Nst),'LineWidth',2,'Color','b');
plot(wl,Iouf2(:,Nst+Nam),'LineWidth',2,'Color','g');
%semilogy(wl,Iouf1(:,Nst),'LineWidth',2,'Color','b');
%semilogy(wl,Iouf2(:,Nst+Nam),'LineWidth',2,'Color','g');
hold off;
legend('P_0', 'Stretcher', 'Amp #1');
xlabel ('Wavelength (nm)','FontSize',14);
ylabel ('Power(W)','FontSize',14);
grid on;
set(gca,'FontSize',12);
axis([1520 1580 0 Inf]);
set(gca,'xtick',[1520 1530 1540 1550 1560 1570 1580])
title('(1km SMF) & (10m Amp #1)','FontSize',14)
drawnow



function [Uout, numFFT] = Propagation_SMF(imp,tol,long,alpha,beta,gamma,ts,fR,hR)
% %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% This function solves the GNLE x
% Using RK4IP x
% Hult JLT 2007 and Sinkin JLT 2003. x
% x
%-------------------------------------------------------------------------x
% INPUT
%
% imp - impulses, starting field amplitude (vector)
% tolerance - local error limitation
% long - length of fiber
% alpha - linear loss coefficient, ie, Pz=P0*exp(-alpha*z)
% betap - dispersion polynomial coefs, [beta_1 ... beta_m], or beta(w)
% gamma - nonlinearity coefficient
% t_shock - self steepening term
% fR - fraction of raman effect
% hR - the response function of Raman effect in freq domain
% L_save - distances along the fiber where we save the data
%-------------------------------------------------------------------------
% OUTPUT
%
% disp - n series of temporal output signal in function of distances
% dispf - n series of spectral output signal in function of distances
% Uout - field at the output
% numFFT - number of performed FFTs
%------------------------------------------------------------------------
global f t
omega = 2*pi*f;
E_z = imp;
fprintf(1, '\nCaculation in process... ');
dz = long/10000; % the first step is chosen random
alpha0 = alpha;
beta0 = beta;
gamma0 = gamma;
Z_prop = 0;
ii = 0;
while Z_prop<long
    ii=ii+1;
    if Z_prop+dz>long
        dz = long-Z_prop;
    end
Z_prop = Z_prop + dz;
% change parameters as Z_propagation change
m = 0;
alpha = alpha0*(1 + m*Z_prop);
beta = beta0*(1 + m*Z_prop);
gamma = gamma0*(1 + m*Z_prop);
%---------
Uf = rk4ip(rk4ip(E_z,omega,dz/2,beta,alpha,gamma,ts,fR,hR),...
omega,dz/2,beta,alpha,gamma,ts,fR,hR);
Uc = rk4ip(E_z,omega,dz,beta,alpha,gamma,ts,fR,hR);
error = sqrt(sum(abs(Uf-Uc).^2))/sqrt(sum(abs(Uf).^2));
factor = tol/error;
% discard the solution and recaculate with a half of old stepsize
if error > 2*tol
Z_prop = Z_prop-dz;
else
E_z = 16/15*Uf-1/15*Uc; % improve the accuracy of the solution
end
dz=dz*factor^(1/5);
fprintf(1, '\b\b\b\b\b\b%5.2f%%', Z_prop* 100.0 /long);
end
Uout = E_z;
numFFT = 3*16*ii;
end


function Uout = Propagation_EDF(Uin,tol,L,g0,Es,alpha,beta,gam,ts,fR,hR)
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% This function solves the GNLE x
% Using RK4IP x
% Hult JLT 2007 and Sinkin JLT 2003. x
%-------------------------------------------------------------------------x
% INPUT
%
% Ei - impulses, starting field amplitude (vector)
% tol - local error limitation
% long - length of fiber
% alpha - linear loss coefficient, ie, Pz=P0*exp(-alpha*z)
% betap - dispersion polynomial coefs, [beta_1 ... beta_m], or beta(w)
% gam - nonlinearity coefficient
% ts - self steepening term (t_shock)
% fR - fraction of raman effect
% hR - the response function of Raman effect in freq domain
%-------------------------------------------------------------------------
% OUTPUT
% Uout - field at the output
%------------------------------------------------------------------------
global f t
delta_t = abs(t(2)-t(1));
omega = 2*pi*f;
Efield = Uin;
fprintf(1, '\nCaculation in process... ');
dz = L/10000; % trial for 1st step
alpha0 = alpha;
beta0 = beta;
gamma0 = gam;
Z_prop = 0;
ii = 0;
while Z_prop<L
ii=ii+1;
if Z_prop+dz>L
dz = L-Z_prop;
end
Z_prop = Z_prop + dz;
% change parameters as Z_propagation change
tt = 40e-15; % T2 in Agrawal 1990 PTL
% g = g0*exp(-sum(abs(E_z).^2)*delta_t/E_sat) (Agrawal 1991)
% or approximation as
g = g0./(1 + sum(abs(Efield).^2)*delta_t/Es);
alpha = alpha0 - g; % loss and gain in NLSE
beta(2) = beta0(2) + 1i*tt*tt*g; % from eq. (2) Agrawal PTL 1991
m = 0; %
gam = gamma0*(1 + m*Z_prop);
%---------optimization of local errors by Sinkin JLWT 2003 ------------
Uf = rk4ip(rk4ip(Efield,omega,dz/2,beta,alpha,gam,ts,fR,hR),...
omega,dz/2,beta,alpha,gam,ts,fR,hR);
Uc = rk4ip(Efield,omega,dz,beta,alpha,gam,ts,fR,hR);
error = sqrt(sum(abs(Uf-Uc).^2))/sqrt(sum(abs(Uf).^2));
factor = tol/error;
if error > 2*tol % discard the solution and recaculate
Z_prop = Z_prop-dz; % with a half of old stepsize
else
Efield = 16/15*Uf-1/15*Uc; % improve the accuracy of the solution
end
dz=dz*factor^(1/5);
end
Uout = Efield;
end

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% Function rk4ip
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
function Eout=rk4ip(A0,omega,dz,beta,alpha,gamma,ts,fR,hR)
Efield = A0;
% Aip: interaction presentation
Aip = D_op(Efield,dz/2,beta,alpha,omega);
k1= D_op(N_op(Efield,dz,gamma,ts,fR,hR),dz/2,beta,alpha,omega);
k2= N_op(Aip + k1/2,dz,gamma,ts,fR,hR);
k3= N_op(Aip + k2/2,dz,gamma,ts,fR,hR);
k4= N_op(D_op(Aip+k3,dz/2,beta,alpha,omega),dz,gamma,ts,fR,hR);
Eout = D_op(Aip+k1/6+k2/3+k3/3,dz/2,beta,alpha,omega)+k4/6;
end

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% Function D_op: D operator
% Caculate "exp(hD){A}", in Nonlinear Schrodinger equation
%
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
function Aout = D_op(Ain,dz,beta,alpha,omega)
% Dispersion term in NLSE.
Dw =zeros(1,length(omega));
% Dw = sum(i*beta(i)*(-w)^i/(i!))
for jl = 1 : length(beta)
Dw = Dw - 1i.*beta(jl)/factorial(jl).*(-omega).^(jl);
end
Aout = ifft(exp(dz*(-Dw - alpha/2)).*fft(Ain));
end

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% Function oper_N
% Caculate "hN{A}.*A" in Nonlinear Schrodinger Equation
%
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
function Aout = N_op(Ain,dz,gamma,t_shock,fR,hR_w)
global t
dt = abs((t(2)-t(1)));
% convolution term int_0^inf(hR(t)|A(z,T-t)|^2 dt
% ie. Eq.(7) Hult 2007, or ifft(hR_w.*fft(abs(A)^2))
hR_A2 = ifft(hR_w.*fft(abs(Ain).^2))*dt;
% 1st term
N1 = 1i*gamma*((1-fR).*Ain.*abs(Ain).^2 + fR*Ain.*hR_A2);
% 2nd term.Can be wrote as gradient(F,dt) = gradient(F)/dt
N2 = -gamma*t_shock.*gradient(((1-fR).*Ain.*abs(Ain).^2 ...
+ fR.*Ain.*hR_A2),dt);
Aout = dz*(N1+N2);
end