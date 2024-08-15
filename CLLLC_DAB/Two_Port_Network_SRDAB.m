%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     CLLLC DAB                                       %%%
%%%       Author: Mafu Zhang                                            %%%
%%%       latest modified Date: April.1.2022                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
%%
addpath('FFT\')
addpath('Vector_Coding\')

%% Parameters
% Simulation parameters
Specs.k=1;                                                % Plot the lossmap curve fitting or not
% Ohter para
Specs               =   getSpecs(Specs);
%% Input
d1=0.79;
d2=0.69;
d3=0.15;
fsw = Specs.fsw;
Deg=90;


%% Initialation
OperatingPoint.vac     =   abs(Specs.Vac_amp*sin(2*pi*Specs.fac*Deg/(360*Specs.fac)));
OperatingPoint.iac     =   abs(Specs.Iac_amp*sin(2*pi*Specs.fac*Deg/(360*Specs.fac)));     % Input current level at this instance
OperatingPoint.vdc     =   Specs.Vdc;
OperatingPoint.Ppfc    =   abs(Specs.Vac_amp*Specs.Iac_amp*...
        sin(2*pi*Specs.fac*Deg/(360*Specs.fac))^2);  
%% Iterating waveform generator
[Info,waveform]  = srDABinfo_ftps_FFT_Two_Port(d1,d2,d3,fsw,OperatingPoint,Specs);

%% plot

fig1=figure();
title('Voltages and iLr after iteration converge')
sf1=subplot(2,1,1);
yyaxis left
plot(waveform.t,waveform.Vai);
hold on
plot(waveform.t,waveform.Vbi);
ylim([-240,240])
legend('Primary side pulse','Secondary side pulse')
grid on
sf2=subplot(2,1,2);
plot(waveform.t,waveform.iLrp);
hold on
plot(waveform.t,-waveform.iLrs)
plot(waveform.t,-waveform.iLm);
legend('primary side current','secondary side current','Magnitizing current')

grid on
movegui(fig1,[600 300]);
linkaxes([sf1,sf2],'x')
datacursormode(fig1,'on')
xlim([0 max(waveform.t)])
% 
% fig2=figure();
% title('Cr voltage')
% plot(waveform.t,waveform.vCrp)
%%
 power=Info.Power
 Isw_p1 = Info.Isw.dc1
 Isw_p2 = Info.Isw.dc2
 Isw_s1 = Info.Isw.ac1
 Isw_s2 = Info.Isw.ac2
 Irms_p = Info.Irms_p
 Irms_s = Info.Irms_s