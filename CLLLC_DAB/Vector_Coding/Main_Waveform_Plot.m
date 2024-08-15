%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     CLLLC DAB waveform plot                         %%%
%%%       Author: Mafu Zhang                                            %%%
%%%       latest modified Date: April.1.2022                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
%%
addpath(genpath([cd(cd('..\..\..\')),'\CLLLC_AC-DC_DAB']))
addpath(genpath([cd(cd('..\..\')),'\Optimization']))
addpath(genpath([cd(cd('..\..\')),'\Optimization\Optimizer']))
addpath(genpath([cd(cd('..\..\')),'\Optimization\components']))
% addpath(genpath([cd(cd('..\')),'\Optimization\components\Semiconductor']))
load('Device_fake_high_Crr.mat');
load('Transformers_new.mat');
%% System Specifications
Specs.inspScope     =   0;      %0:swiching period; 1:line period

%% Simulation control
% figure sitting
fig.on = 1;                 %plot results (0/1)
% analyzer or optimizer
Specs.AorO = 0;  % modes : 0:line period waveform with buck modulation
%         1:Sweeper of triple phase shift
%         2:TBD: Modulation Optimizer
Specs.SweepOrPlot = 1; %=1 will enter vector computation mode, for power sweep.
% circuit specifications

Specs.SWp                 =   Device(7);
Specs.SWs                 =   Device(7);

% Circuit para
Specs.Lr    =   30e-6;
Specs.n     =   1;
Specs       =   getSpecs(Specs);
fs=Specs.fsw;
%% ZVS and loss analysis (对每个角度sweep，记录并可视化三个量，一个是power，另一个是zvs，再一个是conduction loss占比)
%control varibles domain
tic
Ang_swp=5;
ps1 = 0.3;
ps2 = 0.4;
ps3 = 0.1;

%opeartion data obtain

deg=45;
[waveSys_swp,InfoSys_swp]  = System_Wave_Vec(deg,fs,ps1,ps2,ps3,Specs);

%%

fig1=figure();
title('Voltages and iLr after iteration converge')
sf1=subplot(2,1,1);
plot(waveform.t,waveform.Vai,'k-');
hold on
plot(waveform.t,waveform.Vbi,'b-');
plot(waveform.t,waveform.vLm,'r--');
ylim([-500,500])
grid on
sf2=subplot(2,1,2);
plot(waveform.t,waveform.iLrp,'k-');
hold on
plot(waveform.t,waveform.iLm,'r--');
plot(waveform.t,waveform.iLrs,'b-')
ylim([-70,70])
grid on
movegui(fig1,[600 300]);
linkaxes([sf1,sf2],'x')
datacursormode(fig1,'on')
% 
% fig2=figure();
% title('Cr voltage')
% plot(waveform.t,waveform.vCrp)

