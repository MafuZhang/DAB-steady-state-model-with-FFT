%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%             CLLLC SRDAB power sweeper                               %%%
%%%             Author: Mafu Zhang @ UTaustin SPEC                      %%%
%%%             latest modified Date: 29.7.2023                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
%%
addpath(genpath([cd(cd('..\..\..\')),'\CLLLC_AC-DC_DAB']))
addpath(genpath([cd(cd('..\..\')),'\Optimization']))
addpath(genpath([cd(cd('..\..\')),'\Auxiliary_Func']))
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
Specs               =   getSpecs(Specs);
fs=Specs.fsw;
%% ZVS and loss analysis (对每个角度sweep，记录并可视化三个量，一个是power，另一个是zvs，再一个是conduction loss占比)
%control varibles domain
tic
Ang_swp=5;
% ps1 = linspace(0.001, 0.497, 20);
% ps2 = linspace(0.001, 0.497, 20);
% ps3 = linspace(0.001, 0.494, 50);
ps1 = linspace(0, 0.5, 2);
ps2 = linspace(0, 0.5, 2);
ps3 = linspace(0, 0.5, 1);
[PS1,PS2,PS3] = meshgrid(ps1,ps2,ps3);   % generate phase shift matrix

ps1_HR = linspace(0, 0.5, 1);    % Higher resolution by interpolation
ps2_HR = linspace(0, 0.5, 1);
ps3_HR = linspace(0, 0.5, 2);
[PS1_m,PS2_m,PS3_m] = meshgrid(ps1_HR,ps2_HR,ps3_HR);

degCount=1;
Deg=linspace(0,90,degCount);
[DATA_power,DATA_Irms,DATA_ZVS,DATA_Ion_ACa,DATA_Ion_ACb,DATA_Ion_DCa,DATA_Ion_DCb]...
    =deal(zeros(numel(PS3),degCount));
N_Sample=numel(PS3);
ps1_vec=PS1(:)';
ps2_vec=PS2(:)';
ps3_vec=PS3(:)';
%opeartion data obtain
% for i=1:degCount 
%     Deg(i)
%     [waveSys_swp,InfoSys_swp]  = System_Wave_Vec(Deg(i),fs,ps1_vec,ps2_vec,ps3_vec,Specs);
%     DATA_power(:,i)=InfoSys_swp.P_LocalAvg;
%     DATA_Irms(:,i)=InfoSys_swp.Irmsp;
%     DATA_ZVS(:,i)=InfoSys_swp.flagZVS;
%     DATA_Ion_ACa(:,i)=InfoSys_swp.SW_Info_s.HB_a.Ion;
%     DATA_Ion_ACb(:,i)=InfoSys_swp.SW_Info_s.HB_b.Ion;
%     DATA_Ion_DCa(:,i)=InfoSys_swp.SW_Info_p.HB_a.Ion;
%     DATA_Ion_DCb(:,i)=InfoSys_swp.SW_Info_p.HB_b.Ion;
% end
[zvsSweep,powerSweep]...
    =deal(zeros(numel(PS3),degCount));
volCount = 1;
Vs_sweep=linspace(200,200,volCount);
for i=1:volCount 
    i
    [waveSys_swp,InfoSys_swp]  = System_Wave_Vec(Specs.Vdc,Vs_sweep(i),fs,ps1_vec,ps2_vec,ps3_vec,Specs);
    zvsSweep(:,i)=InfoSys_swp.flagZVS;
    powerSweep(:,i)=InfoSys_swp.P_LocalAvg;
end
% %% ZVS range plot

zvsSweep_nan=zvsSweep(:,1);
index_NZVS=find(~zvsSweep_nan);
ps1_vec_NaN=ps1_vec';
ps1_vec_NaN(index_NZVS)=[];
ps2_vec_NaN=ps2_vec';
ps2_vec_NaN(index_NZVS)=[];
ps3_vec_NaN=ps3_vec';
ps3_vec_NaN(index_NZVS)=[];
powerSweep_NaN=powerSweep(:,1);
powerSweep_NaN(index_NZVS)=[];

zvsFig=figure();
hAxes = axes;
scHandle=scatter3(ps1_vec_NaN,ps2_vec_NaN,ps3_vec_NaN,'filled');
set(scHandle,'CData',powerSweep_NaN);
scHandle.SizeData=7;
% plot3(ps1_vec_NaN,ps2_vec_NaN,ps3_vec_NaN,'O-')
xlim([0,0.5])
ylim([0,0.5])
zlim([0,0.5])
xlabel('DC side phase shift');
ylabel('AC side phase shift');
zlabel('Main phase shift');

zvsFig.Position=[1000 606 732 632];
pos = get(gcf, 'Position'); %// gives x left, y bottom, width, height
vievAng=[-60.8689 11.6179];
hAxes.View=vievAng;
zvsFig.Position=[1000 606 732 632];

l = light('Position',[-1.2 -0.7 0.5]);
lightangle(-45,30)
shading interp
colorbar EastOutside

%% power plot
% plot raw data
DATA_Power_deg=powerSweep(:,1);
DATA_Power_mtx=reshape(DATA_Power_deg,size(PS3));
DCslice = [0.499];
ACslice = [0.499];   
Mainslice = [0.15];
%power
figure()
pSlice=slice(PS1,PS2,PS3,DATA_Power_mtx,DCslice,ACslice,Mainslice);
xlabel('DC side phase shift');
ylabel('AC side phase shift');
zlabel('Main phase shift');
% axis vis3d
xlim([0,0.5])
ylim([0,0.5])
zlim([0,0.5])
lighting phong
l = light('Position',[-1.2 -0.7 0.5])
lightangle(-45,30)
 %none%phong%gouraud
shading interp
colorbar EastOutside




%%

% %% plot 
% nn=3000;
% [waveSys_swp,InfoSys_swp]  = System_Wave_Vec(Specs.Vdc,Specs.Vdc,fs,ps1_vec,ps2_vec,ps3_vec,Specs);
% figure
% hold on
% plot(waveSys_swp.vp(:,nn))
% plot(waveSys_swp.vs(:,nn))
% plot(waveSys_swp.iLrp(:,nn))
% zvsflag=InfoSys_swp.flagZVS(nn)
% pave=sum(waveSys_swp.iLrp(:,nn).*waveSys_swp.vp(:,nn))/Specs.Resolution
% power_local=InfoSys_swp.P_LocalAvg(nn)

% %% ZVS-power%volt
% zvspowerSweep = powerSweep;
% indexNzvs=find(~zvsSweep);
% zvspowerSweep(indexNzvs)=NaN;
% XaxisMat = repmat(Vs_sweep,length(PS3(:)),1);
% 
% figure()
% hold on
% for i=volCount
%     scatter(XaxisMat,zvspowerSweep)
% end

%% temp plots
 %%nn=4758; %figure compare
 % 
 % nn=593 %figure 400DAB A
 % nn=82  %figure 400DAB B
 % nn=545 %figure 400DAB C

 nn=693 %figure 200DAB A
 nn=174  %figure 200DAB A
% nn=545 %figure 200DAB C



i
 figure()
 for i = 1:4
     switch i
         case 1
             sp(1)=subplot(6, 1, 1)
             hold on
             plot(waveSys_swp.gp1(:,nn),'b')
             plot(waveSys_swp.gp2(:,nn),'b--')
             ylim([-0.5,1.5])
         case 2
             sp(2)=subplot(6, 1, 2)
             hold on
             plot(waveSys_swp.gs1(:,nn),'r')
             plot(waveSys_swp.gs2(:,nn),'r--')
             ylim([-0.5,1.5])
         case 3
             sp(3)=subplot(6, 1, 3:4)
             hold on
             plot(waveSys_swp.vp(:,nn))
             plot(waveSys_swp.vs(:,nn))
             plot(waveSys_swp.vLm(:,nn))
             ylim([-500,500])
         case 4
             sp(4)=subplot(6, 1, 5:6)
             hold on
             plot(waveSys_swp.iLrp(:,nn))
             plot(-waveSys_swp.iLrs(:,nn))
             plot(waveSys_swp.iLrp(:,nn)+waveSys_swp.iLrs(:,nn))
             ylim([-70 70])
         case 5
     end

     if i == 3 || i == 4
         set(gca, 'xtick', []);

         p = get(gca, 'Position');
         % Increase the height of the first and third subplots by 10%
         p_diff = p(4) * 0.06;
         % Increase the height of the subplot, but this will keep the
         % bottom in the same place
         p(4) = p(4) + p_diff;
         % So also move the subplot down to decrease the gap to the next
         % one.
         p(2) = p(2) - p_diff;
         set(gca, 'Position', p);
     end
 end
 linkaxes(sp,'x')
 xlim([0 Specs.Resolution])
 Phsft1 = ps1_vec(nn)
 Phsft2 = ps2_vec(nn)
 Phsft3 = ps3_vec(nn)
 power=InfoSys_swp.P_LocalAvg(nn)
 Isw_p1 = InfoSys_swp.SW_Info_p.HB_a.Ion(nn)
 Isw_p2 = InfoSys_swp.SW_Info_p.HB_b.Ion(nn)
 Isw_s1 = InfoSys_swp.SW_Info_s.HB_a.Ion(nn)
 Isw_s2 = InfoSys_swp.SW_Info_s.HB_b.Ion(nn)
 Irms_p = InfoSys_swp.Irms_p(nn)
 Irms_s = InfoSys_swp.Irms_s(nn)
 %%
%  % loss and efficiency calculation
%  DATA_Isw.ACa=DATA_Ion_ACa;
%  DATA_Isw.ACb=DATA_Ion_ACb;
%  DATA_Isw.DCa=DATA_Ion_DCa;
%  DATA_Isw.DCb=DATA_Ion_DCb;
%  VACa = (cos(Deg/180*pi)+1)*0.5*Specs.Uacpk;
%  VACb = (cos(Deg/180*pi+pi)+1)*0.5*Specs.Uacpk;
%  VDC = Specs.Udc*ones(1,degCount);
%  DATA_V.DC = repmat(VDC,N_Sample,1);
%  DATA_V.ACa =repmat(VACa,N_Sample,1);
% DATA_V.ACb = repmat(VACb,N_Sample,1);
% 
% [Ptotal,Pp_cond,Pp_sw,Ps_cond,Ps_sw] = Loss_calculation(fs,DATA_V,DATA_Irms,DATA_Isw,Specs);
% DATA_Eff = (abs(DATA_power)-Ptotal)./abs(DATA_power);
% 
% %% interpolation for optimization
% 
% k0=[PS1_m(:),PS2_m(:),PS3_m(:)];
% k=[PS1(:),PS2(:),PS3(:)];
% 
% for i = 1:degCount
%     Deg(i)
%     Power_HR=griddatan(k,DATA_power(:,i),k0);
%     Eff_HR = griddatan(k,DATA_Eff(:,i),k0);
%     DATA_P.(sprintf('degCount%i',Deg(i)))=reshape(Power_HR,size(PS1_m));
%     DATA_E.(sprintf('degCount%i',Deg(i)))=reshape(Eff_HR,size(PS1_m));
% end
% %%
% % plot raw data
% DATA_Power_deg=DATA_power(:,1);
% DATA_Power_mtx=reshape(DATA_Power_deg,size(PS3));
% DCslice = [0.25 0.4];
% ACslice = [-0.1 -0.4];   
% Mainslice = [];
% %power
% figure()
% slice(PS1,PS2,PS3,DATA_Power_mtx,DCslice,ACslice,Mainslice)
% xlabel('DC side phase shift');
% ylabel('AC side phase shift');
% zlabel('Main phase shift');
% %efficiency
% DATA_Eff_deg=DATA_Eff(:,1);
% DATA_Eff_mtx=reshape(DATA_Eff_deg,size(PS3));
% 
% figure()
% slice(PS1,PS2,PS3,DATA_Eff_mtx,DCslice,ACslice,Mainslice)
% xlabel('DC side phase shift');
% ylabel('AC side phase shift');
% zlabel('Main phase shift');
% 
% % plot interpolared data
% DCslice = [0.25 0.4];
% ACslice = [-0.1 -0.4];   
% Mainslice = [];
% %power
% figure()
% slice(PS1_m,PS2_m,PS3_m,DATA_P.degCount2,DCslice,ACslice,Mainslice)
% xlabel('DC side phase shift');
% ylabel('AC side phase shift');
% zlabel('Main phase shift');
% %efficiency
% figure()
% slice(PS1_m,PS2_m,PS3_m,DATA_E.degCount2,DCslice,ACslice,Mainslice)
% xlabel('DC side phase shift');
% ylabel('AC side phase shift');
% zlabel('Main phase shift');
% %% Optimization
% Ind_maxEff = ones(degCount,1);
% Sub_maxEff = ones(degCount,3);
% max_eff = zeros(degCount,1);
% for i=1:degCount
% 
%     PFC     =   0.2*abs(Specs.Uacpk*Specs.Iacpk*cos(2*pi*Deg(i)/360)^2)+100;
%     LowerBound = max(0.92*PFC,PFC-50);
%     UpperBound = min(1.08*PFC,PFC+50);
% 
%     Ind_pfcScreen=find(DATA_P.(sprintf('degCount%i',Deg(i)))<=UpperBound & DATA_P.(sprintf('degCount%i',Deg(i)))>=LowerBound);
% 
%     Power_pfcScreen.(sprintf('degCount%i',Deg(i))) = DATA_P.(sprintf('degCount%i',Deg(i)))(Ind_pfcScreen);
%     Eff_pfcScreen.(sprintf('degCount%i',Deg(i))) = DATA_E.(sprintf('degCount%i',Deg(i)))(Ind_pfcScreen);
% 
%     [~,Ind_maxEff(i)]=max(Eff_pfcScreen.(sprintf('degCount%i',Deg(i))));
% 
%     ind_maxEff = Ind_pfcScreen(Ind_maxEff(i));
%     max_eff(i)=DATA_E.(sprintf('degCount%i',Deg(i)))(ind_maxEff);
% 
%     Sub_maxEff(i,:)=[PS1_m(ind_maxEff),PS2_m(ind_maxEff),PS3_m(ind_maxEff)];
% 
%     [Ind1, Ind2, Ind3]=ind2sub(size(PS2_m),Ind_pfcScreen);
%     Ind_RawPFC.(sprintf('degCount%i',Deg(i)))=[Ind1 Ind2 Ind3];
%     % Ind_maxEff.
% end
% Sub_maxEff
% max_eff
% figure()
% for i=1:degCount
%     hold on
%     len_screen = length(Power_pfcScreen.(sprintf('degCount%i',Deg(i))));
%     scatter(i*ones(1,len_screen),Power_pfcScreen.(sprintf('degCount%i',Deg(i))));
% end
% figure()
% plot(Sub_maxEff)
