function [waveAC,SW_Info] = ACside_SW_Wave_Vec(Angle,psInnerAC,fs,Specs)
%AC side three phase waveform calculator
%   Calculate the DC side waveform under Space Vector modulation in one 
%   line frequency period.
%   Specs: Converter specifications
%
%% Calculate waveforms for one switching period  
N_swp=length(psInnerAC);
t_samp       = linspace(0, 1/fs, Specs.Resolution);     %time vector
% t = repmat(t_samp',1,N_swp);
duty        = zeros(1,length(t_samp));     %syncronously operating the AC side three HB
Carrier     = -sawtooth(t_samp*2*pi*fs+pi/2, 0.5);   %carrier for primary side, syncronous modulation

Vc_a       = (cos(Angle/180*pi)+1)*0.5*Specs.Uacpk;
sw  = duty >= Carrier;                      %switch position of three AC HB

sw_a    =   repmat(sw',1,N_swp);
rank_Set= 1:Specs.Resolution;
rank_Mtx= int32(repmat(rank_Set',1,N_swp));
shift_Set=int32(psInnerAC*Specs.Resolution);
shift_Mtx=mod((rank_Mtx-1)-shift_Set,Specs.Resolution)+1;
sw_b   = sw(shift_Mtx);

waveAC.vNode2M.HB_a   = Vc_a*sw_a;
waveAC.sw.HB_a = sw_a;

Vc_b       = (cos(Angle/180*pi+pi)+1)*0.5*Specs.Uacpk;
waveAC.vNode2M.HB_b   = Vc_b*sw_b;
waveAC.sw.HB_b = sw_b;

waveAC.vA2B=waveAC.vNode2M.HB_a-waveAC.vNode2M.HB_b;
SW_Info.vBL = sum(waveAC.vA2B)/Specs.Resolution;
waveAC.vTx=waveAC.vA2B-SW_Info.vBL;

SW_Info.HB_a.ton=ones(1,N_swp);
% SW_Info.HB_a.toff=Specs.Resolution/2*ones(1,Specs.Resolution);
SW_Info.HB_b.ton=mod(1+int32(psInnerAC*Specs.Resolution)-1,Specs.Resolution)+1;
% SW_Info.HB_b.toff=mod(Specs.Resolution/2+int32(psInnerAC*Specs.Resolution)-1,Specs.Resolution)+1;
end