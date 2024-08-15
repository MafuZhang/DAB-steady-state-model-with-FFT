function [waveDC,SW_Info] = Hbridge_Wave_Vec(psInner,Vdc,fs,Specs)
%DC side three phase waveform calculator
%   Calculate the DC side waveform under Space Vector modulation in one 
%   line frequency period.
%
%% Calculate space vector modulation duty cycle look up table
N_swp=length(psInner);
phi     =   psInner;
theta0  =   (1-2*phi)/4;

%% Calculate the switch signal of each phase
waveDC.t    = linspace(0, 1/fs, Specs.Resolution);     %time vector
duty        = zeros(1,length(waveDC.t));     %syncronously operating the AC side three HB
Carrier     = -sawtooth(waveDC.t*2*pi*fs+pi/2, 0.5);   %carrier for primary side, syncronous modulation
sw  = duty >= Carrier;                      %switches state before inner phase shift

rank_Set= 1:Specs.Resolution;
rank_Mtx= int32(repmat(rank_Set',1,N_swp));

shift_Set_a=int32(theta0*Specs.Resolution);
shift_Mtx_a=mod((rank_Mtx-1)-shift_Set_a,Specs.Resolution)+1;
sw_a   = sw(shift_Mtx_a);

shift_Set_b=int32((theta0+phi)*Specs.Resolution);
shift_Mtx_b=mod((rank_Mtx-1)-shift_Set_b,Specs.Resolution)+1;
sw_b   = sw(shift_Mtx_b);

waveDC.sw.HB_a = sw_a;
waveDC.sw.HB_b = sw_b;

waveDC.vNode2P.HB_a   = Vdc*sw_a;
waveDC.vNode2P.HB_b   = Vdc*sw_b;

SW_Info.HB_a.ton    = mod(int32(theta0*Specs.Resolution)-1,Specs.Resolution)+1;
SW_Info.HB_a.toff   = mod(int32((theta0+0.5)*Specs.Resolution)-1,Specs.Resolution)+1;
SW_Info.HB_b.ton    = mod(int32((phi+theta0)*Specs.Resolution)-1,Specs.Resolution)+1;
SW_Info.HB_b.toff   = mod(int32((phi+theta0+0.5)*Specs.Resolution)-1,Specs.Resolution)+1;

waveDC.vTx= waveDC.vNode2P.HB_a-waveDC.vNode2P.HB_b;

end 