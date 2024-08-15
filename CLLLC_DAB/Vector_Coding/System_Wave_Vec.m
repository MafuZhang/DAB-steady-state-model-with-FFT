function [waveform,Info]  = System_Wave_Vec(VPrim,VSec,fs,ps1,ps2,ps3,Specs)
%System waveform calculator
%   Calculate inductor current of three phase, calculated the local average
%   power, and
%   line frequency period.
%

%% Initialation
N_swp=length(ps1);
Ts     =   1/fs;
Vp     =   VPrim;
Vs     =   VSec;

% Generate the switching pulse of two full bridge
t_deg       =   linspace(0,Ts,Specs.Resolution);
[wave_p,SW_Info_p] = Hbridge_Wave_Vec(ps1,Vp,fs,Specs);
[wave_s,SW_Info_s] = Hbridge_Wave_Vec(ps2,Vs,fs,Specs);
%% Apply mian phase shift to secondary side waveform
rank_Set= 1:Specs.Resolution;
rank_Mtx= int32(repmat(rank_Set',1,N_swp));
shift_Set_AC=int32(ps3*Specs.Resolution);
columeAdd = int32(Specs.Resolution*(0:(N_swp-1)));
shift_Mtx_AC=mod((rank_Mtx-1)-shift_Set_AC,Specs.Resolution)+1+columeAdd;
wave_s.vTx = wave_s.vTx(shift_Mtx_AC);
wave_s.sw.HB_a = wave_s.sw.HB_a(shift_Mtx_AC);
wave_s.sw.HB_b = wave_s.sw.HB_b(shift_Mtx_AC);
% waveAC.vNode2M.HB_a = waveAC.vNode2M.HB_a(shift_Mtx_AC);
% waveAC.vNode2M.HB_b = waveAC.vNode2M.HB_b(shift_Mtx_AC);
% waveAC.sw.HB_a = waveAC.sw.HB_a(shift_Mtx_AC);
% waveAC.sw.HB_b = waveAC.sw.HB_b(shift_Mtx_AC);
SW_Info_s.HB_a.ton=...
    mod(int32(SW_Info_s.HB_a.ton)+int32(ps3*Specs.Resolution)-1,Specs.Resolution)+1;
SW_Info_s.HB_b.ton=...
    mod(int32(SW_Info_s.HB_b.ton)+int32(ps3*Specs.Resolution)-1,Specs.Resolution)+1;

vp      =   wave_p.vTx;
vs      =   wave_s.vTx;
waveform.vp=vp;
waveform.vs=vs;
waveform.gp1 = wave_p.sw.HB_a;
waveform.gp2 = wave_p.sw.HB_b;
waveform.gs1 = wave_s.sw.HB_a;
waveform.gs2 = wave_s.sw.HB_b;
% Frequency domain axis initiation
[Vp_mag, Vp_ang, fh] = VEC_Fourier(vp, t_deg);                                %Fourier tranformation for Resonant tank votlage source
Wh      = 2*pi.*fh;                                                         %Angular Velocity
Vp_fq = Vp_mag.*exp(1i.*Vp_ang);                                         % Frequency domain expression of Resonant tank votlage source

[Vs_mag, Vs_ang, ~] = VEC_Fourier(vs, t_deg);
Vs_fq = Vs_mag.*exp(1i.*Vs_ang); 

%% Waveforms in frequency domain 
% Three Y parameter in T circuit
Lrs_p = Specs.Lrs*Specs.n*Specs.n;
Crs_p = Specs.Crs/(Specs.n*Specs.n);
Rrs_p = Specs.R*Specs.n*Specs.n;

Za = 1./(Specs.Crp*(1i.*Wh))+Specs.Lrp*(1i.*Wh)+Specs.R.*ones(size(Wh));
Zb = 1./(Crs_p*(1i.*Wh))+Lrs_p*(1i.*Wh)+Rrs_p.*ones(size(Wh));
Zc = Specs.Lm*(1i.*Wh);

Zab = Za+Zb+Za.*Zb./Zc;
Zbc = Zb+Zc+Zb.*Zc./Za;
Zca = Zc+Za+Zc.*Za./Zb;

Yab=1./Zab;
Ybc=1./Zbc;
Yca=1./Zca;

Y11 =   Yca+Yab;
Y12 =   -Yab;
Y21 =   -Yab;
Y22 =   Yab+Ybc;

%circuit cal
ILrp_fq     =   Vp_fq.*Y11+Vs_fq.*Y12;
ILrs_p_fq   =   Vp_fq.*Y21+Vs_fq.*Y22;
ILm_fq      =   ILrp_fq+ILrs_p_fq;
% VLrp_fq     =   Specs.Lrp.*ILrp_fq.*(1i.*Wh);
% VLrs_fq     =   Lrs_p.*ILrs_p_fq.*(1i.*Wh);
VLm_fq      =   Specs.Lm.*(1i.*Wh).*ILm_fq;
% VCrp_fq     =   1./(Specs.Crp.*(1i.*Wh)).*ILrp_fq;
% VCrs_p_fq   =   1./(Crs_p.*(1i.*Wh)).*ILrs_p_fq;

% iLrp
ILrp_fq_cal = ILrp_fq;
ILrp_fq_cal(1,:)  =   0;
ILrp_mag = abs(ILrp_fq_cal);
ILrp_ang = angle(ILrp_fq_cal);
iLrp     =   VEC_Fourier_Inverse(ILrp_mag,ILrp_ang,t_deg);
waveform.iLrp    =   iLrp;

% iLrs
ILrs_p_fq_cal = ILrs_p_fq;
ILrs_p_fq_cal(1,:)  =   0;
ILrs_mag = abs(ILrs_p_fq_cal);
ILrs_ang = angle(ILrs_p_fq_cal);
iLrs     =   VEC_Fourier_Inverse(ILrs_mag,ILrs_ang,t_deg);
waveform.iLrs    =   iLrs;

% % iLm
% ILm_fq_cal = ILm_fq;
% ILm_fq_cal(1,:)   =   0;
% ILm_mag = abs(ILm_fq_cal);
% ILm_ang = angle(ILm_fq_cal);
% iLm     =   VEC_Fourier_Inverse(ILm_mag,ILm_ang,t_deg);
% 
% %vLrp
% VLrp_fq_cal=VLrp_fq;
% VLrp_fq_cal(1,:) = 0;
% VLrp_mag = abs(VLrp_fq_cal);
% VLrp_ang = angle(VLrp_fq_cal);
% vLrp     =   VEC_Fourier_Inverse(VLrp_mag,VLrp_ang,t_deg);
% 
% %vLrs
% VLrs_fq_cal=VLrs_fq;
% VLrs_fq_cal(1,:) = 0;
% VLrs_mag = abs(VLrs_fq_cal);
% VLrs_ang = angle(VLrs_fq_cal);
% vLrs     =   VEC_Fourier_Inverse(VLrs_mag,VLrs_ang,t_deg);
% 
%vLm
VLm_fq_cal=VLm_fq;
VLm_fq_cal(1,:) = 0;
VLm_mag = abs(VLm_fq_cal);
VLm_ang = angle(VLm_fq_cal);
vLm     =   VEC_Fourier_Inverse(VLm_mag,VLm_ang,t_deg);
waveform.vLm=vLm;
% %vCrp
% VCrp_fq_cal=VCrp_fq;
% VCrp_fq_cal(1,:) = 0;
% VCrp_mag = abs(VCrp_fq_cal);
% VCrp_ang = angle(VCrp_fq_cal);
% vCrp     =   VEC_Fourier_Inverse(VCrp_mag,VCrp_ang,t_deg);
% 
% %vCrs
% VCrs_p_fq_cal=VCrs_p_fq;
% VCrs_p_fq_cal(1,:) = 0;
% VCrs_mag = abs(VCrs_p_fq_cal);
% VCrs_ang = angle(VCrs_p_fq_cal);
% vCrs     =   VEC_Fourier_Inverse(VCrs_mag,VCrs_ang,t_deg);

%power
P_instant   =   iLrp.*vp;
%% function output

% Switching current
% AC side
%switching point location matrix
%[swLocMtx_ACa,swLocMtx_ACb,swLocMtx_DCa,swLocMtx_DCb] = deal(zeros(size(iLr)));
swLocIndx_pa = sub2ind(size(iLrp),SW_Info_p.HB_a.ton,1:N_swp);
swLocIndx_pb = sub2ind(size(iLrp),SW_Info_p.HB_b.ton,1:N_swp);

swLocIndx_sa = sub2ind(size(iLrs),SW_Info_s.HB_a.ton,1:N_swp);
swLocIndx_sb = sub2ind(size(iLrs),SW_Info_s.HB_b.ton,1:N_swp);

% Switching current
Info.SW_Info_p.HB_a.Ion  = iLrp(swLocIndx_pa);
Info.SW_Info_p.HB_b.Ion  = -iLrp(swLocIndx_pb);

Info.SW_Info_s.HB_a.Ion  = iLrs(swLocIndx_sa);
Info.SW_Info_s.HB_b.Ion  = -iLrs(swLocIndx_sb);


% % waveform result output port
% waveform.vTx.ACside    =   waveAC.vTx;
% waveform.vTx.DCside    =   waveDC.vTx;
% waveform.vTx.System    =   waveAC.vTx-waveDC.vTx;
% waveform.P_instant     =   P_instant;

% Other information result output port

Info.P_LocalAvg     =   sum(P_instant)/Specs.Resolution;
Info.Irms_p           =   rms(iLrp);
Info.Irms_s           =   rms(iLrs);
Info.flagZVS        =   ...
    ((Info.SW_Info_s.HB_a.Ion<0)&(Info.SW_Info_s.HB_b.Ion<0))&...
    ((Info.SW_Info_p.HB_a.Ion<0)&(Info.SW_Info_p.HB_b.Ion<0));

end