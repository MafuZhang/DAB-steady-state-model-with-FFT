function [Ptotal,Pp_cond,Pp_sw,Ps_cond,Ps_sw] = Loss_calculation(fsw,DATA_V,DATA_Irms,DATA_Isw,Specs)
%Semiconductor loss calculation
%   Detailed explanation goes here
Irms_s=Irms_p*Specs.n;
%% Switch losses
%% AC side
% Conduction loss
Pp_cond     =   DATA_Irms.^2*(2*Specs.SWp.Rds);
% Turn off loss 
E_ACa     =   Turn_off_loss(DATA_V.ACa,DATA_Isw.ACa,Specs.SWp);
E_ACb     =   Turn_off_loss(DATA_V.ACb,DATA_Isw.ACb,Specs.SWp);

P_ACa      =   E_ACa*fsw;
P_ACb      =   E_ACb*fsw;
Pp_sw           =   (P_ACa+P_ACb)*2;
%% DC side
% Conduction loss
Ps_cond         =   DATA_Irms^2*((2*Specs.SWs.Rds)/Specs.n_SWs_parallel);
% Turn off loss          
E_DCa      =   Turn_off_loss(Vdc,Isw.dc1,Specs.SWs);
E_DCb      =   Turn_off_loss(Vdc,Isw.dc2,Specs.SWs);

P_DCa      =   E_DCa*fsw;
P_DCb      =   E_DCb*fsw;
Ps_sw          =   (P_DCa+P_DCb)*2;


%% total loss
Ptotal = Pp_sw+Ps_sw+Pp_cond+Ps_cond;

end