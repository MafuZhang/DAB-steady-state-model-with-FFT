function [X_mag,X_ang,f_axes] = VEC_Fourier(x,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         Name: Vector Form Fourier Transform        	%%%
%%%                         Date: 05.07.2023                            %%%  
%%%                         Author: Mafu Zhang @ UT Austin              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N,Nk] = size(x);                      %Number of Samples
k_temp = 0:N-1;                      	%FST Harmonic Order
k=repmat(k_temp',1,Nk);
%Maximum Subvector of DFT that Contains Information
if mod(N,2)==0
    k_max = N/2 + 1;                
elseif mod(N,2)==1
    k_max = (N+1)/2;
end            

T = t(end)-t(1);                    %time domain frequency
X = fft(x);                         %DFT Transform
X(k_temp>k_max-1,:) = '';                	%Keep Spectral Information
k(k_temp>k_max-1,:) = '';

X_mag      = 2/N.*abs(X);    	    %DFT <---> FST: Magnitutde    
X_mag(1,:)   = X_mag(1,:)/2;            %DC Component Correction
X_ang      = angle(X);            	%DFT <---> FST: Angles
f_axes     = k./T;               	%DFT <---> FST: Frequency Axes 


end
