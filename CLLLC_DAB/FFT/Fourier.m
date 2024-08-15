function [X_mag,X_ang,f_axes] = Fourier(x,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         Name: Fourier Transform                    	%%%
%%%                         Date: 02.04.17                              %%%  
%%%                         Author: Michael Antivachis                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = length(x);                      %Number of Samples
k = 0:N-1;                      	%FST Harmonic Order

%Maximum Subvector of DFT that Contains Information
if mod(N,2)==0
    k_max = N/2 + 1;                
elseif mod(N,2)==1
    k_max = (N+1)/2;
end            

T = t(end)-t(1);                    %time domain frequency
X = fft(x);                         %DFT Transform
X(k>k_max-1) = '';                	%Keep Spectral Information
k(k>k_max-1) = '';

X_mag      = 2/N.*abs(X);    	    %DFT <---> FST: Magnitutde    
X_mag(1)   = X_mag(1)/2;            %DC Component Correction
X_ang      = angle(X);            	%DFT <---> FST: Angles
f_axes     = k./T;               	%DFT <---> FST: Frequency Axes 


end
