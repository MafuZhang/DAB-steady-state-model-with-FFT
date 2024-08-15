function [x] = Fourier_Inverse(X_mag,X_ang,t,f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         Name: Fourier Inverse Transform           	%%%
%%%                         Date: 02.04.17                              %%%  
%%%                         Author: Michael Antivachis                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('f')==0
    
    N = length(t);                                  %Length of Time Vector
    
    X = N/2.*X_mag.*exp(1i.*X_ang);              	%FST <---> DFT: Phasor
    
    %Reconstruct Complete DFT Signal
    if mod(N,2)==0
        X = [2*X(1), X(2:end-1), X(end), flip(conj(X(2:end-1)))];
    elseif mod(N,2)==1
        X = [2*X(1), X(2:end), flip(conj(X(2:end)))];
    end
    
    x = real(ifft(X));    %Inverse DFT Transform
    
    %Check Arguments Validity
    if abs(N - 2*length(X_mag))>2
        error('Inverse Fourier');
    end

elseif exist('f')==1
    
    x = zeros(1,length(t));
    for i=1:length(f)
        x = x + X_mag(i).*cos(2*pi*f(i).*t + X_ang(i));
    end
    
end



end

