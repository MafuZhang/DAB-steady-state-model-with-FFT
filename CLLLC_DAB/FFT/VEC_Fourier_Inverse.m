function [x] = VEC_Fourier_Inverse(X_mag,X_ang,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%              Name: Vector Form Fourier Inverse Transform           	%%%
%%%              Date: 05.07.2023                                       %%%
%%%              Author: Mafu Zhang @ UT Austin                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = length(t);                                  %Length of Time Vector

X = N/2.*X_mag.*exp(1i.*X_ang);              	%FST <---> DFT: Phasor

%Reconstruct Complete DFT Signal
if mod(N,2)==0
    X = [2*X(1,:); X(2:end,:); flip(conj(X(2:end-1,:)))];
elseif mod(N,2)==1
    X = [2*X(1), X(2:end), flip(conj(X(2:end)))];
end

x = real(ifft(X));    %Inverse DFT Transform

%Check Arguments Validity
if abs(N - 2*size(X_mag,1))>2
    error('Inverse Fourier');
end

end

