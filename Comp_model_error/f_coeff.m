%This function computes the Fourier coefficients (Sine & Cosine)
%for a given
%sig- Conductivity angular discritization
%th- Angular discretization

function f_coe=f_coeff(sig,th)
n=length(th);
 for i=1:length(th)
     %compute the cosine fourier coefficients
     f_c(i)=(2/(n))*(sum(sig(1:end).*cos(i*th(1:end))));
     %compute the sine fourier coefficients
     f_s(i)=(2/(n))*(sum(sig(1:end).*sin(i*th(1:end))));    
 end
 %compute zero freqency coefficient
 f_0=(1/(n))*sum(sig(1:end));
 %Arrange the coefficients in the order of [cosine,sine] frequencies
 f_coe=[f_0,f_c,f_s];
end
