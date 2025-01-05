function [R,D,G]=mat_build(i, sig, th_disc, n)
%n size of NtD matrix
%sig- is a vector containing conductivity values at
%space discretization points of current ring
%n_th- contains the number of discretization points 
%at all the rings
%th_disc- contains the angular discretization values at all the rings



if min(sig)==0
    epss=10^(-6);
else
    epss=10^(-8)*min(sig);
end

rho=1./(sig + epss);%1/conductivity values





m=length(sig);%length of given conductivity values

f_coe=f_coeff(sig,th_disc{i});
%compute the fourier coefficients of discretized angular conductivity values

f_trunccoe=[f_coe(1),f_coe(2:1+floor(m/2)),f_coe(m+2:m+1+floor(m/2))];
%take half of the coefficients (Nyquest sampling)

[fr_coe]=f_coeff(rho,th_disc{i});
%for 1/conductivity


fr_trunccoe=[fr_coe(1),fr_coe(2:1+floor(m/2)),fr_coe(m+2:m+1+floor(m/2))];
%Nyquest sampling


R=build_R(fr_trunccoe,n);
%R matrix in matrix layer stripping equation 
D=build_D(n);
G=build_G(f_trunccoe,n);
%G matrix in matrix layer stripping equation 

end