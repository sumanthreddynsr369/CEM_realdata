%this function generates the radially symmetric smooth
%log(conductivity/background conductivity)
%samples
%%input
%n_num-number of discretization points in space
%num_sample- number of required prior samples
%gamma-variance in sample
%L-second order smoothness prior matrix

%%output
%Xi- matrix with rows as xi- value at space discretization points
%and columns as different samples
%conductivity
%Xi=log(sigma/sigma_back)
function [Xi,a]=xi_sample_radial(n_num,num_sample,a_in,a_fin,gamma,alp)

[a,L]=fin_diff_mat_1d(a_in,a_fin,n_num);
Xi=zeros(num_sample,n_num+1);
for i=1:num_sample
    w=randn(n_num+1,1);%draw sample from white noise
    x=(L + (1/alp^2).*eye(size(L)))\(gamma*w);%generate the smooth sample
    %x=(L)\(gamma*w);
    Xi(i,:)=x;
end
end
    
    
