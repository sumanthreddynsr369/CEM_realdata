%this program generates the sample of NtD matrix diagonal entries for each
%and every conductivity sample at the boundaries of space discretizaion
%points

function [NtD_samp]=NtDmatrix_sample_radial(n_freq,Sig_samp,a)


%Sig_samp is a matrix with rows asdifferent samples 
%and columns as conductivity value at space discretization points

% a--- is space discretization
n_space=size(a,2);%number of space discretization points including origin (n_s+1)
n_samp=size(Sig_samp,1);%number of conductivity samples
nn=n_freq;%number of sine and cosine frequencies 

M=@ (t) (1-t)./(1+t); %mobius map

%in which NtD matrix is represented
RR=(1./[1:1:nn]');
Y_in=zeros(nn,n_samp);
for i=1:n_samp
NtD_0=(1/Sig_samp(i,1))*RR;%computing the initial NtD diagonal entries for each 
%created sample
Y_in(:,i)=diag(1:1:nn)*NtD_0;
end

NtD_samp=zeros(nn*n_samp,n_space-1);
%this matrix saves all the NtD components at each discretization points 
%as concatenated columns for each sample 

for j=1:n_samp
  Y_old=Y_in(:,j);
 for i=2:n_space-1    
   Y_new =(1/Sig_samp(j,i))*M(ydiag(a(i)/a(i+1),nn)*M(Sig_samp(j,i).*Y_old));
   Y_old=Y_new;  
    W(i-1,:)=Y_new;
 end
 NtD_samp(nn*(j-1) + 1 :nn*j,2:end)=W';
end
% 
NtD_samp(:,1)=Y_in(:);

end

