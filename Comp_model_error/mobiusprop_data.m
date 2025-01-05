
%clear all;
clc;




load VoltageMesh_12 Vertices Topol 
load Real_data_matrix_homo sig_CEM
V=Vertices;%used to generate a random sample 
n_samp=250;%number of samples
n_space = 30;       %number of radial discretization points (excluding the center)
n_freq = 30;        %number of sine/cosine frequencies in the NtD matrix, hence size of NtD matrix
f=1;%sparsity controller for prior space discretization
          %as f--->0.1, the prior space discretization becomes coarser.   
  visualize='no';%visualize the prior space discretization
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%radial discretization 
a_in=0;             %center of the disc
a_fin=1;            %
a=[a_in:(a_fin-a_in)/n_space:a_fin];%radial discretization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
[Vertices,theta,n_th]=prior_space_discre(a,f,visualize);%this program generates the
%discretization on disc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Construct a test conductivity distribution

 visualize_cond='no';%yes/no, to visualize the true conductivity
 
 n=2*n_freq;
 n_NtD=n;%size of NtD matrix
  I_2=eye(2*n_NtD);%will be used in second order mobius update
  S=zeros(n_NtD);%will be used in second order mobius update
 sig_bg=1.5;
 for nn=1:250
 sig_true=true_cond(Vertices, visualize_cond,V,nn);%this file creates a cell structure of true conductivity 
      %values from the outer most ring to the inner ring
  sig_true_init=sig_true{end};%true conductivity values in the deepest ring
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %constructing the initial NtD matrix on the boundary of innermost ring
 RR=[1./[1:1:n_freq]';1./[1:1:n_freq]'];
cond_in=(1/size(sig_true_init,1)).*(sum(sig_true_init));%the constant conductivity 
 W=(1/cond_in).*diag(RR);%computing the initial NtD matrix, assuming
 % the conductivity value is constant in the first ring 
 W_0=W(:);%vectorize NtD matrix
 
 %collect all the conductivity values in a vector
 sig_true_ent=[];%collect all the conductivities
 for i=1:n_space-1
     sig_true_ent=[sig_true_ent;sig_true{i}];
 end
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%matrices required for Mobius propagator
 
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i=n_space-1;%initial NtD map is available on the 30th ring
[W_1]=mobius_prop(i,a,log(sig_true_ent./sig_bg),W_0,I_2,S,n_NtD,sig_bg,n_th,theta);
load AnalyticMatricesFile z L n Y D M Phi Scale
  NtD=reshape(W_1,n,n);

Lambda=NtD\(eye(size(NtD)));

% Scaling
% The scaling matrix
%Scale = diag([1./sqrt(1:n/2),1./sqrt(1:n/2)]);
%Scale_ext = [[1,zeros(1,n)];[zeros(n,1),Scale]];



% Scaling Lambda, Y and M
%Lambdasc = Scale*Lambda*Scale;
%Msc = Scale_ext*M*Scale_ext;
%Ysc = Scale_ext*Y;


% Testing the discrepancy 

Lambda_ext = [[0,zeros(1,n_NtD)];[zeros(n_NtD,1),Lambda]];
Lambda_ext = Scale*Lambda_ext*Scale;
Msc = Scale*M*Scale;
Ysc = Scale*Y;
B = Ysc*D*Phi;
C = Phi'*D*Phi;

LL = (C - B'*((Lambda_ext+Msc)\B)); 
b_mob(:,nn)=LL(:);
nn
 end
 save Mobius_data b_mob
