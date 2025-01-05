clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%PARAMETERS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%inhomogeneous experiment number
ddd=[1,1,2,4,5,3,3];sss=[1,2,3,3,2,5,3];
for xp=1:length(ddd)
    dd=ddd(xp);ss=sss(xp);
%dd=2;ss=1;
sig_bg=1.5;         %back ground conductivity
[LL]=R_matrix_SN_inhomo(dd,ss,sig_bg);
%Global parameters:
n_space = 30;       %number of radial discretization points (excluding the center)
n_sample = 1000;    %number of NtD samples considered
L = 30;        %number of sine/cosine frequencies in the NtD matrix, hence the 
                        %size of NtD matrix is (2*n_freq)*(2*n_freq)
n_NtD=2*L;%size of NtD matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Parameters used for generting prior NtD sample (based on radially
%symmetric conductivites):
eta=0.2;   %variance in generative Gaussian model for radially symmetric conductivity sample
alp=2;     %this controls the correlation between the adjacent conductivity values
view_radial='no';%visualize the radially symmetric conductivity samples based on which prior 
                      %NtD sample is generated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Parameters used for generating prior space discretization on which the
%EnKF algorithm is performed    
  f_t=1;%sparsity controller for prior space discretization
          %as f--->0.1, the prior space discretization becomes coarser f\in[0.1,1].   
  view_disc='no';%visualize the prior space discretization yes/no
  save_plot='no';%do you want to save the figure? yes/no
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 %Parameters used for generating the conductivity sample in the first ring:  
eta_first=0.004;%variance in the Gaussian generative model for the first ring
alph_first=0.01;%correlation relaxation parameter
view_pr_1='yes';%visualize random five draws from the generated sample for the first ring
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 %Parameters used for loading and generating noisy data:
 data=1; %data=1 for loading FEM data and data=2 for loading data computed using forward mobius solver
  %noise_data_level=10^(-7); %the relative noise in the data for additive noise =0.03%
 
 %noise_data_level=3.6*10^(-7);%0.06%
 noise_data_level=10^(-6);%0.1%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 %Parameters used for generating conditional prior:
 n_neigh=4;%number of neighbors to be considered in creating graph Laplacian
  view_neigh='no';%visualize the neighbors (scatter plot)
 eta_cond=0.009;%variance in the generative Gaussian model for conditional smoothness prior
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 %Parameter used in propagation step:
  eta_mob=0.5;%safeguard factor in estimating numerical error induced by
            %propagator
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Radial discretization: 
a_in=0;             %center of the disc
a_fin=1;            %radius of the disc
a=[a_in:(a_fin-a_in)/n_space:a_fin];%radial discretization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Generating prior space discretization:  
[Vertices,theta,n_th]=prior_space_discre(a,f_t,view_disc,save_plot);%this program generates the 
            %prior space discretization details
    center=[0;0];%add origin to the vertices list
%collect all the vertices/prior discretization points 
Vert_all=[];
for i=1:size(Vertices,2)
Vert_all=[Vert_all,Vertices{i}];
end  
Vert_all=[Vert_all,center];%all the prior space discretization points
dist=sqrt(Vert_all(1,:).^2+Vert_all(2,:).^2);%distances of vertices from the center
tri=delaunay(Vert_all(1,:),Vert_all(2,:));%create a delaunay mesh, it is used in visualizing the conductivity
n_vert=size(Vert_all,2);%total number of vertices including center (used in plotting the conductivity)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Construct a test conductivity distribution:
view_true_cond='no';
 sig_true=true_cond(Vertices, view_true_cond);%this file creates a cell structure of true conductivity 
      %values from the outer most ring to the inner most ring
  sig_true_first=sig_true{1};%extract true conductivity values in the first ring to evaluate
                                %the initialization of the algorithm
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generating prior NtD particles:              
[~,Sig_samp_radial]=cond_sample_radial(n_space,n_sample,eta,alp,sig_bg,view_radial);
%this program generates the sample of 
              %smooth radially symmetric conductivity samples based on Aristotelian boundary condition
[NtD_samp]=NtDmatrix_sample_radial(L,Sig_samp_radial,a);%this program generates the prior NtD sample
               %based on radially symmetric smooth conductivities generated in the
               %previous step (cond_sample_radial) at the boundaries of all
               %the rings. Size---(L*L)\times(L) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load the data:
%load the data computed using FEM if data=1
%the data has been generated by using NtD_comp_smooth_NS2.m under the
      %folder generate data.
 if data==1 
  load Real_data_matrix_homo LL_diff R_CEM
  %R_inv_CEM=R_CEM\eye(size(R_CEM));
  %LL_CEM=R_inv_CEM(:);
  load_data_file = sprintf('Real_data_matrix_inhomo_%d_%d', dd, ss)  ;
  load (  load_data_file, 'LL') 
  NtD_data=LL;
    %the size of NtD_data matrix is (2*L)*(2*L)
  load Model_err_FEMcoarsevsMobius err_mean err_cov
 else
%load the data computed using forward mobius propagator if data=2
%the data has been generated by using mobiusprop_data.m under the
      %folder generate data_FEM.
  load Input_data_mobius NtD_mobius
  NtD_data=NtD_mobius;
 end
eta_data=(noise_data_level*(max(max(abs((NtD_data))))^2));%std^2 of Gaussian noise added to each component of NtD data matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generating prior conductivities for the first ring:
samp_first=ring1_sample(Vertices,n_sample,eta_first,alph_first,view_pr_1);%this function is used to produce a smooth 1d sample of 
    %conductivitivity profiles for the first ring (outermost ring)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Matrices needed for drawing Gaussian conditional smoothness prior:
   L2=graph_laplace(Vertices,Vert_all,n_neigh,view_neigh);%this file creates the 
   %graph Laplacian matrix
   S_cond=L2'*L2;%Precision matrix for generative Gaussian model
   n_pr=0;%used to plot the generated conditional prior samples
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matrices required for Mobius propagator:
 I=eye(n_NtD);%will be used in first order mobius update
  I_2=eye(2*n_NtD);%will be used in second order mobius update
  S=zeros(n_NtD);%will be used in second order mobius update
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matrices required for EnKF update:
n_data=size(LL,1);
B1=speye(n_data^2);%used in determining observation matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Generating the noisy data:
 D_add=LL_diff(:)*LL_diff(:)' + eta_data .*eye(n_data^2);%observation model covariance matrix for additive noise
 %adding noise to the data
 ll=size(D_add,1);

   b_noisy=NtD_data ;%noisy data 
   %b_noisy= (b_noisy + b_noisy')/2;%symmetrize the data after adding the noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   %%Generating ensemble of observations:
   D= n_space.*(D_add+err_cov);%Covariance matrix for noise in the data additive+ modeling error
   %and adjusted to compensate the re-usage of the data
    R1=chol(D);%R1'*R1=D
  for k=1:n_sample    
     obs_ens =-err_mean - LL_diff(:) + (b_noisy(:)) + (R1'*randn(n_data^2,1));%generate 
%     ensemble of observation particles
     Obs_ens=reshape(obs_ens,n_data,n_data);
     %Obs_ens=(Obs_ens+Obs_ens')/2;%symmetrize the observation ensemble
     b(:,k)=Obs_ens(:);%collect all the observation ensemble
  end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
   
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   
   
   
 %EnKF algorithm:
 
 lambda_pr=[];%prior conductivity distribution
 lambda_post=[];%posterior conductivity distribution
 for i= 1:n_space-1%run the Kalman update from first layer till the last layer 
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %step 1: Generate the prior NtD cloud of same size as data
     Y_i=(reshape(NtD_samp(:,n_space-i),n_NtD/2,n_sample));
   %extract prior particle cloud of NtD matrices 
   %at internal boundary of current ring (the size is L*1 (diagonal))
   %we extrapolate it by copying the same diagonal entries which makes the prior NtD
   %sample of size (2*L)*(2*L)
          for k=1:n_sample
               W_R=diag([Y_i(:,k);Y_i(:,k)]);
               W_i(:,k)=W_R(:);
          end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %step 2: Generate the prior conductivities
         if i==1
            lambda_pr_gen{i}=samp_first;% generated prior sample for the first ring
            lambda_pr{i}=lambda_pr_gen{i};
         else
             
             [lambda_pr_gen{i},lambda_pr{i}]=conditional_pr(i,n_th,lambda_post{i-1},S_cond,n_sample,eta_cond);
             %generated prior sample for the current ring        
         end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %step 3: Propagate the NtD prior cloud to the boundary/prediction step
         [R_inv,W_0,diff_solver]=mobius_prop(i,a,lambda_pr{i},W_i,I,I_2,S,n_NtD,sig_bg,n_th,theta,n_sample,eta_mob);
         %this file gives the propagated particles W_1
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %step 4: Analysis step
      
       Z=[R_inv;(lambda_pr{i})];%augmented predicted ensemble concatenated with conductivity values from current ring till the
                        %boundary.   
   Z_mean=(1/n_sample).*sum(Z,2);%mean of the predicted ensemble 
   Z_c= Z - (Z_mean*ones(1,n_sample));%centralizing the predicted ensemble
   C_0=(1/(n_sample-1)).*(Z_c*Z_c'); %covariance of augmented prediction ensemble
 %taking care of rank deficiency of covariance matrix

         if cond(C_0)>10^(5)
           ep=eigs(C_0,1)*10^(-5);
           C_0= C_0 + ep.*eye(size(C_0));%to make sure covariance matrix is positive definite and
                                         %well-conditioned
         end
    B=[B1,zeros(n_data^2,size(lambda_pr{i},1))];%Observation matrix
    %Notice that its size changes according to the size of parameter    
  %Kalman gain matrix
   K=(C_0*B')/((B*C_0*B')+D); 
 %generate the posterior states
   Z_1 = Z + K*(b - B*Z);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
  lambda_post{i}=Z_1((n_data^2)+1:end,1:end);%extract the updated posterior conductivities  
  
  sig_app{i}=sig_bg*exp((1/n_sample)*sum(lambda_post{i},2)); %average of the
                         % posterior conductivity samples 
  
       sig_post_1=sig_app{i};%extract the posterior update in the first ring from cell structure                  
      sig_app_first=sig_post_1(1:n_th(1),:);
      %average of the samples of log conductivities in the first ring


%plot the true and approximated conductivities in the first ring
figure(5)
plot(sig_app_first,'LineWidth',2)
hold on

set(gca,'FontSize',15)
title('Angular conductivity values in the first ring')
legend('approximate','Location','southeast')
hold off
pause(0.05)
hold on
                
%plot the approximated posterior conductivity distribution for each ring
%update
figure(5+i)
%compute the number of vertices from current updated ring to the boundary
sig_app_plot=[sig_app{i};NaN(n_vert-length(sig_app{i}),1)];%the unknown conductivities set to NaN 

trisurf(tri,Vert_all(1,:),Vert_all(2,:),sig_app_plot);%plot the approximated conductivities
 view(2)
 caxis([0 5])
 title(['Layer number= ',num2str(i)]) 
 shading interp
axis('off')
    colormap('jet')
 axis equal;
 colorbar;
 
 
if i>1%except the first layer
 n_pr=n_pr+n_th(i-1);%length of approximated/posterior conductivity values      
 ind_post=randi(1000,1,9);%draw 9 generated random prior conditional samples
 tri_pr=delaunay(Vert_all(1,n_pr+1:end),Vert_all(2,n_pr+1:end));%create a mesh on vertices 
                                                      %those correspond to
                                                    %drawn conditional prior conductivity values
    figure(40+i) %plot the conditional prior samples  
    lambda_pr_plot=lambda_pr_gen{i};%extract the generated prior samples from cell structure
    
    %plot the generated prior conditional samples
  for ll=1:9 
      ind_plot=ind_post(ll);
   subplot(3,3,ll)
  trisurf(tri_pr,Vert_all(1,n_pr+1:end),Vert_all(2,n_pr+1:end),sig_bg*exp(lambda_pr_plot(:,ind_plot)));%plot the approximated conductivities
 view(2) 
 shading interp
axis('off')
    colormap('jet')
 axis equal;
 colorbar;
 set(gca,'FontSize', 15)
  end
  suptitle(['Layer number= ',num2str(i)])
end
i
 end
 filename = sprintf('main1_result_inhomo_%d_%d.mat', dd,ss);
save(filename, 'sig_app', 'lambda_post', 'lambda_pr_gen', 'Vert_all', 'n_neigh', 'eta_cond', 'eta_data')

close all;
[c_max]=plot_results2(dd,ss);%plot the final result

end

 