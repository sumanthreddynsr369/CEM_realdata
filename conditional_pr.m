
%this function is used to generate the Gaussian condictional smoothness
%prior samples
%%%%%%%%%%%%%%%%%%%
%Input:
%i---numbering of the current ring
%Vert_all--all the prior discretization points including the center
%B---generative Gaussian precision matrix
%n_Sample---number of conditional samples to be generated
%eta_cond---variance in the generative Gaussian model

%Output:
%lambda_pr_gen---the generated prior conductivity sample inside (entire region)
                 %the current boundary of the ring
 %lambda_pr---prior conductivity sample for the next update where we attach
              %the generated current ring conductivity values to the previous update 
                    

function [lambda_pr_gen,lambda_pr]=conditional_pr(i,n_th,lambda_post,B,n_sample,eta_cond)
   

%the B matrix is already ordered according to known and unknown conductivity
                                     %indices
                                     
  n_known=0;%index of known conductivity values upon which we are going to condition
  for j=1:i-1
      n_known=n_known+n_th(j);
  end
  B11=B(n_known+1:end,n_known+1:end);%unknown block
  B12=B(n_known+1:end,1:n_known);%known cross-block
  
  n_c=size(B11,1);%size of unknown conductivities
  R=chol(B11);%cholesky factor of the unknown block 
  
  lambda_pr_gen=B11\(-B12*lambda_post) + eta_cond*(R\randn(n_c,n_sample));%generate the prior conditioned
  %realizations in the space inside the current ring (entire region)
  lambda_pr_ring=lambda_pr_gen(1:n_th(i),:);%extract the portion of conductivities corresponding
  %to current ring
  lambda_pr=[lambda_post;lambda_pr_ring];%attach the generated prior samples to the current posterior