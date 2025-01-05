%this function propagates the NtD particles on the boundary of current ring
%to the boundary of the unit disc where the data is available

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input: 
%i---numbering of the current ring
%a---radial discretization of the disc (from boundary to center of the disc
          %(1-->0))
%lambda_pr---prior conductivities used in propagation
%W_0---prior NtD cloud at current ring
%n_NtD--size of NtD matrix
%I---used in first order Mobius propagator
%I_2---used in second order mobius propagator
%S---used in second order mobius propagator
%sig_bg---background conductivity
%n_th--a vector containing number of angular discretization points in each
      %ring
%theta--a cell structure containing the angular discretization location of
          %all the rings
%n_sample---number of prior samples considered in the algorithm
%eta_mob---safeguard factor in estimating numerical error induced by
            %propagator

%Output:
%W_1---the propagated NtD prior cloud 
%diff_solver---componentwise difference between the numerical integrator
               %solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [R_inv,W_1,diff_solver,eig_DSDW]=mobius_prop(i,a,lambda_pr,W_0,I,I_2,S,n_NtD,sig_bg,n_th,theta,n_sample,eta_mob)

 n_ring=size(n_th,1);%total number of rings
 
                    
for j=1:n_sample
  W_old=reshape(W_0(:,j),n_NtD,n_NtD);%prior NtD sample at current ring
  W_o1=W_old;%initial particles at current ring for first order 
  %mobius solver
  W_o2=W_old;%initial particles at current ring for second order 
  %mobius solver
    a_c= [a(n_ring-i+1:end)];%Current radius of ring to the boundary 
    n_th_c=flip(n_th(1:i),1);%extract the number of angular discretization points 
                      %from current ring to the boundary
    n_iter=0;%this is used to identify the indexing of 
    %conductivity values that are used to propagate NtD matrix one step
    for l=1:length(a_c)-1
        h_r=a_c(l+1)-a_c(l);
        sig=sig_bg*exp(lambda_pr(end-n_iter-n_th_c(l)+1:end-n_iter,j));%extract the conductivity
        %values for the 1 step propagation (we do this till the boundary of unit disc)
        
    [R,D,G]=mat_build(i-l+1,sig,theta,n_NtD);%build the matrices required for 
    %layer stripping equation
    W_o1= (I*W_o1 + (h_r/a_c(l))*R)/(-(h_r/a_c(l))*(D*(G*(D*W_o1))) + I);
    %First order mobius solver
    eig_DSDW(:,l)=eigs(D*(G*(D*W_o1)),n_NtD);
    A_1=[S, (1/a_c(l))*R; (-1/a_c(l))*D*(G*D),S];%required for second order mobius solver
     Gamma= I_2 + h_r*A_1 + (h_r^2/2)*A_1^2;
     alpha=Gamma(1:n_NtD,1:n_NtD);
     beta=Gamma(1:n_NtD,n_NtD+1:2*n_NtD);
     gamma=Gamma(n_NtD+1:2*n_NtD,1:n_NtD);
     delta=Gamma(n_NtD+1:2*n_NtD,n_NtD+1:2*n_NtD);
     W_o2= (alpha*W_o2 + beta)/(gamma*W_o2 +delta);
     %Second order mobius solver 
     
     n_iter=n_iter+n_th_c(l);%used to extract conductivities required 
     %for the update of NtD matrices to next ring..
     %Notice that we do this all the way to the boundary
     
    end
    
% %     
    W_mobius_1(:,j)=W_o1(:);%collect first order propogated NtD matrices at the boundary
    
    W_mobius_2(:,j)=W_o2(:);%collect the second order propogated NtD matrices
    
    diff_solver(i,j)=norm(abs(W_mobius_2(:,j)-W_mobius_1(:,j)));%
    %compare the discrepancy in the propogation using different accuracy solvers
    
    gam=diag(eta_mob*abs(W_mobius_2(:,j)-W_mobius_1(:,j)));
    %Covariance in the predicted ensemble is found by using difference between 
    %first and second order Mobius solvers
    
    W_1(:,j)=W_mobius_1(:,j) + gam*randn(n_NtD^2,1);%Updated NtD samples 
    %along with innovation 
  load AnalyticMatricesFile z L n Y D M Phi Scale
  NtD=reshape(W_1(:,j),n_NtD,n_NtD);

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

R = (C - B'*((Lambda_ext+Msc)\B)); 
R_inv(:,j)=R(:);
end


end
