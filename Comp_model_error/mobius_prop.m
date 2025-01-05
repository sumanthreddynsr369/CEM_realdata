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



function [W_1]=mobius_prop(i,a,lambda_pr,W_0,I_2,S,n_NtD,sig_bg,n_th,theta)

 n_ring=size(n_th,1);%total number of rings
  W_o2=reshape(W_0,n_NtD,n_NtD);%initial particles at current ring for second order 
  %mobius solver
    a_c= [a(n_ring-i+1:end)];%Current radius of ring to the boundary 
    n_th_c=flip(n_th(1:i),1);%extract the number of angular discretization points 
                      %from current ring to the boundary
    n_iter=0;%this is used to identify the indexing of 
    %conductivity values that are used to propagate NtD matrix one step
    for l=1:length(a_c)-1
        h_r=a_c(l+1)-a_c(l);
        sig=sig_bg*exp(lambda_pr(end-n_iter-n_th_c(l)+1:end-n_iter));%extract the conductivity
        %values for the 1 step propagation (we do this till the boundary of unit disc)
        
    [R,D,G]=mat_build(i-l+1,sig,theta,n_NtD);%build the matrices required for 
    %layer stripping equation
    
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

W_1=W_o2(:);

end
