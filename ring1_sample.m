
%this function is used to produce a smooth 1d sample of 
    %conductivitivity profiles for the first ring (outermost ring)
    
%%%%%%%%%%%%
%Input: 
%Vert--- the cell structure of discretization points on the 
       %entire disk
%n_samp---- the number of samples to be considered
%gam_first---the variance in the stochastic approximation
%alp_first---the correlation relaxation parameter
%%%%%%%%%%%%
%Output:
%samp_first----sample of smooth one dimensional 
          %conductivities generated for the first ring (outermost ring)

function samp_first=ring1_sample(Vert,n_samp,eta_first,alp_first,visualize)
V=Vert{1};%extract the vertices in the first ring (outer most ring)
n_num=size(V,2);%number of angular discretization/vertices in the first ring
%create the augmented second order finite difference matrix
aux = zeros(n_num,1);
aux(1) = 2;%main diagonal entries
aux(2) = -1;%off diagonal entries
L_D  =  toeplitz(aux);%second order finite difference matrix

%generate the specified number of smooth one dimensional conductivity 
 %profiles (function of angle only) for the first ring 
 for l=1:n_samp
 samp_first(:,l)=(L_D+alp_first*eye(size(L_D)))\(eta_first.*randn(n_num,1));
 end
 %plot five random draws from the generated sample
 r_draw=randi([1 1000],5,1);
if strcmp(visualize,'yes')
    figure(4)
     for i=1:5
         t=r_draw(i);
   plot(samp_first(:,t),'LineWidth',2);
   xlabel('Angular discretization index')
   ylabel('log(\sigma/\sigma_b)')
   title('Samples for the first ring (\lambda_1)')
   set(gca,'FontSize',18)
%   
    hold on
     end
        
end
