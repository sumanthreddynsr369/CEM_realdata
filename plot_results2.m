%the results are crated by using the file Enkf_smooth.m


function [c_max]=plot_results2(dd,ss)
  filename = sprintf('main1_result_inhomo_%d_%d.mat', dd,ss);
  load (filename)      %load the result

%pr_update--- contains the updated posterior mean conductivities 
%pr_update--size 
%(number of prior space discretizaton points*number of layers)
%Hence each column contains the updated conductiviy values for each layer
%The first layer update is the last column, second layer update is 
%the second column from the last and the orering follows
% at the final ring we dont have any prior to update hence we can use the
% average of conductivity values in the second last updated conductivity
% values

%pr_update--- contains the updated posterior mean conductivities 
%pr_update--size 
%(number of prior space discretizaton points*number of layers)
%Hence each column contains the updated conductiviy values for each layer
%The first layer update is the last column, second layer update is 
%the second column from the last and the orering follows
% at the final ring we dont have any prior to update hence we can use the
% average of conductivity values in the second last updated conductivity
% values

c_max=max(sig_app{29});
tri=delaunay(Vert_all(1,:),Vert_all(2,:));
n_vert=size(Vert_all,2);%total number of vertices including center (used in plotting the conductivity)
sig_bg=1.5;
n_elem = size(tri,1);
% Plotting an nnxmm array of reconstructions
gap = 0.2; % Gap between the images
nn = 2;
mm = 3;
% Image centers
cx = (2+gap)*(0:mm-1);
cy = (2+gap)*(nn-1:-1:0);
c  = [kron(ones(1,nn),cx);kron(cy,ones(1,mm))];

% Selecting the images to be plotted

%ii = [30,27,24,21,18,15,12,9,6,3,2,1];

ii=[2,10,15,20,24,29];


for ind = 1:length(ii)
    i = ii(ind);
    sig=[sig_app{i};NaN(n_vert-length(sig_app{i}),1)];


h=figure(1);
 for j = 1:n_elem
       xj = [Vert_all(1,tri(j,1))+c(1,ind),Vert_all(1,tri(j,2))+c(1,ind), ...
                Vert_all(1,tri(j,3))+c(1,ind),Vert_all(1,tri(j,1))+c(1,ind)];
       yj = [Vert_all(2,tri(j,1))+c(2,ind),Vert_all(2,tri(j,2))+c(2,ind), ...
                Vert_all(2,tri(j,3))+c(2,ind),Vert_all(2,tri(j,1))+c(2,ind)];
       uj = [sig(tri(j,1)),sig(tri(j,2)),sig(tri(j,3)),sig(tri(j,1))];
       hpp= patch(xj,yj,uj);
       set(hpp,'EdgeColor','interp');
       hold on
 end

 text(c(1,ind)-1,c(2,ind)+1,num2str(i),'FontSize',15)
% axis('square')
% colorbar('location','SouthOutside')
% set(gca,'FontSize',20)
% axis('off')
 
 
end



tri=delaunay(Vert_all(1,:),Vert_all(2,:));%create delaunay mesh
n_elem = size(tri,1);


 
caxis([0 4])
%shading interp
axis('off')
colormap('jet')
axis('equal')
colorbar
set(gca,'FontSize',15)
axis('off')
hold off
suptitle(['\eta_c=',num2str(eta_cond),', k=',num2str(n_neigh)]);
set(gca,'FontSize',15)
saveas(h,sprintf('main1_result_%d_%d.png',dd,ss));
end

%plot the generated prior samples




