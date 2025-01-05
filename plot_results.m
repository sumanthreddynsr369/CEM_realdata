
clear all
clc


load highnoise_Neumann_data2 sig_app lambda_post lambda_pr_gen Vert_all n_neigh eta_cond std_noise_data %load the result

n_step=size(sig_app,2);%number of steps in the algorithm
tri=delaunay(Vert_all(1,:),Vert_all(2,:));%create a delaunay mesh, it is used in visualizing the conductivity
n_vert=size(Vert_all,2);%total number of vertices including center (used in plotting the conductivity)


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:n_step

figure(i)
%compute the number of vertices from current updated ring to the boundary
sig_app_plot=[sig_app{i};NaN(n_vert-length(sig_app{i}),1)];%the unknown conductivities set to NaN 

trisurf(tri,Vert_all(1,:),Vert_all(2,:),sig_app_plot);%plot the approximated conductivities
 view(2)
 title(['Layer number= ',num2str(i)]) 
 shading interp
axis('off')
caxis([0 5])
    colormap('jet')
 axis equal;
 colorbar;
 pause(0.05)
end