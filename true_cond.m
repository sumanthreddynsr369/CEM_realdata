%%%%%%%%%%%%%%%%%%%
%this file creates a cell structure of true conductivity 
      %values from the inner most ring to the outer ring
%Input: 

%V: a cell structure of discretization on disc, the vertices vertices a
      %are arranged in the order from deepest ring to the boundary ring.
%plot: 'yes' to visualize the true conductivity
%%       
%output: a cell structure of the true 
          %conductivity values at the discretization points defined by
          %"V"
%%%%%%%%%%%%%%%%%%%%
function cond_vals=true_cond(V, plot)
n_rings=size(V,2);%number of rings/radial discretization points 
                     %(excluding center)

 for i=1:n_rings
    cond_vals{i}=EvalsmoothConductivity(V{i});%construct the 
    %true smooth conductivity values at the prescribed discretization
    %points
 end
 
 
 %to plot the true conductivity distribution on prior space discretization
 
 Vert=[];%collect all the discretization points here
 cond_plot=[];%collect all the conductivity values at the 
       %discretization points
 if strcmp(plot,'yes')%if plot is yes
   for i=1:n_rings
       Vert=[V{i},Vert];%collect all the vertices to define delaunay mesh on it
   end
   cond_plot=EvalsmoothConductivity(Vert);%construct the conductivities 
       %at the corresponding space discretization
   tri=delaunay(Vert(1,:),Vert(2,:));%delaunay triangles
   trisurf(tri,Vert(1,:),Vert(2,:),cond_plot);
     
   title('True') 
view(2)
shading interp
axis('off')
    colormap('jet')
 axis equal;
  set(gca,'FontSize',18)
 colorbar;
 print -dpng ../Latex/figs/true_conductivity.png 
 end
 