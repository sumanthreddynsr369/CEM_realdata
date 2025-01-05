%This file discribes the prior space discretization on which EnKF algorithm
%is carried out.

%%Input---
%a--- radial discretization
%f---you can control the sparsity using "f"  (0.1<f<=1)
     %AS F--->0.1 THE PRIOR DISCETIZATION BECOMES COARSER
%plot---'yes', if you want to visualize the prior space discretization

%%Output
%Verices---a cell containing array of the angular discretization points at in the
%rings which are defined by angular discretization "a"
%th_disc---a cell containing array of the angular position of the angular discretization points
%n_th----a vector containing the number of angular discretization points in
%each ring

function [Vertices,th_disc,n_th]=prior_space_discre(a,f_t,Visualize,save_plot)

n_space = length(a)-1;%number of rings

h = 1/(n_space);%angular resolution in each ring
n_th = zeros(n_space,1);
% Offset the rings by Fibonacci angle to maximize incoherence
phi = (1+sqrt(5))/2;%off set angle
th_F = 2*pi/phi^2;
th_last = 0;
r=1;
Vert_plot=[];%collect all the vertices to plot
h_n=1/(f_t*n_space);%this is used to define the number of angular discretization points in each ring
                   %used in the variable nj

for j = 1:n_space
    nj = round(2*pi*r/h_n); % Number of nodes in the jth ring
    n_th(j) = nj;%save the number of nodes in j^th ring
    thj = 2*pi/nj*[0:nj-1] + (th_last+th_F);%offset the rings by using golden ratio
    th_last = thj(1);
    Vj = [r*cos(thj);r*sin(thj)];
    th_disc{j}=thj';%save the angular location of discretized points in each ring
    %saves as a cell structure
    Vertices{j} = Vj;%save the discretization points in each ring
    %saves as a cell structure
     Vert_plot= [Vert_plot,Vj];%used to plot the discretization
     %saves all the discretization points in one array
    r = r - h;
    %f--controls the sparsity
end

%if plot is 'yes', then plot the prior discretization
 if strcmp(Visualize,'yes')
%plot the prior space discretization 
for i=1:n_space-1
    %create a tile plot on the rings, connect the angular discretization
    %points in each ring from inner boundary to the outer boundary
    V1=[a(i+1)*cos(th_disc{n_space-i}),a(i+1)*sin(th_disc{n_space-i})];%
    V2=[a(i+2)*cos(th_disc{n_space-i}),a(i+2)*sin(th_disc{n_space-i})];
    for j=1:size(V1,1)
   plot([V1(j,1) V2(j,1)],[V1(j,2) V2(j,2)],'LineWidth',1.5)%tiles
    hold on
    end
    hold on 
 center=zeros(size(a,2),2);%center for circles
 viscircles(center,a,'Color','k','LineWidth',1.5);%plot the concentric circles/rings
end
set(gca,'FontSize', 15)
title(['Prior space discretization, f_t=',num2str(f_t)])
axis equal;
   set(gca,'FontSize',20)
   hold off
   if strcmp(save_plot,'yes')%do you want to save figure to folder Latex/figs
         print -depsc ../Latex/figs/LS_space_discretization.eps %save the figure in the folder
   end
 end

 end


