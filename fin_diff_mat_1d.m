%This function gives the space discretization for
%prior radially symmetric conductivity sample
% and gives the smoothness promoting prior second order finite difference
% matrix
%with Aristotelian boundary condition (p.g.no 154 in Ten lectures on
%subjective computing by DC and ES)
%
%%input
%a_in- initial point
%a_fin- final point
%n_num- number of discretization points
%%output
%a-gives the discretization
%L-second order finite difference matrix

function [a,L]=fin_diff_mat_1d(a_in,a_fin,n_num)
a=[a_in:(a_fin-a_in)/n_num:a_fin];%radial discretization
%create the 2nd order finite difference matrix
aux = zeros(n_num+1,1);
aux(1) = 2;
aux(2) = -1;
L_D  =  toeplitz(aux);% second order finite difference
% matrix
 L_Dinv=L_D\eye(size(L_D));
 L_Dinv_T=L_Dinv';
del= 1./sqrt(sum(L_Dinv_T(:,floor(n_num/2)).^2));%augment the smooth ness matrix 
%using  Aristotelian boundary condition
L=0.5.*[2*del,zeros(1,n_num);L_D(2:end-1,1:end);zeros(1,n_num),2*del];
end
