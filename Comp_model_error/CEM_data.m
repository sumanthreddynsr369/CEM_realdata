
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
clear all;
clc;
% Fine mesh to generate data
%load VoltageMesh_15 Vertices Topol ElectrodeNodes
% Coarser mesh
load VoltageMesh_12 Vertices Topol ElectrodeNodes_ext V_extra Topol_2
load Real_data_matrix_homo z
n_samp=250;
n_nodes_2=size(V_extra,2);%total number of nodes 
n_elems = size(Topol_2,1);%number of elements
% Boundary nodes
radii2 = sum(V_extra.^2,1);
I_bdry_2 = find(radii2 >1 - 1e-6);%Boundary nodes
I_int  = setdiff(1:n_nodes_2,I_bdry_2);%interior nodes

[L,p] = size(ElectrodeNodes_ext);
p = (p-1)/2; % Number of elements under each electrode



% Contact impedances and conductivity

%conductivity = @(z) 1.5*ones(size(z,1),1); % Constant conductivity for test purposes

% Basis for the voltage patterns. Here the trig pattern. Observe: the sum
% along columns is zero (Kirchhoff's law)

freq = 2*pi/L*(1:L/2);
el_ind = (0:L-1)';
Volt = [cos(el_ind*freq),sin(el_ind*freq(1:L/2-1))];

% Normalize the basis of electrode voltages and currents

nVolt    = sqrt(L/2)*ones(1,15);
nVolt(8) = sqrt(L);
Volt     = Volt*diag(1./nVolt);

% Voltage potential block
% Compute the stiffness matrix

A=[p_eval(0,0);p_eval(1,0);p_eval(0,1);p_eval(1/2,0);p_eval(1/2,1/2);p_eval(0,1/2)];
%second order shape functions at the 6 nodes

Coeff=A\eye(size(A));%Coefficients for shape functions

 
% 
 xw = TriGaussPoints(6);
 nodes = xw(:,1:2);
 weights = xw(:,3);
for ss=1:n_samp
    S11 = sparse(n_nodes_2,n_nodes_2);%Stiffness matrix
% 
 for ell = 1:n_elems
   %step 1:
    I = Topol_2(ell,:);
     V = V_extra(:,I);  % Vertices
     
     
     bdry_elemns=intersect(I,I_bdry_2);
     
          Kloc=zeros(size(V,2));
          
    if length(bdry_elemns)>2%if the elements are on the boundary
       Z= secondorder_map(nodes',V,Coeff);%Z=F_b(t,s) the isoparametric map into physical domain
     else
         %if the elements are in the inner domain
         F = [V(:,2) - V(:,1), V(:,3) - V(:,1)];
         Z = V(:,1) + F*nodes';
    end
          
             CondVals=EvalConductivity(Z,Vertices,ss);
          
       for k=1:size(nodes,1)%over the quadrature points
       %step 2:    
        if  length(bdry_elemns)>2 
            
          J=V*Dphi(nodes(k,:));%Jacobian if element is on the boundary
        else
            J=F;%if the element is in the internal domain
        end
          Grad_phi=Dphi(nodes(k,:));%gradient of shapefunctions
          %step 3:
          G= (J')\(Grad_phi');%required in the integral computation
          B=G'*G; %step 4
          d=abs(det(J));
          
          Kloc= Kloc + 0.5*weights(k)*CondVals(k)*d*B;  %Compute the local stiffness matrix and 
          %loop over all the quadrature nodes
       end
    S11(I,I)= S11(I,I) + Kloc;

 end
% Save this matrix for DtN computation

S0 = S11;

% Add te contribution from the electrodes

for ell = 1:L
    % One electrode at a time
    I = ElectrodeNodes_ext(ell,:);
   % Observe that the nodes
    % under electrode are ordered in counterclockwise order
    
    for j = 1:p
        Ij = I(2*(j-1) + 1:2*j + 1);
       
        EdgeLength=sqrt(sum((V_extra(:,Ij(3)) - V_extra(:,Ij(1))).^2,1));
        Kloc = (1/z(ell))*(EdgeLength/30) *[4 2 -1;2 16 2;-1 2 4];
        S11(Ij,Ij) = S11(Ij,Ij) + Kloc;
    end
end

% Save the electrode lengths
ElLengths = zeros(1,L);

S12 = sparse(n_nodes_2,L-1);
for ell = 1:L
    % One electrode at a time
    I = ElectrodeNodes_ext(ell,:);
    v_el = V_extra(:,I);
    % Lengths of the p edges under the electrode. Observe that the nodes
    % under electrode are ordered in counterclockwise order
    EdgeLengths = sqrt(sum((v_el(:,2:(2*p)+1) - v_el(:,1:2*p)).^2,1));
    ElLengths(ell) = sum(EdgeLengths);
    % Integrals of basis functions under the electrode
    ints=zeros(length(I),1);
       for j = 1:p
        Ij = I(2*(j-1) + 1:2*j + 1);
       
        EdgeLength=sqrt(sum((V_extra(:,Ij(3)) - V_extra(:,Ij(1))).^2,1));
       
      ints(2*(j-1)+1)= ints(2*(j-1)+1)+EdgeLength/6;
        
           
      ints(2*(j-1)+2)=(2*EdgeLength)/3;
  
      ints(2*(j-1)+3)=  ints(2*(j-1)+3)+EdgeLength/6;
       end
      
    Ints = (ints');
    Kloc = -(1/z(ell))*Ints'*Volt(ell,:);
    S12(I,1:L-1) = S12(I,1:L-1) + Kloc;
end

% Electrode block

S22 = Volt'*diag(ElLengths./z)*Volt;
S = [[S11,S12];[S12',S22]];

% The inverse of the resistance matrix (call it conductance matrix)
% Is the Schur complement S/S22
% this is the data
LL = (S22 - S12'*(S11\S12));
b_CEM(:,ss)=LL(:);
ss
end
save CEM_data b_CEM

