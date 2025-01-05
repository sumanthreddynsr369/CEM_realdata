% Testing the performance of the inverse solver in the case of
% smooth conductivity. The DtN matrix is not analytically
% known. The performance is tested with different truncation levels 
% of the CGLS solver. This is using the scaled matrices
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
%clear all;
clc;
% Selections: 
% (a) Centered vs. non-centered sample
%centered = 'yes';
centered = 'no';
% (b) regularization: CGLS or truncated PCA
%method = 'CGLS';
method = 'tPCA';
% (c) conductivity
cond_model = 'constant';
%cond_model = 'smooth';
%cond_model = 'jump'; % Not yet available
noiselevel = 0.001;
%noiselevel = 0;
%--------------------------------------------------------------------------

% Fine mesh to generate data
%load VoltageMesh_15 Vertices Topol ElectrodeNodes
% Coarser mesh
load VoltageMesh_12 Vertices Topol ElectrodeNodes_ext V_extra Topol_2
load Real_data_matrix_homo z

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

 S11 = sparse(n_nodes_2,n_nodes_2);%Stiffness matrix
% 
% 
 xw = TriGaussPoints(6);
 nodes = xw(:,1:2);
 weights = xw(:,3);

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
          if strcmp(cond_model,'smooth')
              CondVals=EvalConductivity(Z);
          elseif strcmp(cond_model,'constant')
              CondVals=1.5*ones(size(weights));
          end
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
       EdgeLength=sqrt(sum((V_extra(:,Ij(2)) - V_extra(:,Ij(1))).^2,1))+...
                   sqrt(sum((V_extra(:,Ij(3)) - V_extra(:,Ij(2))).^2,1));
        %EdgeLength=sqrt(sum((V_extra(:,Ij(3)) - V_extra(:,Ij(1))).^2,1));
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
       
      ints(2*(j-1)+1)= ints(2*(j-1)+1)+ EdgeLength/6;
        
           
      ints(2*(j-1)+2)=(2*EdgeLength)/3;
  
      ints(2*(j-1)+3)=  ints(2*(j-1)+3)+ EdgeLength/6;
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
 R = (LL)\eye(size(LL));
 

% Compute the DtN numerically

S011 = S0(I_bdry_2,I_bdry_2);
S012 = S0(I_bdry_2,I_int);
S021 = S0(I_int,I_bdry_2);
S022 = S0(I_int,I_int);
DtN = S011 - S012*(S022\S021); % DtN in roodtop basis

% Transformation into the cosine-sine basis

load AnalyticMatricesFile z L n Y D M Phi Scale
x_bdry = V_extra(1,I_bdry_2);
y_bdry = V_extra(2,I_bdry_2);
th = atan2(y_bdry,x_bdry)';
G = 1/sqrt(pi)*[cos(th*(1:n/2)),sin(th*(1:n/2))];

Lambda = G'*DtN*G;
NtD=Lambda\(eye(size(Lambda)));

% Scaling
% The scaling matrix
%Scale = diag([1./sqrt(1:n/2),1./sqrt(1:n/2)]);
%Scale_ext = [[1,zeros(1,n)];[zeros(n,1),Scale]];



% Scaling Lambda, Y and M
%Lambdasc = Scale*Lambda*Scale;
%Msc = Scale_ext*M*Scale_ext;
%Ysc = Scale_ext*Y;


% Testing the discrepancy 

Lambda_ext = [[0,zeros(1,n)];[zeros(n,1),diag(diag(Lambda))]];
Lambda_ext = Scale*Lambda_ext*Scale;
Msc = Scale*M*Scale;
Ysc = Scale*Y;
B = Ysc*D*Phi;
C = Phi'*D*Phi;

ModError = (C - B'*((Lambda_ext +Msc)\B)) - LL; 
RelModError = norm(ModError,'fro')/norm(LL,'fro');
ModErrorNorm = norm(ModError(:));
imagesc(log(abs(ModError)))
colorbar;
title(max(max(abs(ModError))));
%save data_32_electrode LL Msc B C 
