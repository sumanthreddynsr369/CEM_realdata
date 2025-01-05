% Compute analytically the matrices appearing in the matrix equation
% These matrices are not scaled by 1/sqrt(j), but the scaling matrix is
% computed. For proper scaling, one needs to multiply
%
% Y --> Scale*Y
% M --> Scale*M*Scale
%--------------------------------------------------------------------------
%
%
%--------------------------------------------------------------------------
%clear all;
clc

load VoltageMesh_12 Vertices Topol ElectrodeNodes_ext V_extra Topol_2
load Real_data_matrix_homo z
L=16;
z=0.0017*ones(1,L);
% The mesh is needed only for extracting the electrode geometry, the
% meshing is irrelevant
L=size(ElectrodeNodes_ext,1);


n = 60; % n/2 cosines, n/2 sines, one constant basis vector

% Finding the angles that define the electrodes
[L,kL] = size(ElectrodeNodes_ext); % L = number of electrodes

ElectrodeAngles = NaN(2,L);
for ell = 1:L
    jmin = ElectrodeNodes_ext(ell,1);
    jmax = ElectrodeNodes_ext(ell,kL);
    vmin = V_extra(:,jmin);
    vmax = V_extra(:,jmax);
    thmin = atan2(vmin(2),vmin(1));
    thmax = atan2(vmax(2),vmax(1));
    ElectrodeAngles(:,ell) = [thmin,thmax];
end

%Checking
for ell = 1:L
   thmin = ElectrodeAngles(1,ell);
   thmax = ElectrodeAngles(2,ell);
   fill([0,cos(thmin),cos(thmax)],[0,sin(thmin),sin(thmax)],'b')
   axis([-1,1,-1,1])
   axis('square')
   hold on
   pause(0.5)
end
hold off

% 1. The matrix Y

Y = NaN(n+1,L);
width = NaN(1,L);
for ell = 1:L
    thmin = ElectrodeAngles(1,ell);
    thmax = ElectrodeAngles(2,ell);
    if thmax>thmin
       width(ell) = thmax - thmin;
    else
      % the angles are at different sides of the negative x-axis
      width(ell) = (thmax + 2*pi) - thmin;
    end
    Y(1,ell) = 1/sqrt(2*pi);
    for j = 1:n/2
        Y(j+1,ell) = (1/(sqrt(pi)*j*width(ell)))*(sin(j*thmax)-sin(j*thmin));
        Y(n/2+j+1,ell) = (1/(sqrt(pi)*j*width(ell)))*(cos(j*thmin)-cos(j*thmax));
    end
end

% 2. The matrix D

D = diag(width./z);

% 3. The matrix M

M = NaN(n+1,n+1);
M(1,1) = (1/(2*pi))*trace(D);
thmin = ElectrodeAngles(1,:);
thmax = ElectrodeAngles(2,:);
for k = 1:n/2
    M(1,k+1) = (1/(sqrt(2)*pi*k))*sum((sin(k*thmax) - sin(k*thmin))./z);
    M(1,n/2+k+1) = (1/(sqrt(2)*pi*k))*sum((cos(k*thmin) - cos(k*thmax))./z);
    M(k+1,1) = M(1,k+1);
    M(n/2+k+1,1) = M(1,n/2+k+1);
end

for k = 1:n/2
    for j = 1:n/2
        aux1 = 1/(j+k)*sum((sin((j+k)*thmax) - sin((j+k)*thmin))./z);
        if j ~= k
            aux2 = 1/(j-k)*sum((sin((j-k)*thmax) - sin((j-k)*thmin))./z);
        else
            aux2 = sum(width./z);
        end
        M(1+j,1+k) = (1/(2*pi))*(aux1 + aux2);
        M(1+n/2+j,1+n/2+k) = (1/(2*pi))*(aux2 - aux1);
        aux3 = 1/(j+k)*sum((cos((j+k)*thmin) - cos((j+k)*thmax))./z);
        if j ~= k
            aux4 = 1/(j-k)*sum((cos((j-k)*thmin) - cos((j-k)*thmax))./z);
        else
            aux4 = 0;
        end
        M(1+j,n/2+1+k) = (1/(2*pi))*(aux3 + aux4);
        M(n/2+1+k,1+j) = (1/(2*pi))*(aux3 + aux4);
    end
end

% 4. The matrix Phi

cc = cos(2*pi/L*(0:L-1)'*(1:L/2));
ss = sin(2*pi/L*(0:L-1)'*(1:L/2-1));
Phi = [cc,ss];
% Normalizing
normPhi2 = L/2*ones(1,L-1);
normPhi2(L/2) = L;
Phi = Phi*diag(1./sqrt(normPhi2));

Scale = diag([1,1./sqrt(1:n/2),1./sqrt(1:n/2)]);

save AnalyticMatricesFile z L n Y D M Phi Scale







