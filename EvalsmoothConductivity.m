function CondVals = EvalsmoothConductivity(xi)

% This function program evaluates the generated conductivity at required
% points.
% Input: xi  - 2xn array, coordinates at which the value is needed
%        CondVertices - array of vertices where the values are given
%        CondTopol    - topology matrix
%        sigma        - conductivity values at vertices

% Define the conductivity
sigma_bg = 1.5; % 1

center1 = [0.5;-0.1];%middle
width1 = 0.2;
ampl1 = 2.0; 

center2 = [-0.4;-0.7]; %bottom blob
width2 = 0.15;
ampl2 = -1.2; % 0.6;

center3 = [-0.4;0.6]; %[0.4,-0.4];top
width3 = 0.2;
ampl3 =  3.0; % 0.6;

aux = (xi(1,:) - center1(1)).^2 + (xi(2,:) - center1(2)).^2;
blob1 = ampl1*exp(-1/(2*width1^2)*aux)';

aux = (xi(1,:) - center2(1)).^2 + (xi(2,:) - center2(2)).^2;
blob2 = ampl2*exp(-1/(2*width2^2)*aux)';

aux = (xi(1,:) - center3(1)).^2 + (xi(2,:) - center3(2)).^2;
blob3 = ampl3*exp(-1/(2*width3^2)*aux)';

CondVals = sigma_bg + blob1 + blob2 + blob3;



