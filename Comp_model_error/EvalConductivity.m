
%This program generates a random conductivity sample and evaluates the
%conductivity value at enquiry xi point

%input:
%xi---enquire point
%Vertices---the discretization of vertices in unit disc
%i---sample number, used to fix the random seeding

%output:
%CondVals---conductivity value of ith sample at xi
function CondVals = EvalConductivity(xi,Vertices,i)
sigma_bg=1.5;
n_total=size(Vertices,2);
rng(i);
n=randi(4);%number of centers
rng(i);
c_indices=randi(n_total,n,1);
centers=Vertices(:,c_indices);

if n==1
center1 = centers(:,1);
rng(i)%fix the random number generator
width1 = 0.1 + 0.2*rand ;%width of the inclusion
rng(i+1)
ampl1 = 3.2*rand - 0.5; %amplitude of the inclusion
aux = (xi(1,:) - center1(1)).^2 + (xi(2,:) - center1(2)).^2;
blob1 = ampl1*exp(-1/(2*width1^2)*aux)';
blob2=zeros(size(blob1));
blob3=zeros(size(blob1));
blob4=zeros(size(blob1));
elseif n==2
   
center1 = centers(:,1);
 rng(i)
width1 = 0.1 + 0.2*rand ;
rng(i+1)
ampl1 = 3.2*rand - 0.5; 
aux = (xi(1,:) - center1(1)).^2 + (xi(2,:) - center1(2)).^2;
blob1 = ampl1*exp(-1/(2*width1^2)*aux)';
center2 = centers(:,2);
rng(i+2)
width2 = 0.1 + 0.2*rand ;
rng(i+3)
ampl2 = 3.2*rand - 0.5; 
aux = (xi(1,:) - center2(1)).^2 + (xi(2,:) - center2(2)).^2;
blob2 = ampl2*exp(-1/(2*width2^2)*aux)';
blob3=zeros(size(blob1));
blob4=zeros(size(blob1));
elseif n==3
    center1 = centers(:,1);
    rng(i)
width1 = 0.1 + 0.2*rand ;
rng(i+1)
ampl1 = 3.2*rand - 0.5; 

center2 = centers(:,2);
rng(i+2)
width2 = 0.1 + 0.2*rand ;
rng(i+3)
ampl2 = 3.2*rand - 0.5; 

center3 = centers(:,3);
rng(i+4)
width3 = 0.1 + 0.2*rand ;
rng(i+5)
ampl3 = 3.2*rand - 0.5; 

aux = (xi(1,:) - center1(1)).^2 + (xi(2,:) - center1(2)).^2;
blob1 = ampl1*exp(-1/(2*width1^2)*aux)';

aux = (xi(1,:) - center2(1)).^2 + (xi(2,:) - center2(2)).^2;
blob2 = ampl2*exp(-1/(2*width2^2)*aux)';

aux = (xi(1,:) - center3(1)).^2 + (xi(2,:) - center3(2)).^2;
blob3 = ampl3*exp(-1/(2*width3^2)*aux)';
blob4=zeros(size(blob1));
elseif n==4
 center1 = centers(:,1);
 rng(i)
width1 = 0.1 + 0.2*rand ;
rng(i+1)
ampl1 = 3.2*rand - 0.5; 

center2 = centers(:,2);
rng(i+2)
width2 = 0.1 + 0.2*rand ;
rng(i+3)
ampl2 = 3.2*rand - 0.5; 

center3 = centers(:,3);
rng(i+4)
width3 = 0.1 + 0.2*rand ;
rng(i+5)
ampl3 = 3.2*rand - 0.5; 

center4 = centers(:,4);
rng(i+6)
width4 = 0.1 + 0.2*rand ;
rng(i+7)
ampl4 = 3.2*rand - 0.5; 

aux = (xi(1,:) - center1(1)).^2 + (xi(2,:) - center1(2)).^2;
blob1 = ampl1*exp(-1/(2*width1^2)*aux)';

aux = (xi(1,:) - center2(1)).^2 + (xi(2,:) - center2(2)).^2;
blob2 = ampl2*exp(-1/(2*width2^2)*aux)';

aux = (xi(1,:) - center3(1)).^2 + (xi(2,:) - center3(2)).^2;
blob3 = ampl3*exp(-1/(2*width3^2)*aux)';

aux = (xi(1,:) - center4(1)).^2 + (xi(2,:) - center4(2)).^2;
blob4 = ampl4*exp(-1/(2*width4^2)*aux)';

end
CondVals = sigma_bg + blob1 + blob2 + blob3 + blob4;
end
