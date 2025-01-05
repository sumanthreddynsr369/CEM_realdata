
%this program generates the sample of radially symmetric conductivities

%%Input:
%n_space---number of radial discretization points
%n_sample---number of samples to be generated
%gamma---std in the generative Gaussian model
%alp---correlation parameter
%sig_bg---background conductivity

%%Output:
%a---radial discretization
%Sig_samp_radial---radially symmetric conductivity sample (size ((n_space+1)*n_samp))

function [a,Sig_samp_radial]=cond_sample_radial(n_space,n_sample,gamma,alp,sig_bg,visualize)

a_in=0;%domain discretization initial value
a_fin=1;%radius of the disc
%space discretization (a_fin-a_in)/n_num



[Xi,a]=xi_sample_radial(n_space,n_sample,a_in,a_fin,gamma,alp);
%generated log of the radially symmetric smooth conductivity samples
std_Xi=std(Xi);%mean of the sample
mean_Xi=mean(Xi);%std of the sample

%visualize three random draws of generated sample 
if strcmp(visualize,'yes')
figure(1)

%%
XX=[a,fliplr(a)];        %create continuous x value array for plotting
YY=[(mean_Xi- 2*std_Xi),fliplr((mean_Xi+ 2*std_Xi))];      %create y values for out and then back
fill(XX,YY,1,....
        'facecolor','blue', ...
        'edgecolor','none', ...
        'facealpha', 0.3);
hold on
  plot(a,mean_Xi,'LineWidth',2);
  xlabel('Radius')
  ylabel('log(\sigma/\sigma_b)')
  title('Sample of log(\sigma/\sigma_b) value')
  set(gca,'FontSize',18)
  hold on

for j=1:3
    i=randi(num_sample);
    pp(j)=i;
  plot(a,Xi(i,:),'LineWidth',2);
  set(gca,'FontSize',18)
  
hold on
end
  
legend('Two standard deviation range','Mean')
%%
hold off
end


  Sig_samp_radial=sig_bg.*exp(Xi);%the positive conductivities
% 

mean_sig=mean(Sig_samp_radial);%mean
std_sig=std(Sig_samp_radial);%std

% std_sig=sqrt((exp(std_Xi.^2)-1).*exp(2.*mean_Xi + (std_Xi.^2)));
if strcmp(visualize,'yes')
figure(2)
YY=[(mean_sig - 2*std_sig),fliplr((mean_sig + 2*std_sig))];      %create y values for out and then back
fill(XX,YY,1,....
        'facecolor','blue', ...
        'edgecolor','none', ...
        'facealpha', 0.3);
hold on
  plot(a,mean_sig,'LineWidth',2);
  xlabel('Radius')
   ylabel('Conductivity value')
   title('Sample of Conductivities')
  set(gca,'FontSize',18)
  hold on

for j=1:3
     i=pp(j);
  plot(a,Sig_samp_radial(i,:),'LineWidth',2);
  set(gca,'FontSize',18)
  
hold on
end
  
legend('Two standard deviation range','Mean')
hold off
end
end