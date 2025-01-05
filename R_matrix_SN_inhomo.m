% Load the current/voltage experimental data
function [LL]=R_matrix_SN_inhomo(d,s,sig_back_CEM)
[CurrentPattern, Uel]=load_data(d,s);
load Real_data_matrix_homo
sig_opt=sig_CEM;
r=14;%radius of the tank in the experimental setup

sig_back=0.3;%background conductivity in the measurements
amp=2;%amplitude of the applied currents

I = CurrentPattern;
dV = Uel;
sig_back_CEM=1.5;
% The data are given in terms of voltage differences between adjacent
% electrodes. Solve for electrode voltages

L = 16;
V = zeros(L,L-1);
for k = 1:L-1
    % loop over current feeds
    V_cum = zeros(L,1);
    for j = 2:L
        V_cum(j) = V_cum(j-1) + dV(j-1,k);
    end
    V1 =   1/L*sum(V_cum);
    V(:,k)= V1*ones(L,1)- V_cum;
end
%V = -V;

scale=sig_opt/(sig_back_CEM);

V_sc=((sig_back*r*scale)/(sig_back_CEM*amp))*V;%scaled voltages
% Trigonometric current patterns

L = 16;
aux = (0:L-1)';
J = [cos(2*pi/L*aux*(1:L/2)),sin(2*pi/L*aux*(1:L/2-1))];
dd = sqrt(2/L)*ones(1,L-1);
dd(8) = sqrt(1/L);
J = J*diag(dd);

%plot the real measurements-- currents and voltages
% for i=1:L-1
%         
%      figure(2)
%      subplot(4,4,i)
%      plot(V(:,i),'LineWidth',2.5)
%      hold on 
%      plot(I(:,i),'LineWidth',2.5)
%      hold off
%      legend('V','I')
%      title(['V_{', num2str(i),'}'])
% end
%  suptitle('Voltage and current measurements')


 
coeff_inv=J'*I;%inverse of the coefficients
V_T=V/coeff_inv;%V_(trig)=V_a*(Phi'*I_a)^{-1}

V_T_sum=sum(V_T);%Recheck whether the sum is zero


R=(J'*V_T);%resistance matrix on the model tank
%of radius '14' and background conductivity 0.3

V_T_sc=V_sc/coeff_inv;%V_(trig)=V_a*(Phi'*I_a)^{-1}

V_T_sc_sum=sum(V_T_sc);%Recheck whether the sum is zero


R_sc=(J'*V_T_sc);%resistance matrix on the model tank
%of radius '1' and background conductivity 1.5


LL=R_sc\eye(size(R_sc));

% figure(1)
% plot(diag(R),'k.','MarkerSize',25)
% set(gca,'FontSize',15)
% set(gca,'Xtick',1:15)
% grid

 

%  inc=0:0.001:0.2;
%  for i=1:length(inc)
% z = (0.001+inc(i))*ones(1,L); 
% R_CEM=est_cont_imp(z);
% err_norm(i)=norm(R_sc-R_CEM,'fro')/norm(R_CEM,'fro');
%  end
%  figure(4)
%  plot(inc,err_norm,'k.','MarkerSize',25)
% set(gca,'FontSize',15)
% grid
% [~,ind]=min(err_norm)
% z=(0.001+inc(ind))*ones(1,L);%finding the contact impedance where 
% %Frobenious norm between R_CEM and R_data is minimum

filename = sprintf('Real_data_matrix_inhomo_%d_%d.mat', d,s);
save(filename, 'J', 'V_T', 'V_sc', 'R_sc', 'LL')
end


    


