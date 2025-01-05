
clear all;
clc;


load CEM_data b_CEM%load the data computed by using CEM model
load Mobius_data b_mob


n_s=size(b_CEM,2);%number of samples
err=b_CEM-b_mob;%modeling error
err_mean=(1/n_s).*(sum(err,2));%modeling error mean 
err_c= err - (err_mean*ones(1,n_s));%centralizing the sample

err_cov=(1/(n_s-1)).*(err_c*err_c');%modeling error covariance
figure(3)
imagesc(reshape(log(abs(err_mean)),15,15))
colorbar;
title('Mean of the error with coarse FEM mesh vs Mobius')
set(gca,'FontSize',15)
figure(4)
imagesc(log(abs(err_cov)))
colorbar
title('Covariance of the error with coarse FEM mesh vs Mobius')
set(gca,'FontSize',15)
save Model_err_FEMcoarsevsMobius err_mean err_cov

