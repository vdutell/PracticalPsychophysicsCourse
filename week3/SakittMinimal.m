%SakittMinimal
clear all;close all;clc
%z2p= @(z)(1-erf(z/sqrt(2)))/2;
disp('***********Analysis of Sakitt''s data***************')
data=[566 192 33 9 0 0 0; 83 104 78 87 36 11 1 ;70 75 66 109 63 12 5]'
info=1; %=1 for first 6 params to be criteria, =.1 for params to be diff of criteria
params0=[1 2 3 4 5 6 1 2 .5]*info  %random guesses for initial parameters
[params,chisq,~,~,~,~,j] = lsqnonlin('SetRoc_Sakitt',params0,[],[],[], data,info);
j=full(j);  %jacobian is really nifty, very close to linreg design matrix
cov=inv(j'*j);
params
SE=sqrt(diag(cov))'
correl=cov./(SE'*SE)
figure(1); PlotSakittROC_data   %plots the ROC curves on p and z axes
figure(3)
if info==.1, params(1:6)=cumsum(params(1:6));end
PlotSakittGaussians(params, SE)