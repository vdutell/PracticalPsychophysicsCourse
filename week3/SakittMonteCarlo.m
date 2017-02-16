%SakittMonteCarlo
clear all;format compact;clc;close all
%p2z= @(p) erfinv(2*p-1)*sqrt(2); %ignore 'ProbExact' name
disp('****************Rating STD analysis*******************')
%LLorChisq=2; %2 for chisq and 1 for LL
info=.1;  %=1 for first 6 params to be criteria, =.1 for params to be diff of criteria 
crit=[.5:5.5]*info ; 
dprime=[0:2]; sigma=1+.25*dprime;
expected2_17(crit, dprime, sigma,1);
data= (ones(7,1)*[800 400 400]).*expected2_17(crit, dprime, sigma,1);
params0=[crit dprime(2:3) .25];
% findzero=sum((sum(data')==0)) %count how many rows are all zero
% if findzero>0; data=data+.0001;end  %add .0001 if there are all zero rows
%Do a Monte Carlo simulation for the blank stimulus
data=data+.00001  %This is done to avoid ALL the rows being zero. 
Nsim=100;  %number of Monte Carlo simulations
Nresp=7; %number of categories
Nstim=3; %number of stimuli (B, W S)
%Make Nsim datasets (idealized Sakitt data)
    MCdataB=randn(800,Nsim);
    for k=1:Nresp;tot(k,:)=sum(MCdataB<(k-.5));end
    difftotB=diff([zeros(1,Nsim); tot]);
    MCdataW=dprime(2)+randn(400,Nsim)*(1+.25*dprime(2));
    for k=1:Nresp;tot(k,:)=sum(MCdataW<(k-.5));end
    difftotW=diff([zeros(1,Nsim); tot]);
    MCdataS=dprime(3)+randn(400,Nsim)*(1+.25*dprime(3));
    for k=1:Nresp;tot(k,:)=sum(MCdataS<(k-.5));end
    difftotS=diff([zeros(1,Nsim); tot]);
float=-1; %-1 is for floating ROC slope 0 fixed slope.info(2)is the only item we may want to modify
options=optimset('display','off');%iter
for kdat=1:Nsim;
    data=[difftotB(:,kdat) difftotW(:,kdat) difftotS(:,kdat)]+.001;
    CumSum=cumsum(data);
    Total=CumSum(end,:); %total number of trials at each level
    CumProb=1-CumSum./(ones(Nresp,1)*Total);%this make probs to RIGHT of criterion
    CumProb=[ones(1,Nstim); CumProb]; %the extra row is for pretty plot
    [params,chisq,~,~,~,~,j] = lsqnonlin('SetRoc_Sakitt',params0,[],[],options, data, info);%options
    paramsAll(kdat,:)=params;
    j=full(j);  %jacobian is really nifty, very close to linreg design matrix
    cov=inv(j'*j);
    SE=sqrt(diag(cov))';
    params;
    correlAll(kdat,:,:)=cov./(SE'*SE);
    SEAll(kdat,:)=SE;
    chisqAll(kdat)=chisq;
end
Nsim
params0
paramsAll;
meanParams=mean(paramsAll)
meanSE_LSQ=mean(SEAll)
rmsSE_LSQ=sqrt(mean(SEAll.^2))
STDparams_MC=std(paramsAll)
SEparams_MC=STDparams_MC/sqrt(Nsim)
meanCorrel=squeeze(mean(correlAll));
SECorrel=squeeze(std(correlAll));
for k1=1:9;for k2=k1+1:9  %this is the cute way of also showing SE
     meanCorrel(k2,k1)=SECorrel(k1,k2);
    end, end
SE_meanCorrel_MC=meanCorrel
info   %To indicate whether criteria or delta criteria were used
figure(1)
PlotSakittROC_data
title('Plotted data is from final sample')
figure(3)
if info==.1, params(1:6)=cumsum(params(1:6));end
PlotSakittGaussians(params, SE)