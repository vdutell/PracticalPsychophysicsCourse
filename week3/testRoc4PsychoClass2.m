clear all;format compact;clc;close all
p2z= @(p) erfinv(2*p-1)*sqrt(2); %This isn't used, but it is interesting
z2p= @(z)(1-erf(z/sqrt(2)))/2;
disp('****************Rating STD analysis*******************')
LLorChisq=2; %2 for chisq and 1 for LL
datatype=3;
params0= [0 1 1.5 2 2.5 3  1 2 ]; %6 criteria and 2 dprimes for Sakitt data
if datatype==1,   %4x4 data with d'=1,2,3
    data=[84 50 16  3;13 34 34 13;3 13 34 34;0  3 16 50]
    %Note the 84 is #trials for blank that are BELOW the 1st criterion
    params0= [-1 0 1 1 2 3];  %3 criteria and 3 dprimes
elseif datatype==2,%2x4 data with 4 responses
    data=[84 50 16  3;16 50 84 97]
    params0= [0  1 2 3]; %1 criterion and 3 dprimes
elseif datatype==3,%Sakitt BS
    data=[566 192 33 9 0 0 0; 83 104 78 87 36 11 1 ;70 75 66 109 63 12 5]'
elseif datatype==4;%Sakitt data LF
    data= [585 163 32 18 2 0 0;133 83 78 58 38 9 1;91 85 81 83 40 19 1]'
elseif datatype==5;%Sakitt KD
    data=[679 78 19 6 0 0 0;149 74 94 81 2 0 0; 105 74 112 106 3 0 0]'
elseif datatype==6; %data close to Sakitt's but weaker weak stimulus
    crit=[.5:5.5]; dprime=[0:2]; sigma=1+.25*dprime;
    data= (ones(7,1).*[800 400 400]).*expected2_17(crit, dprime, sigma,1)
end
findzero=sum((sum(data')==0)) %count how many rows are all zero
if findzero>0; data=data+.000001;end  %add .000001 if there are all zero rows 
%Do a Monte Carlo simulation for the blank stimulus

[Nresp Nstim]=size(data)
CumSum=cumsum(data);
Total=CumSum(end,:); %total number of trials at each level
CumProb=1-CumSum./(ones(Nresp,1)*Total);%this make probs to RIGHT of criterion
CumProb=[ones(1,Nstim); CumProb] %the extra row is for pretty plot

if Nresp==2   %This is for non-Sakitt data
    prob=data(2,:)./CumSum(end,:)
    z=erfinv(2*prob-1)*sqrt(2)
    dprime2resp=diff(z)
    disp('***   Palamedes 1AFC but note order of HF is reversed ***')
    [dp C lnBeta Pc]=PAL_SDT_1AFC_PHFtoDP(prob(1:2));
    DprimePAL=dp'
end
float=-1; %-1 is for floating ROC slope 0 fixed slope.info(2)is the only item we may want to modify
info(1)=6;          %6 means separate d' values (I think)
info(2)=float;      %-1 Float ROC slope   0 leave slope at 0
info(3)=LLorChisq;  %LLorChisq
info(4)=1;          %dprimeType  =1 for float d'
offset=0;           %for some sort of transducer option
if info(2)==-1; params0=[params0 .1];end %the 0.1 is an initial guess for 9th parameter
if LLorChisq==1, %1 is max likelihood,  2 for chi square (typically we do chi square)
    disp('this is LogLikelihood')
    params = fminsearch('SetRocSimpler',params0,[], offset, data,info)
    %   [LL EstSimpleLL]=SetRocSimpler(params, offset, data,info);
    [LL EstSimpleLL]=SetRoc_w_sigma(params, offset, data,info);
end
%[params,chisq1,f,EXITFLAG,OUTPUT,LAMBDA,j] = lsqnonlin('SetRocSimpler',params0,[],[],[], offset, data,info);
[params,chisq1,f,EXITFLAG,OUTPUT,LAMBDA,j] = lsqnonlin('SetRoc_w_sigma',params0,[],[],[], offset, data,info);
disp('Minimizing chisq')
chisq1
params
dprime=[0 params(Nresp:Nresp+Nstim-2)]
sigma=1+params(Nresp+Nstim-1)*dprime

j=full(j);  %jacobian wasn't known to Kendrick Kay or to Goef Boynton and most others
cov=inv(j'*j);
SE=sqrt(diag(cov))'   %But one should always check with Monte Carlo or Bootstrap
correl=cov./(SE'*SE)
%[err, expectSimple]=SetRocSimpler(params, offset, data,info);
%[err, expectSimple]=SetRoc_w_sigma(params, offset, data,info);
[a,a,expectSimple]=SetRoc_w_sigma(params, offset, data,info);
expect=num2str(expectSimple,3)
data
chisq=sum(sum((expectSimple-data).^2./expectSimple))
DegFree=(Nresp-1)*Nstim-length(params0)

figure(1); 
PlotSakittROC_data
figure(3)
PlotSakittGaussians(params, SE)