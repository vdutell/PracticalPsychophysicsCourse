%RepeatedMeasureAnova_effectSize
%https://statistics.laerd.com/statistical-guides/repeated-measures-anova-statistical-guide-2.php
clear all;clc
F2p  = @(F,df1,df2) betainc(1 ./(1+F*df1/df2),df2/2,df1/2 ); %for 2tail t2p wth df1=1
data=[45 42 36 39 51 44; 50 42 41 35 55 49; 55 45 43 40 59 56]
[Nrows Ncolumns] = size(data);
params6=mean(data)   %the 6 params for the 6 subjects
disp('These are the parameters for the 6 subject means')
Expect6=ones(Nrows,1)*params6;  %fit with 6 parameters for the 6 subjects%
params2=mean((data-Expect6)')  %Two more parameters (mean=0)
disp('These are params for the two treatments. Total means were already subtracted.')
disp('Notice how easy it was to get params, but the challenge is to get F')
Treatment=params2'*ones(1,6);
Expect8=Expect6+ Treatment; %subtract treatment
error=data-Expect8;
DegFree=Nrows*Ncolumns-Nrows-Ncolumns+1 %6 subject and two treatments
SSTot=sum(sum(data))
SSE6param=sum(sum((data-Expect6).^2)) %SSE for 6 param fit (subject differences)
SSE8param=sum(sum(error.^2))   %SSE for 8 param fit
Var=SSE8param/DegFree
SD= sqrt(Var)
dif=(Expect6-Expect8)/SD
disp('This is the item that is needed for lsqnonlin')
disp('It is quite easy to generate dif using lsqnonlin')
chisqF2=sum(sum(dif.^2))
F=chisqF2/(Nrows-1)  %The -1 is because the mean was already removed
p = F2p(F,2, 10)  %the 2 is degfree of the difference (numerator), 10 is DF of denom
SSEdif=SSE6param-SSE8param; %this is the SSE difference due to treatment params
EffectSize=SSEdif/SSE6param
EffectSize=SSEdif/(SSE8param+SSEdif)  %this is a common way of writing it

%********************************
disp(' ')
disp('Now do it with lsqnonlin. The first call to lsqnonlin is to get SD.')
Nparams=Nrows+Ncolumns-1; %the -1 is becausee the mean was removed with subjects
params0=randn(1,Nparams); %just a totally random initial guess
options=optimset('display','off');%iter
SD1=1;
info=1;  %Info=1 Use separate levels for the two treatments
         %Info=2 Use mean and separation for the two treatments
[params,SSE,~,~,~,~,j] = ...
    lsqnonlin('RepeatedMeasureAnovaFun', params0,[],[],options, data, SD1,info);
params
SD=sqrt(SSE/DegFree)   %Knowing SSE for SD=1 enables us to get true SD. 
for iter=1:2,
    disp(' ')
    disp(['Iteration #' num2str(iter)])
    if iter==2, disp('Now do it without the treatment params');
        params=params(1:6);  %the first iter was the FULL model.         
    end  %for no effect of treatment
    [params,chisqF,~,~,~,~,j] = ...
        lsqnonlin('RepeatedMeasureAnovaFun', params,[],[],options, data, SD,info);
    j=full(j);  %jacobian is really nifty, very close to linreg design matrix
    cov=inv(j'*j);
    params
    SE=sqrt(diag(cov))'
    chisqFAll(iter)=chisqF;
    if iter==1, 
        endparams=params(Ncolumns+1:end);
        treatment=[-sum(endparams) endparams ];
        disp(['The 3 treatments were: ' num2str(treatment,2)])
        correl=cov./(SE'*SE), 
    end
end
disp('Notice the correlations involving subjects is zero, but not between treatments.')
disp('Is there a better way of parameterizing the treatments?')
disp('All those zeros mean the params and SE of subjects is same across iterations')
chisqFAll
disp('Note that the full model always has chisqF = DF = 10.')
disp('Thus there is no chisq goodness of fit if SD is estimated from data!')
disp(' '); disp('Now take the difference of the two iterations to get chisqF and F')
disp('F = chisqF divided by degfreedom that is max #good params - #of model params')
chisqF2=diff(chisqFAll)  %The difference between 8 and 6 param fits
F= chisqF2/(Nrows-1) %where Nrows is Total number of treatments (3 in our example)
p = F2p(F,2, 10)  %2 & 10 are the DF of numerator & SD (denominator)
EffectSize = diff(chisqFAll)/chisqFAll(2)
disp('This is the same as the standard way in the ANOVA tutorial')




