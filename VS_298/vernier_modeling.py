from pylab import *
import scipy.io
from sklearn import datasets, linear_model

f = open('data51by14.txt')
data = []
for line in f:
    if not ('Data' in line):
        line_data = line.split()
        for ld in line_data:
            ld = ld.replace('\r\n','').replace('];','')
            data.append(float(ld))
f.close()
data = array(data)
data = data.reshape((51,14))

pres = []
posts = []
pprs = []

for subj_num in range(data.shape[1]):
    subj_data = data[:,subj_num]
    pre = mean(subj_data[1:5]); pres.append(pre)
    post = mean(subj_data[-4:]); posts.append(post)
    ppr = pre/post; pprs.append(ppr)

    x = array([1,2,3,4,47,48,49,50])
    y = concatenate((subj_data[1:5],subj_data[-4:]))
    x = x.reshape((len(x),1))
    y = y.reshape((len(y),1))

    regr = linear_model.LinearRegression(fit_intercept=True)
    regr_xint = linear_model.LinearRegression(fit_intercept=True)

    regr.fit(x,y)
    regr_xint.fit(y,x)

    slope = regr.coef_[0,0]
    slope_xint = 1.0/regr_xint.coef_[0,0]

    y_inter = regr.intercept_[0]
    x_inter = regr_xint.intercept_[0]

    
    xs = arange(len(subj_data))
    ys1 = slope*xs + y_inter
    ys2 = slope_xint*(xs - x_inter)
    #ys3 = slope*(xs - x_inter)

    plot(xs,ys1,'g-', 
         xs,ys2,'b-',
         xs,subj_data,'ro')
    show()

'''
for i=Nsubj:Nsubj  %We set Nsubj=1 for this demo program
    Data=DataTraining(start:end,i)';
    Pre=mean(Data(2:6)); Post(i)=mean(Data(end-4:end));
    PPR(i)=Pre/Post(i);
    params0=[PPR(i) Post(i) 15];
    
    LB=[.5 2 4]; UB=[5 20 40]; FixParam=[0 0 0];
    if whichPost==1, iter=9; % This is for post at 50.
    else             iter=8; % This is for post at infinity.
    end
    SD=1;
    [params,SSE,~,~,~,~,j]=lsqnonlin(@ExponFun3,params0,LB,UB,options, xfit,Data,SD,iter,doLog);
    %   on this initial iteration three is no need to get covar
    Weber(i)=sqrt(SSE/DegFree)%start=2 to avoid 1st staircase
    SD3=SD*Weber(i); %SD=1 Was the initial guess
    
    [params,chisqF,~,~,~,~,j]=lsqnonlin(@ExponFun3,params,LB,UB,options, xfit,Data,SD3,iter,doLog,FixParam);
    chisqFAll(i)=chisqF
    j=full(j)' %the 'full' command converts sparse matrices into regular ones
    alpha=j*j';     %this is NR (Numerical Recipes) 15.4.8  (Note our j is NR's A)
    covar=inv(alpha) %see NR and simulations to clarify covariance matrix
    se=sqrt(diag(covar))';   %Standard error
    paramsAll(i,:)=params
    seAll(i,:)=se
    SDall(i,:)=DataTraining(:,i)'*Weber(i);  %to include the first datum
    correl=covar./(se'*se)
    correlAll(i,:)=[correl(1,2) correl(1,3) correl(2,3)];
    [dif,E]=ExponFun3(params,x,DataTraining(:,i)',SDall(i,:),iter,doLog,FixParam); %this x is full 51
end
'''
