%For plotting Sakitt data
z2p= @(z)(1-erf(z/sqrt(2)))/2;
testz2p=z2p(2)  %Just checking to see if erf is okay
[Nresp Nstim]=size(data);
CumSum=cumsum(data);
Total=CumSum(end,:); %total number of trials at each level
CumProb=1-CumSum./(ones(Nresp,1)*Total);%this make probs to RIGHT of criterion
CumProb=[ones(1,Nstim); CumProb] %the extra row is for pretty plot
dprime=[0 params(Nresp:Nresp+Nstim-2)]
sigma=1+params(Nresp+Nstim-1)*dprime
float=1;
datatype=6;
findzero=sum((sum(data')==0)); %count how many rows are all zero

subplot(1,2,1)
plot(CumProb(:,1),CumProb(:,2:Nstim),'*', [0 1],[0 1], '-k');hold on
x=-6:.1:6;
Nparams=length(params);
clr='br';
for k=2:3,  %remember that k=1 is the blank
    if Nparams==9; y=(x-dprime(k))./(1+params(9)*dprime(k));
    else y=(x-dprime(k));
    end
    py=z2p(y);
    px=z2p(x);
    plot(px,py,clr(k-1)); hold on
end

%CumProb(:,1),CumProb(:,2:Nstim),'-',
z=erfinv(2*CumProb-1)*sqrt(2);  %each z row is for one criterion5
%dprime=z-z(:,1)*ones(1,Nstim); %recall the 1st column is for blanks
xlabel('probability of responding to blank')
ylabel('probability of responding to signal')
text(.1,.1,['dprime=' num2str(dprime,2)])
text(.1,.05,['sigma =' num2str(sigma,2)])

subplot(1,2,2)
%plot(z(:,1), z(:,2:Nstim),'-',z(:,1), z(:,2:Nstim),'*');hold on
plot(z(:,1), z(:,2:Nstim),'*');hold on
plot([-4 0],[-4 0],'k-',[-6 0],[0 0],'k-'); hold on
if float==0, params(9)=0; end
for k=0:1,  %plot the two d's
    dprime=params(Nresp+k);
    y1=(-6+dprime)/(1+params(9)*dprime);
    y2=(-0+dprime)/(1+params(9)*dprime);
    plot([-6 0], [y1 y2],clr(k+1))
end
if datatype==6; zz=-5.9;
else zz=-2.9;
end
xlabel('z-score of responding to blank')
ylabel('z-score of responding to signal')
text(zz, 1.75,['dprime for W=' num2str(params(7),3)  '+-' num2str(SE(7),2)])
text(zz, 1.5,['dprime for S=' num2str(params(8),3)  '+-' num2str(SE(8),2)])
text(zz, 1.25,['slope param=' num2str(params(9),3)  '+-' num2str(SE(9),2)])
text(zz, 1,['chi square  =' num2str(chisq,3) '+-' num2str(sqrt(2*chisq),2)])
text(zz, .75,['deg freedom = 3*6-9=9 '])
if findzero>0, text(zz,.5,['correct df=' num2str(3*findzero) ' smaller']);end
if datatype==6, text(-3.5,-.0,'dprime= ');end
for k=1:2 text(-k,0, num2str(k));end
axis([zz 0 -4 2])