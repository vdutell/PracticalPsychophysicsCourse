function PlotGaussians(params, SE);
%plot Sakitt
x=-3:.1:6;
for k=1:3;
    if k==1, gauss= exp(-x.^2/2);
        sigma=1;
    else
        sigma=1+params(9)*params(5+k);
        gauss= exp(-(x-params(5+k)).^2/2/sigma^2)/sigma;
        text(params(5+k)-.3,1/sigma+.02,[num2str(params(5+k),2) '+-' num2str(SE(5+k),1)])
    end
    plot(x, gauss,'b'); hold on
end
for k=1:6;
    plot(params(k)*[1 1], [0 1.1], 'k'); hold on
    text(params(k)-.4,1.13,[num2str(params(k),2) '+-' num2str(SE(k),1)])
end
xlabel('Internal response (Sakitt''s counting quanta)')
ylabel('probability of response')