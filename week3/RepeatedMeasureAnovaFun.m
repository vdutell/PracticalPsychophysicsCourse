function dif=RepeatedMeasureAnovaFun(params, data, SD, option)
%Used by lsqnonlin to getthe difference between expected and observed
[nrows ncolumns] = size(data);
nparams=length(params);
if nparams>ncolumns,
    if option==1, endparams=params(ncolumns+1:end);
        else endparams=params(ncolumns+1)*[1 1]+params(ncolumns+2)*[1 -1];
    end
    treatment=[-sum(endparams) endparams ];
else treatment=zeros(1,nrows);
end
expect= ones(nrows,1)*params(1:ncolumns)+ treatment'*ones(1, ncolumns);
dif= (expect- data)./SD;%the deviation (SE) units for leastsq