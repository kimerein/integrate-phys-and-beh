function combineReachPlotDatasets(plotReturn1,whichTrial1,plotReturn2,whichTrial2)

if ~ismember(whichTrial1,[1 2])
    error('whichTrial1 must be 1 or 2');
end
if ~ismember(whichTrial2,[1 2])
    error('whichTrial2 must be 1 or 2');
end

% Standard error of raw reaching was calculated as
% std(raw reach rate) / sqrt(n trials)

% 1. Get back to standard deviation
temp1=plotReturn1.(['data' num2str(whichTrial1) '_se']);
sd1=temp1{1}.*sqrt(plotReturn1.n);
temp2=plotReturn2.(['data' num2str(whichTrial2) '_se']);
sd2=temp2{1}.*sqrt(plotReturn2.n);

% 2. Variance from standard deviation
var1=sd1.^2;
var2=sd2.^2;

% 3. Combine variances
varcombo=combineVar(var1,var2,plotReturn1.n,plotReturn2.n);

% 4. Take sqrt to get combined standard error
combo_se=sqrt(varcombo);

% 5. Get trial-weighted combined mean
temp1=plotReturn1.(['data' num2str(whichTrial1) '_mean']);
temp2=plotReturn2.(['data' num2str(whichTrial2) '_mean']);
mecombo=combineMean(temp1{1},temp2{1},plotReturn1.n,plotReturn2.n);

end

function mecombo=combineMean(me1,me2,n1,n2)

mecombo=(n1.*me1 + n2.*me2) ./ (n1+n2);

end

function varcombo=combineVar(var1,var2,n1,n2)

varcombo=((n1-1).*var1 + (n2-1).*var2)./( (n1-1) + (n2-1) );

end