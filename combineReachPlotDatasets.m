function [mecombo,combo_se]=combineReachPlotDatasets(plotReturn1,whichTrial1,plotReturn2,whichTrial2,linecol)

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

% 4. Take sqrt to get combined standard deviation
combo_sd=sqrt(varcombo);

% 5. Divide through by sqrt(n trials) to get se
combo_se=combo_sd./sqrt(plotReturn1.n + plotReturn2.n);

% 5. Get trial-weighted combined mean
temp1=plotReturn1.(['data' num2str(whichTrial1) '_mean']);
temp2=plotReturn2.(['data' num2str(whichTrial2) '_mean']);
mecombo=combineMean(temp1{1},temp2{1},plotReturn1.n,plotReturn2.n);

% 6. Plot combo
if length(plotReturn1.time_for_x)~=length(plotReturn2.time_for_x)
    error('Assumed times were the same for plotReturn1 and plotReturn2');
end
figure();
[n,x]=cityscape_hist(mecombo,plotReturn1.time_for_x);
plot(x,n,'Color',linecol);
hold on;
[n,x]=cityscape_hist(mecombo-combo_se,plotReturn1.time_for_x);
plot(x,n,'Color',linecol);
[n,x]=cityscape_hist(mecombo+combo_se,plotReturn1.time_for_x);
plot(x,n,'Color',linecol);
xlabel('Time (sec)'); ylabel(['Reach rate (reaches per sec) ' num2str(plotReturn1.n+plotReturn2.n) ' trials']);

end

function mecombo=combineMean(me1,me2,n1,n2)

mecombo=(n1.*me1 + n2.*me2) ./ (n1+n2);

end

function varcombo=combineVar(var1,var2,n1,n2)

varcombo=((n1-1).*var1 + (n2-1).*var2)./( (n1-1) + (n2-1) );

end