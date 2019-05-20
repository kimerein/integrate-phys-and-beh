function plotRTs_contingentOnTrialType(tbt,whichReach,useAsCue,metadata,contingency1,contingency2,contingencyName)

% if contingency2 is empty, will not consider this contingency
% will just plot contingency1==true trials

nbins=200; % for histograms
zscore_RTs=1; % if 1, will use z-scored instead of raw reaction times

reactionTimes=getPairedReactionTimes(tbt,whichReach,useAsCue,metadata,zscore_RTs,[],[]);

plotHist(reactionTimes,nbins,['Compare RTs dependent on: ' contingencyName],contingency1,contingency2);
plotCDF(reactionTimes,nbins,['Compare RTs dependent on: ' contingencyName],contingency1,contingency2);
testRanksum(reactionTimes,contingency1,contingency2);

end

function testRanksum(reactionTimes,testcond1,testcond2)

if isempty(testcond1) || isempty(testcond2)
    disp('testcond is empty');
    return
else
    if size(testcond1,1)>1
        % turn into row vector
        testcond1=testcond1';
    end
    if size(testcond2,1)>1
        % turn into row vector
        testcond2=testcond2';
    end
end

[p,h]=ranksum(reactionTimes(testcond1==true),reactionTimes(testcond2==true));
disp('p-value');
disp(p);

end

function plotHist(reactionTimes,bins,tit,testcond1,testcond2)

if size(testcond1,1)>1
    % turn into row vector
    testcond1=testcond1';
end
if ~isempty(testcond2)
    if size(testcond2,1)>1
        % turn into row vector
        testcond2=testcond2';
    end
end

[n,x]=histcounts(reactionTimes(testcond1==true),bins);
x_backup=x;
[n,x]=cityscape_hist(n,x);
figure();
plot(x,n./nansum(n),'Color','k');
xlabel('RT (sec)');
ylabel('Count');
title(tit);

% hold on;
% [parmhat,parmci]=lognfit(reactionTimes(testcond1==true & ~isnan(reactionTimes)));
% % plot fit
% Y=lognpdf(nanmin(x):0.01:nanmax(x),parmhat(1),parmhat(2));
% plot(nanmin(x):0.01:nanmax(x),Y*(x(2)-x(1))*0.55,'Color','c');

if ~isempty(testcond2)
    [n,x]=histcounts(reactionTimes(testcond2==true),x_backup);
    [n,x]=cityscape_hist(n,x);
    hold on;
    plot(x,n./nansum(n),'Color','r');
    leg={'testcond 1','testcond 2'};
    legend(leg);
end

end

function plotCDF(reactionTimes,bins,tit,testcond1,testcond2)

if size(testcond1,1)>1
    % turn into row vector
    testcond1=testcond1';
end
if ~isempty(testcond2)
    if size(testcond2,1)>1
        % turn into row vector
        testcond2=testcond2';
    end
end

[n,x]=histcounts(reactionTimes(testcond1==true),bins);
x_backup=x;
x_mids=nanmean([x(1:end-1); x(2:end)],1);
cond1_cdf=accumulateDistribution(n);
figure();
plot(x_mids,cond1_cdf./nanmax(cond1_cdf),'Color','k');
xlabel('CDF of change in RT');
ylabel('Count');
title(tit);

if ~isempty(testcond2)
    [n,x]=histcounts(reactionTimes(testcond2==true),x_backup);
    hold on;
    cond2_cdf=accumulateDistribution(n);
    plot(x_mids,cond2_cdf./nanmax(cond2_cdf),'Color','r');
    leg={'testcond 1','testcond 2'};
    legend(leg);
end


end

function cdf=accumulateDistribution(data)

cdf=nan(size(data));
for i=1:length(data)
    cdf(i)=nansum(data(1:i));
end

end