function [reactionTimes,RT_pairs]=getRTpairs_contingentOnTrialType(tbt,whichReach,useAsCue,metadata,contingency1,contingency2,prevOrCurrTrial1,prevOrCurrTrial2,contingencyName,selectTrials,nApart,doPlot)

% if contingency2 is empty, will not consider this contingency
% will just plot contingency1==true trials

nbins=100; % for histograms
jitter=0.03; % for scatter plots
addSlopeLine=0; % for scatter plots, will add fit to slope if 1
zscore_RTs=0; % if 1, will use z-scored instead of raw reaction times

[reactionTimes,RT_pairs,tbt]=getPairedReactionTimes(tbt,whichReach,useAsCue,metadata,zscore_RTs,selectTrials,nApart);
if isempty(nApart)
    nApart=1;
end

if length(selectTrials)>1
    % selectTrials is a subset of trials
    contingency1=contingency1(selectTrials==1);
    if ~isempty(contingency2)
        contingency2=contingency2(selectTrials==1);
    end
end
    
% does contingency 1 apply to the first or second trial in the pair?
if prevOrCurrTrial1==1
    % contingency applies to first trial in pair, i.e., previous trial
    contingency1=contingency1(1:end-nApart);
elseif prevOrCurrTrial1==2
    % contingency applies to second trial in pair, i.e., current trial
    contingency1=contingency1(nApart+1:end);
else
    error('prevOrCurrTrial must be 1 or 2');
end
if ~isempty(contingency2)
    % does contingency 2 apply to the first or second trial in the pair?
    if prevOrCurrTrial2==1
        % contingency applies to first trial in pair, i.e., previous trial
        contingency2=contingency2(1:end-nApart);
    elseif prevOrCurrTrial2==2
        % contingency applies to second trial in pair, i.e., current trial
        contingency2=contingency2(nApart+1:end);
    else
        error('prevOrCurrTrial must be 1 or 2');
    end
end

if doPlot==1
    plotHist(RT_pairs,nbins,['Compare RT change from prev to curr trial dependent on: ' contingencyName],contingency1,contingency2);
    plotCDF(RT_pairs,4*nbins,['CDF ' contingencyName],contingency1,contingency2);
    plotScatter(RT_pairs,['Prev vs curr trial RT dependent on: ' contingencyName],contingency1,contingency2,jitter,addSlopeLine);
    if ~isempty(contingency2)
        % only compare if 2 groups
        testRanksum(RT_pairs,contingency1,contingency2);
    end
end

end

function testRanksum(RT_pairs,testcond1,testcond2)

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

disp('median of contingency 1');
disp(nanmedian(RT_pairs.rt_change(RT_pairs.real_rt_pair==true & testcond1==true),2));
disp('median of contingency 2');
disp(nanmedian(RT_pairs.rt_change(RT_pairs.real_rt_pair==true & testcond2==true),2));

[p,h]=ranksum(RT_pairs.rt_change(RT_pairs.real_rt_pair==true & testcond1==true),RT_pairs.rt_change(RT_pairs.real_rt_pair==true & testcond2==true));
disp('p-value');
disp(p);

end

function plotScatter(RT_pairs,tit,testcond1,testcond2,jitter,addSlopeLine)

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

figure();
usePrev=RT_pairs.prev_trial_rt+rand(size(RT_pairs.prev_trial_rt)).*jitter;
useCurr=RT_pairs.curr_trial_rt+rand(size(RT_pairs.curr_trial_rt)).*jitter;
s=scatter(usePrev(RT_pairs.real_rt_pair==true & testcond1==true),useCurr(RT_pairs.real_rt_pair==true & testcond1==true),[],'k','filled');
pause;
m=get(s,'MarkerHandle');
alpha=0.3;
m.FaceColorData=uint8(255*[0;0;0;alpha]);
xlabel('RT previous trial in sec');
ylabel('RT current trial in sec');
title(tit);
if addSlopeLine==1
   dlm=fitlm(usePrev(RT_pairs.real_rt_pair==true & testcond1==true)',useCurr(RT_pairs.real_rt_pair==true & testcond1==true)','Intercept',false);
   xpoints=0:0.001:nanmax(usePrev(RT_pairs.real_rt_pair==true & testcond1==true));
   ypoints=dlm.Coefficients.Estimate.*xpoints;
   hold on;
   plot(xpoints,ypoints,'Color','k');
end

if ~isempty(testcond2)
    hold on;
    s=scatter(usePrev(RT_pairs.real_rt_pair==true & testcond2==true),useCurr(RT_pairs.real_rt_pair==true & testcond2==true),[],'r','filled');
    pause;
    m=get(s,'MarkerHandle');
    alpha=0.3;
    m.FaceColorData=uint8(255*[1;0;0;alpha]);
    leg={'testcond 1','testcond 2'};
    legend(leg);
    if addSlopeLine==1
        dlm=fitlm(usePrev(RT_pairs.real_rt_pair==true & testcond2==true)',useCurr(RT_pairs.real_rt_pair==true & testcond2==true)','Intercept',false);
        xpoints=0:0.001:nanmax(usePrev(RT_pairs.real_rt_pair==true & testcond2==true));
        ypoints=dlm.Coefficients.Estimate.*xpoints;
        hold on;
        plot(xpoints,ypoints,'Color','r');
    end
end

end

function plotCDF(RT_pairs,bins,tit,testcond1,testcond2)

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

[n,x]=histcounts(RT_pairs.rt_change(RT_pairs.real_rt_pair==true & testcond1==true),bins);
x_backup=x;
x_mids=nanmean([x(1:end-1); x(2:end)],1);
cond1_cdf=accumulateDistribution(n);
figure();
plot(x_mids,cond1_cdf./nanmax(cond1_cdf),'Color','k');
xlabel('CDF of change in RT');
ylabel('Count');
title(tit);

if ~isempty(testcond2)
    [n,x]=histcounts(RT_pairs.rt_change(RT_pairs.real_rt_pair==true & testcond2==true),x_backup);
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

function plotHist(RT_pairs,bins,tit,testcond1,testcond2)

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

[n,x]=histcounts(RT_pairs.rt_change(RT_pairs.real_rt_pair==true & testcond1==true),bins);
x_backup=x;
[n,x]=cityscape_hist(n,x);
figure();
plot(x,n./nansum(n),'Color','k');
xlabel('Change in RT');
ylabel('Count');
title(tit);

if ~isempty(testcond2)
    [n,x]=histcounts(RT_pairs.rt_change(RT_pairs.real_rt_pair==true & testcond2==true),x_backup);
    [n,x]=cityscape_hist(n,x);
    hold on;
    plot(x,n./nansum(n),'Color','r');
    leg={'testcond 1','testcond 2'};
    legend(leg);
end

end