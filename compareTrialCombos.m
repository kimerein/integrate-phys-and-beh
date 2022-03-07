function [p_RT_pairs,p_RTs,out]=compareTrialCombos(varargin)

% 3 formats for arguments to this function
% Format 1 (3 args): tbt, trialTypes, metadata
% Format 2 (5 args): tbt, trialTypes, metadata, templateSequence1,
%   templateSequence2
% Format 3 (10 args); tbt, trialTypes, metadata, templateSequence1,
%   nNext1, templateSequence2, nNext2, useAsCue, whichReach, doPlot

if length(varargin)==3
    tbt=varargin{1};
    trialTypes=varargin{2};
    metadata=varargin{3};
    % settings
    whichReach='all_reachBatch';
    useAsCue='cueZone_onVoff';
    doPlot=1;
    nNext1=1;
    nNext2=1;
elseif length(varargin)==5
    tbt=varargin{1};
    trialTypes=varargin{2};
    metadata=varargin{3};
    templateSequence1=varargin{4};
    templateSequence2=varargin{5};
    % settings
    whichReach='all_reachBatch';
    useAsCue='cueZone_onVoff';
    doPlot=0;
    nNext1=1;
    nNext2=1;
elseif length(varargin)==10
    tbt=varargin{1};
    trialTypes=varargin{2};
    metadata=varargin{3};
    templateSequence1=varargin{4};
    nNext1=varargin{5};
    templateSequence2=varargin{6};
    nNext2=varargin{7};
    useAsCue=varargin{8};
    whichReach=varargin{9};
    doPlot=varargin{10};
else
    error('wrong number of arguments to compareTrialCombos');
end

% trial type 1
templateSequence=templateSequence1;
nNext=nNext1;
[reactionTimes1,sequenceMatchStarts1,forPairs_sequenceMatchStarts1,RT_pairs1]=plotTrialCombosAndRTs(templateSequence,metadata,trialTypes,tbt,whichReach,useAsCue,nNext,0);

% out.rt_change=nan(1,length(curr_rt)-1);
% out.prev_trial_num=1:length(curr_rt)-1;
% out.prev_trial_rt=nan(1,length(curr_rt)-1);
% out.curr_trial_num=2:length(curr_rt);
% out.curr_trial_rt=nan(1,length(curr_rt)-1);
% out.real_rt_pair=nan(1,length(curr_rt)-1);
% out.is_video_file_transition=zeros(1,length(curr_rt));


% trial type 2
templateSequence=templateSequence2;
nNext=nNext2;
[reactionTimes2,sequenceMatchStarts2,forPairs_sequenceMatchStarts2,RT_pairs2]=plotTrialCombosAndRTs(templateSequence,metadata,trialTypes,tbt,whichReach,useAsCue,nNext,0);

if length(forPairs_sequenceMatchStarts1)<length(RT_pairs1.real_rt_pair)
    RT_pairs1.real_rt_pair=RT_pairs1.real_rt_pair(1:length(forPairs_sequenceMatchStarts1));
end
if length(forPairs_sequenceMatchStarts2)<length(RT_pairs2.real_rt_pair)
    RT_pairs2.real_rt_pair=RT_pairs2.real_rt_pair(1:length(forPairs_sequenceMatchStarts2));
end

% save output
out.rt_pairs1_contingent=RT_pairs1.rt_change(RT_pairs1.real_rt_pair==true & forPairs_sequenceMatchStarts1==true);
out.rt_pairs2_contingent=RT_pairs2.rt_change(RT_pairs2.real_rt_pair==true & forPairs_sequenceMatchStarts2==true);
out.sequenceMatchStarts1=sequenceMatchStarts1;
out.sequenceMatchStarts2=sequenceMatchStarts2;
out.all_rt1=reactionTimes1;
out.all_rt2=reactionTimes2;
out.rt_pairs1=RT_pairs1.rt_change;
out.rt_pairs2=RT_pairs2.rt_change;
out.real_rt_pair1=RT_pairs1.real_rt_pair;
out.real_rt_pair2=RT_pairs2.real_rt_pair;
out.forPairs_sequenceMatchStarts1=forPairs_sequenceMatchStarts1;
out.forPairs_sequenceMatchStarts2=forPairs_sequenceMatchStarts2;
out.all_rt1_triali=[reactionTimes1(1+nNext1:end) nan(1,nNext1)];
out.all_rt2_triali=[reactionTimes2(1+nNext2:end) nan(1,nNext2)];

if doPlot==1
    % reaction time pairs comparison
%     plotHist(RT_pairs1,RT_pairs2,forPairs_sequenceMatchStarts1,forPairs_sequenceMatchStarts2,250,'Compare Pairs','Change in RT');
%     plotCDF(RT_pairs1,RT_pairs2,forPairs_sequenceMatchStarts1,forPairs_sequenceMatchStarts2,4*200,'CDF Compare Pairs');
%     plotScatter(RT_pairs1,RT_pairs2,forPairs_sequenceMatchStarts1,forPairs_sequenceMatchStarts2,0.03,'Scatter Compare Pairs');
    plotHist(RT_pairs1,RT_pairs2,forPairs_sequenceMatchStarts1,forPairs_sequenceMatchStarts2,200,'Compare Pairs','Change in RT');
    plotCDF(RT_pairs1,RT_pairs2,forPairs_sequenceMatchStarts1,forPairs_sequenceMatchStarts2,4*200,'CDF Compare Pairs');
    plotScatter(RT_pairs1,RT_pairs2,forPairs_sequenceMatchStarts1,forPairs_sequenceMatchStarts2,0.03,'Scatter Compare Pairs');
end
p_RT_pairs=testRanksum(RT_pairs1,RT_pairs2,forPairs_sequenceMatchStarts1,forPairs_sequenceMatchStarts2,doPlot);
    
% reaction times comparison
temp1.rt_change=reactionTimes1;
temp1.real_rt_pair=ones(size(reactionTimes1));
temp2.rt_change=reactionTimes2;
temp2.real_rt_pair=ones(size(reactionTimes2));
f1=find(sequenceMatchStarts1==1);
f2=find(sequenceMatchStarts2==1);
temp_match1=zeros(size(sequenceMatchStarts1));
temp_match1(f1+nNext)=1;
temp_match2=zeros(size(sequenceMatchStarts2));
temp_match2(f2+nNext)=1;
if doPlot==1
    plotHist(temp1,temp2,temp_match1,temp_match2,250,'Compare Non-Paired RTs','RTs');
end
p_RTs=testRanksum(temp1,temp2,temp_match1,temp_match2,doPlot);

end

function p=testRanksum(RT_pairs1,RT_pairs2,testcond1,testcond2,dispStuff)

if dispStuff==1
    disp('median of contingency 1');
    disp(nanmedian(RT_pairs1.rt_change(RT_pairs1.real_rt_pair==true & testcond1==true),2));
    disp('median of contingency 2');
    disp(nanmedian(RT_pairs2.rt_change(RT_pairs2.real_rt_pair==true & testcond2==true),2));
end

if isempty(RT_pairs1.real_rt_pair)
    p=nan;
    return
end

if all(isnan(RT_pairs1.rt_change(RT_pairs1.real_rt_pair==true & testcond1==true))) | all(isnan(RT_pairs2.rt_change(RT_pairs2.real_rt_pair==true & testcond2==true)))
    p=nan;
    return
end
[p,h]=ranksum(RT_pairs1.rt_change(RT_pairs1.real_rt_pair==true & testcond1==true),RT_pairs2.rt_change(RT_pairs2.real_rt_pair==true & testcond2==true));
if dispStuff==1
    disp('p-value');
    disp(p);
end

end

function plotScatter(RT_pairs1,RT_pairs2,testcond1,testcond2,jitter,tit)

figure();
usePrev=RT_pairs1.prev_trial_rt+rand(size(RT_pairs1.prev_trial_rt)).*jitter;
useCurr=RT_pairs1.curr_trial_rt+rand(size(RT_pairs1.curr_trial_rt)).*jitter;

s=scatter(usePrev(RT_pairs1.real_rt_pair==true & testcond1==true),useCurr(RT_pairs1.real_rt_pair==true & testcond1==true),50,'k','filled');
set(gcf,'position',[10,10,1000,1000]);
pause;
m=get(s,'MarkerHandle');
alpha=0.25;
m.FaceColorData=uint8(255*[0;0;0;alpha]);
xlabel('RT previous trial');
ylabel('RT current trial');
title(tit);
% xlim([0 5]);
% ylim([0 9.5]);

hold on;
usePrev=RT_pairs2.prev_trial_rt+rand(size(RT_pairs2.prev_trial_rt)).*jitter;
useCurr=RT_pairs2.curr_trial_rt+rand(size(RT_pairs2.curr_trial_rt)).*jitter;
s=scatter(usePrev(RT_pairs2.real_rt_pair==true & testcond2==true),useCurr(RT_pairs2.real_rt_pair==true & testcond2==true),50,'r','filled');
pause;
m=get(s,'MarkerHandle');
alpha=0.25;
m.FaceColorData=uint8(255*[1;0;0;alpha]);
leg={'testcond 1','testcond 2'};
legend(leg);

end

function plotHist(RT_pairs1,RT_pairs2,testcond1,testcond2,bins,tit,xlab)

[n,x]=histcounts(RT_pairs1.rt_change(RT_pairs1.real_rt_pair==true & testcond1==true),bins);
x_backup=x;
[n,x]=cityscape_hist(n,x);
figure();
plot(x,n./nansum(n),'Color','k');
xlabel(xlab);
ylabel('Count');
title(tit);

[n,x]=histcounts(RT_pairs2.rt_change(RT_pairs2.real_rt_pair==true & testcond2==true),x_backup);
[n,x]=cityscape_hist(n,x);
hold on;
plot(x,n./nansum(n),'Color','r');
leg={'cond 1','cond 2'};
legend(leg);

end

function plotCDF(RT_pairs1,RT_pairs2,testcond1,testcond2,bins,tit)

[n,x]=histcounts(RT_pairs1.rt_change(RT_pairs1.real_rt_pair==true & testcond1==true),bins);
x_backup=x;
x_mids=nanmean([x(1:end-1); x(2:end)],1);
cond1_cdf=accumulateDistribution(n);
figure();
plot(x_mids,cond1_cdf./nanmax(cond1_cdf),'Color','k');
xlabel('CDF of change in RT');
ylabel('Count');
title(tit);

[n,x]=histcounts(RT_pairs2.rt_change(RT_pairs2.real_rt_pair==true & testcond2==true),x_backup);
hold on;
cond2_cdf=accumulateDistribution(n);
plot(x_mids,cond2_cdf./nanmax(cond2_cdf),'Color','r');
leg={'testcond 1','testcond 2'};
legend(leg);

end

function cdf=accumulateDistribution(data)

cdf=nan(size(data));
for i=1:length(data)
    cdf(i)=nansum(data(1:i));
end

end