function plotNthReactionTimeAway(templateSequence,metadata,trialTypes,tbt,whichReach,useAsCue,nNexts)

% have templateSequence take on a form like
% templateSequence{1}=trialsMatchingCondA;
% templateSequence{2}=nan;      nan is a wildcard indicating that 2nd trial
% in sequence can have any value
% templateSequence{3}=trialsMatchingCondB;
% templateSequence{4}=trialsMatchingCondC;

doPlot=1;

f=fieldnames(trialTypes);
% get matches to templateSequence
sequencesInSet=zeros(1,length(trialTypes.(f{1}))-(length(templateSequence)-1));
for i=1:length(templateSequence)
    if ~isnan(templateSequence{i})
        temp=templateSequence{i}==1;
        % make row vector
        if size(temp,1)>1
            temp=temp';
        end
        if size(temp,1)>1
            error('templateSequence{i} must be 1D vector');
        end
        sequencesInSet(1:length(trialTypes.(f{1}))-(length(templateSequence)-1))=sequencesInSet(1:length(trialTypes.(f{1}))-(length(templateSequence)-1))+temp(i:end-((length(templateSequence)-1)-(i-1)));
    else
        sequencesInSet=sequencesInSet+1; % nan in templateSequence indicates wildcard, all trial types acceptable
    end
end
% look for spans of trials where all conditions in templateSequence are
% true sequentially
sequenceMatchStarts=sequencesInSet>=length(templateSequence);
forPairs_sequenceMatchStarts=sequenceMatchStarts;
sequenceMatchStarts=[sequenceMatchStarts zeros(1,length(templateSequence)-1)];

% Iterate through different nNexts (i spacings between trials to get
% reaction time pairs)
currColor=[0    0.1; ...
           0    0.6; ...
           0    0.7];
% currColor=[0    0.6; ...
%            0    0.2; ...
%            0    0.8];
% currColor=[0    0.5     0.9     0.9     0.9     0.6     0.2     0.1     0.1     0.1     0.3; ...
%            0    0       0       0.5     0.8     0.9     0.5     0.6     0.3     0       0; ...
%            0    0       0       0       0       0       0       0.7     0.9     0.7     0.5];
RTpairs=cell(1,length(nNexts));
for i=1:length(nNexts)
    currNext=nNexts(i);
    if i==1
        hist_h=[];
        cdf_h=[];
        scatter_h=[];
    end
    [~,RT_pairs,hist_h,cdf_h,scatter_h]=getRTpairs(tbt,whichReach,useAsCue,metadata,sequenceMatchStarts,[],1,[],num2str(currNext),[],currNext,doPlot,hist_h,cdf_h,scatter_h,currColor(:,i));
    RTpairs{i}=RT_pairs.rt_change(RT_pairs.real_rt_pair==true);
end
p=ranksum(RTpairs{1},RTpairs{2});
disp('ranksum p-val');
disp(p);

end

function [reactionTimes,RT_pairs,tbt]=getPairedRTs(tbt,whichReach,useAsCue,metadata,zscore_RTs,selectTrials,nApart)

longRT_ifNoReach=0;

% cue ind
avcue=nanmean(tbt.(useAsCue),1);
[~,maind]=max(avcue);

% which reaches to get
temp=tbt.(whichReach);
newtemp=nan(size(tbt.(whichReach)));

% get only first reaches
firstreaches=nan(1,size(temp,1));
% note that everything has already been cue-aligned
for i=1:size(temp,1)
    fi=find(temp(i,maind+1:end)>0,1,'first')+maind;
    if ~isempty(fi)
        fi=fi(1);
        firstreaches(i)=fi;
    end
    if longRT_ifNoReach==1 && isempty(fi)
        fi=size(temp,2);
        firstreaches(i)=fi;
    end
    if ~isempty(fi)
        newtemp(i,:)=zeros(size(newtemp(i,:)));
        newtemp(i,fi)=1;
    end
end
reactionTimes=(firstreaches-maind).*mode(diff(nanmean(tbt.times,1)));
tbt.firstReachesAfterCueOnset=newtemp;

if zscore_RTs==1
    reactionTimes=Zscore_by_session(reactionTimes,metadata.sessid);
end

if ~isempty(selectTrials)
    if isempty(nApart) 
        RT_pairs=getRTchange_trialToTrial_varyApart(reactionTimes,metadata,selectTrials);
    else
        RT_pairs=getRTchange_trialToTrial_varyAnd_nApart(reactionTimes,metadata,selectTrials,nApart);
    end
else
    if isempty(nApart)
        RT_pairs=getRTchange_trialToTrial(reactionTimes,metadata);
    else
        RT_pairs=getRTchange_trialToTrial_nApart(reactionTimes,metadata,nApart);
    end
end

end

function out=getRTchange_trialToTrial(curr_rt,metadata)

% set up output
out.rt_change=nan(1,length(curr_rt)-1);
out.prev_trial_num=1:length(curr_rt)-1;
out.prev_trial_rt=nan(1,length(curr_rt)-1);
out.curr_trial_num=2:length(curr_rt);
out.curr_trial_rt=nan(1,length(curr_rt)-1);
out.real_rt_pair=nan(1,length(curr_rt)-1);
out.is_video_file_transition=zeros(1,length(curr_rt));

% find transitions between video files
% this is important, because don't want to consider these trials
% as "neighbors", i.e., previous and current trials
temp=metadata.sessid(2:end)'-metadata.sessid(1:end-1)';
temp(temp~=0)=1;
out.real_rt_pair=temp==0;
out.is_video_file_transition=[temp 0];
f=find(out.is_video_file_transition>0);
out.is_video_file_transition(f+1)=1;

% check whether empty reaction times
if all(isnan(curr_rt(1:end-1)-curr_rt(2:end)))
    return
end

out.rt_change=curr_rt(1:end-1)-curr_rt(2:end); % previous minus current
out.prev_trial_rt=curr_rt(1:end-1);
out.curr_trial_rt=curr_rt(2:end);
out.rt_change(out.real_rt_pair==false)=nan;

end

function out=getRTchange_trialToTrial_nApart(curr_rt,metadata,nApart)

% set up output
out.rt_change=nan(1,length(curr_rt)-nApart);
out.prev_trial_num=1:length(curr_rt)-nApart;
out.prev_trial_rt=nan(1,length(curr_rt)-nApart);
out.curr_trial_num=nApart+1:length(curr_rt);
out.curr_trial_rt=nan(1,length(curr_rt)-nApart);
out.real_rt_pair=nan(1,length(curr_rt)-nApart);
out.is_video_file_transition=zeros(1,length(curr_rt));

% find transitions between video files
% this is important, because don't want to consider these trials
% as "neighbors", i.e., previous and current trials
temp=metadata.sessid(2:end)'-metadata.sessid(1:end-1)';
temp(temp~=0)=1;
out.is_video_file_transition=[temp 0];
f=find(out.is_video_file_transition>0);
out.is_video_file_transition(f+1)=1;
atFileTransition=zeros(1,length(out.prev_trial_num));
for i=1:length(out.prev_trial_num)
    % check for file break between previous and current trials
    if sum(out.is_video_file_transition(out.prev_trial_num(i):out.curr_trial_num(i)))>1
        atFileTransition(i)=1;
    end
end
out.real_rt_pair=atFileTransition==0;

% check whether empty reaction times
if all(isnan(curr_rt(1:end-1)-curr_rt(2:end)))
    return
end

out.rt_change=curr_rt(out.prev_trial_num)-curr_rt(out.curr_trial_num); % previous minus current
out.prev_trial_rt=curr_rt(out.prev_trial_num);
out.curr_trial_rt=curr_rt(out.curr_trial_num);
out.rt_change(out.real_rt_pair==false)=nan;

end

function out=getRTchange_trialToTrial_varyApart(curr_rt,metadata,selectTrials)

select_inds=find(selectTrials==1);

% make select_inds a row vector
if size(select_inds,1)>1
    select_inds=select_inds';
end
if size(select_inds,1)>1
    error('select trials need to be a 1D vector');
end

% set up output
out.rt_change=nan(1,length(select_inds)-1);
out.prev_trial_num=select_inds(1:end-1);
out.prev_trial_rt=nan(1,length(select_inds)-1);
out.curr_trial_num=select_inds(2:end);
out.curr_trial_rt=nan(1,length(select_inds)-1);
out.real_rt_pair=nan(1,length(select_inds)-1);
out.is_video_file_transition=zeros(1,length(curr_rt));

% find transitions between video files
% this is important, because don't want to consider RT changes
% across file breaks
temp=metadata.sessid(2:end)'-metadata.sessid(1:end-1)';
temp(temp~=0)=1;
out.is_video_file_transition=[temp 0];
f=find(out.is_video_file_transition>0);
out.is_video_file_transition(f+1)=1;
atFileTransition=zeros(1,length(out.prev_trial_num));
for i=1:length(out.prev_trial_num)
    % check for file break between earlier and current trials
    if sum(out.is_video_file_transition(out.prev_trial_num(i):out.curr_trial_num(i)))>1
        atFileTransition(i)=1;
    end
end
out.real_rt_pair=atFileTransition==0;

% check whether empty reaction times
if all(isnan(curr_rt(1:end-1)-curr_rt(2:end)))
    return
end

out.rt_change=curr_rt(out.prev_trial_num)-curr_rt(out.curr_trial_num); % previous minus current
out.prev_trial_rt=curr_rt(out.prev_trial_num);
out.curr_trial_rt=curr_rt(out.curr_trial_num);
out.rt_change(out.real_rt_pair==false)=nan;

end

function out=getRTchange_trialToTrial_varyAnd_nApart(curr_rt,metadata,selectTrials,nApart)

select_inds=find(selectTrials==1);

% make select_inds a row vector
if size(select_inds,1)>1
    select_inds=select_inds';
end
if size(select_inds,1)>1
    error('select trials need to be a 1D vector');
end

% set up output
out.rt_change=nan(1,length(select_inds)-nApart);
out.prev_trial_num=select_inds(1:end-nApart);
out.prev_trial_rt=nan(1,length(select_inds)-nApart);
out.curr_trial_num=select_inds(nApart+1:end);
out.curr_trial_rt=nan(1,length(select_inds)-nApart);
out.real_rt_pair=nan(1,length(select_inds)-nApart);
out.is_video_file_transition=zeros(1,length(curr_rt));

% find transitions between video files
% this is important, because don't want to consider RT changes
% across file breaks
temp=metadata.sessid(2:end)'-metadata.sessid(1:end-1)';
temp(temp~=0)=1;
out.is_video_file_transition=[temp 0];
f=find(out.is_video_file_transition>0);
out.is_video_file_transition(f+1)=1;
atFileTransition=zeros(1,length(out.prev_trial_num));
for i=1:length(out.prev_trial_num)
    % check for file break between earlier and current trials
    if sum(out.is_video_file_transition(out.prev_trial_num(i):out.curr_trial_num(i)))>1
        atFileTransition(i)=1;
    end
end
out.real_rt_pair=atFileTransition==0;

% check whether empty reaction times
if all(isnan(curr_rt(1:end-1)-curr_rt(2:end)))
    return
end

out.rt_change=curr_rt(out.prev_trial_num)-curr_rt(out.curr_trial_num); % previous minus current
out.prev_trial_rt=curr_rt(out.prev_trial_num);
out.curr_trial_rt=curr_rt(out.curr_trial_num);
out.rt_change(out.real_rt_pair==false)=nan;

end


function [reactionTimes,RT_pairs,hist_h,cdf_h,scatter_h]=getRTpairs(tbt,whichReach,useAsCue,metadata,contingency1,contingency2,prevOrCurrTrial1,prevOrCurrTrial2,contingencyName,selectTrials,nApart,doPlot,hist_h,cdf_h,scatter_h,currColor)

% if contingency2 is empty, will not consider this contingency
% will just plot contingency1==true trials

nbins=15; % for histograms
jitter=0.01; % for scatter plots
addSlopeLine=0; % for scatter plots, will add fit to slope if 1
zscore_RTs=0; % if 1, will use z-scored instead of raw reaction times

[reactionTimes,RT_pairs,tbt]=getPairedRTs(tbt,whichReach,useAsCue,metadata,zscore_RTs,selectTrials,nApart);
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
    hist_h=plotHist(hist_h,RT_pairs,nbins,['Compare RT change from prev to curr trial dependent on: ' contingencyName],contingency1,contingency2,currColor);
    cdf_h=plotCDF(cdf_h,RT_pairs,4*nbins,['CDF ' contingencyName],contingency1,contingency2,currColor);
    scatter_h=plotScatter(scatter_h,RT_pairs,['Prev vs curr trial RT dependent on: ' contingencyName],contingency1,contingency2,jitter,addSlopeLine,currColor);
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

function h=plotScatter(h,RT_pairs,tit,testcond1,testcond2,jitter,addSlopeLine,currColor)

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

if isempty(h)
    h=figure();
else
    figure(h);
end
usePrev=RT_pairs.prev_trial_rt+rand(size(RT_pairs.prev_trial_rt)).*jitter;
useCurr=RT_pairs.curr_trial_rt+rand(size(RT_pairs.curr_trial_rt)).*jitter;
s=scatter(usePrev(RT_pairs.real_rt_pair==true & testcond1==true),useCurr(RT_pairs.real_rt_pair==true & testcond1==true),[],'k','filled');
pause;
m=get(s,'MarkerHandle');
alpha=0.3;
m.FaceColorData=uint8(255*[currColor; alpha]);
xlabel('RT previous trial in sec');
ylabel('RT current trial in sec');
title(tit);
hold on;
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

function h=plotCDF(h,RT_pairs,bins,tit,testcond1,testcond2,currColor)

persistent x_backup

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

if isempty(h)
    [n,x]=histcounts(RT_pairs.rt_change(RT_pairs.real_rt_pair==true & testcond1==true),bins);
else
    [n,x]=histcounts(RT_pairs.rt_change(RT_pairs.real_rt_pair==true & testcond1==true),x_backup);
end
x_mids=nanmean([x(1:end-1); x(2:end)],1);
cond1_cdf=accumulateDistribution(n);
if isempty(h)
    h=figure();
    x_backup=x;
else
    figure(h);
end
plot(x_mids,cond1_cdf./nanmax(cond1_cdf),'Color',currColor);
xlabel('CDF of change in RT');
ylabel('Count');
title(tit);
hold on; 

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

function h=plotHist(h,RT_pairs,bins,tit,testcond1,testcond2,currColor)

persistent x_backup

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

if isempty(h)
    [n,x]=histcounts(RT_pairs.rt_change(RT_pairs.real_rt_pair==true & testcond1==true),bins);
else
    [n,x]=histcounts(RT_pairs.rt_change(RT_pairs.real_rt_pair==true & testcond1==true & ~isnan(RT_pairs.rt_change)),x_backup);
end
if isempty(h)
    h=figure();
    x_backup=x;
else
    figure(h);
end
[n,x]=cityscape_hist(n,x);
plot(x,n./nansum(n),'Color',currColor);
xlabel('Change in RT');
ylabel('Count');
title(tit);
hold on;

if ~isempty(testcond2)
    [n,x]=histcounts(RT_pairs.rt_change(RT_pairs.real_rt_pair==true & testcond2==true),x_backup);
    [n,x]=cityscape_hist(n,x);
    hold on;
    plot(x,n./nansum(n),'Color','r');
    leg={'testcond 1','testcond 2'};
    legend(leg);
end

end
