function [reactionTimes,RT_pairs,tbt]=getPairedReactionTimes(tbt,whichReach,useAsCue,metadata,zscore_RTs,selectTrials,nApart)

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
