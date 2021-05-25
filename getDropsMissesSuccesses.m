function [successes,drops,misses,totalTrials]=getDropsMissesSuccesses(trialTypes,metadata,useMouseID)
disp('calculating successes, drops, etc. AT ANY TIME IN TRIAL');
useSess=unique(metadata.nth_session(metadata.mouseid==useMouseID));
successes=nan(1,length(useSess));
drops=nan(1,length(useSess));
misses=nan(1,length(useSess));
totalTrials=nan(1,length(useSess));
for i=1:length(useSess)
    successes(i)=nansum(trialTypes.consumed_pellet(metadata.mouseid==useMouseID & metadata.nth_session==useSess(i))==1);
end
for i=1:length(useSess)
    drops(i)=nansum(trialTypes.touched_pellet(metadata.mouseid==useMouseID & metadata.nth_session==useSess(i))==1 & trialTypes.consumed_pellet(metadata.mouseid==useMouseID & metadata.nth_session==useSess(i))==0);
end
for i=1:length(useSess)
    misses(i)=nansum(trialTypes.touched_pellet(metadata.mouseid==useMouseID & metadata.nth_session==useSess(i))==0);
end
for i=1:length(useSess)
    totalTrials(i)=nansum(metadata.mouseid==useMouseID & metadata.nth_session==useSess(i));
end