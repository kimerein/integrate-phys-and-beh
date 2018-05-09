function [unitFx,behFx,supp]=getBehWithUnitEffects(spikes,beh_sweeps,led_window,norm_window,tbt,whichReach,binsize)

temp=unique(spikes.sweeps.trials);
bothSpikesAndBehTrials=temp(ismember(temp,unique(spikes.trials)));
bothLogical=ismember(1:length(beh_sweeps.led),bothSpikesAndBehTrials);

% choose behavior condition to analyze
beh_cond=beh_sweeps.pawOutDuringWheel==0 & beh_sweeps.reachedAfterCue==1 & bothLogical';
beh_cond=beh_cond';

% get spikes with and without opto for this behavior condition
[unitFx,all_frs,supp]=getUnitsEffects_withBehavior(spikes,led_window,norm_window,beh_cond,binsize);

% get reaches for this behavior condition
behFx.x=nanmean(tbt.times(beh_cond,:),1);
behFx.cue=nanmean(tbt.cueZone_onVoff(beh_cond,:),1);
behFx.opto=nanmean(tbt.optoZone(beh_cond,:),1);
temp=tbt.(whichReach);
behFx.reaches_no_led=nanmean(temp(beh_cond & beh_sweeps.led==0,:),1);
behFx.reaches_led=nanmean(temp(beh_cond & beh_sweeps.led==1,:),1);

