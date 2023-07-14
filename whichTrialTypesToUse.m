function [trial1,trial2,trial1_LED,trial2_LED]=whichTrialTypesToUse(alltbt,trialTypes,metadata,whichEventType,timeWindow)

[~,cueindma]=nanmax(nanmean(alltbt.cueZone_onVoff,1));
timestep=mode(diff(nanmean(alltbt.times,1)));
if isempty(timeWindow)
    trialTypes.reachedInTimeWindow=ones(size(trialTypes.led));
    trialTypes.reachedInTimeWindow_1forward=[trialTypes.reachedInTimeWindow(2:end); 1];
    trialTypes.reachedInTimeWindow_1back=[1; trialTypes.reachedInTimeWindow(1:end-1)];
else
    timeWindowInds(1)=floor(timeWindow(1)/timestep);
    timeWindowInds(2)=floor(timeWindow(2)/timestep);
    trialTypes.reachedInTimeWindow=any(alltbt.all_reachBatch(:,cueindma+timeWindowInds(1)-1:cueindma+timeWindowInds(2))>0.05,2);
    trialTypes.reachedInTimeWindow_1forward=[trialTypes.reachedInTimeWindow(2:end); 0];
    trialTypes.reachedInTimeWindow_1back=[0; trialTypes.reachedInTimeWindow(1:end-1)];
end

switch whichEventType
    case 'success'
        trial1='trialTypes.consumed_pellet_1back==1 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.success_in_cued_window_1forward==1 & trialTypes.consumed_pellet_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==0'; 
        trial2='trialTypes.optoGroup~=1 & (trialTypes.led_1forward==1 | trialTypes.led_2forward==1 | trialTypes.led_3forward==1 | trialTypes.led_4forward==1 | trialTypes.led_1back==1)';
        trial1_LED='trialTypes.consumed_pellet_1back==1 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.success_in_cued_window_1forward==1 & trialTypes.consumed_pellet_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==1'; 
        trial2_LED='trialTypes.optoGroup~=1 & (trialTypes.led_1forward==0 | trialTypes.led_2forward==0 | trialTypes.led_3forward==0 | trialTypes.led_4forward==0 | trialTypes.led_1back==0)';
end

end