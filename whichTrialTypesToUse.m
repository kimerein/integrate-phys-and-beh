function [trialTypes,trial1,trial2,trial1_LED,trial2_LED]=whichTrialTypesToUse(alltbt,trialTypes,metadata,whichEventType,timeWindow,whichReachInTimeWindow)

[~,cueindma]=nanmax(nanmean(alltbt.cueZone_onVoff,1));
timestep=mode(diff(nanmean(alltbt.times,1)));
if isempty(timeWindow)
    trialTypes.reachedInTimeWindow=ones(size(trialTypes.led));
    trialTypes.reachedInTimeWindow_1forward=[trialTypes.reachedInTimeWindow(2:end); 1];
    trialTypes.reachedInTimeWindow_1back=[1; trialTypes.reachedInTimeWindow(1:end-1)];
else
    timeWindowInds(1)=floor(timeWindow(1)/timestep);
    timeWindowInds(2)=floor(timeWindow(2)/timestep);
    temp=alltbt.(whichReachInTimeWindow);
    trialTypes.reachedInTimeWindow=any(temp(:,cueindma+timeWindowInds(1):cueindma+timeWindowInds(2))>0.05,2);
    trialTypes.reachedInTimeWindow_1forward=[trialTypes.reachedInTimeWindow(2:end); 0];
    trialTypes.reachedInTimeWindow_1back=[0; trialTypes.reachedInTimeWindow(1:end-1)];
end

trialTypes.isLongITI_1forward=[trialTypes.isLongITI(2:end); 0];
trialTypes.optoGroup_1forward=[trialTypes.optoGroup(2:end); 0];
trialTypes.optoGroup_1back=[0; trialTypes.optoGroup(1:end-1)];
trialTypes.isLongITI_1back=[0; trialTypes.isLongITI(1:end-1)];
trialTypes.optoGroup_2back=[0; 0; trialTypes.optoGroup(1:end-2)];
trialTypes.noReach=~any(alltbt.all_reachBatch>0.05,2);
trialTypes.noReach_1forward=[trialTypes.noReach(2:end); 0];
trialTypes.noReach_1back=[0; trialTypes.noReach(1:end-1)];
trialTypes.reachedBeforeCue=any(alltbt.all_reachBatch(:,1:cueindma-1)>0.05,2);
trialTypes.reachedAfterCue=any(alltbt.all_reachBatch(:,cueindma:end)>0.05,2);
trialTypes.reachToPelletBeforeCue=any(alltbt.reachStarts_pelletPresent(:,1:cueindma-1)>0.05,2);
trialTypes.reachedBeforeCue_1forward=[trialTypes.reachedBeforeCue(2:end); 0];
trialTypes.reachedBeforeCue_1back=[0; trialTypes.reachedBeforeCue(1:end-1)];
trialTypes.reachToPelletBeforeCue_1forward=[trialTypes.reachToPelletBeforeCue(2:end); 0];
trialTypes.reachToPelletBeforeCue_1back=[0; trialTypes.reachToPelletBeforeCue(1:end-1)];
trialTypes.reachedAfterCue_1forward=[trialTypes.reachedAfterCue(2:end); 0];

switch whichEventType
    case 'cued success'
        trial1=['(trialTypes.optoGroup~=1 & trialTypes.consumed_pellet_1back==1 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==0) | ' ...
                '(trialTypes.optoGroup~=1 & trialTypes.isLongITI_1forward==1 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==0)']; 
        trial2='trialTypes.optoGroup~=1 & (trialTypes.led_1forward==1 | trialTypes.led_2forward==1 | trialTypes.led_3forward==1 | trialTypes.led_4forward==1 | trialTypes.led_1back==1)';
        trial1_LED=['(trialTypes.optoGroup~=1 & trialTypes.consumed_pellet_1back==1 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==1) | ' ...
                    '(trialTypes.optoGroup~=1 & trialTypes.isLongITI_1forward==1 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==1)']; 
        trial2_LED='trialTypes.optoGroup~=1 & (trialTypes.led_1forward==0 | trialTypes.led_2forward==0 | trialTypes.led_3forward==0 | trialTypes.led_4forward==0 | trialTypes.led_1back==0)';
    case 'cued success accumulate'
        trial1=['(trialTypes.optoGroup~=1 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==0)'];
        trial2='trialTypes.optoGroup~=1 & (trialTypes.led_1forward==1 | trialTypes.led_2forward==1 | trialTypes.led_3forward==1 | trialTypes.led_4forward==1 | trialTypes.led_1back==1)';
        trial1_LED=['(trialTypes.optoGroup~=1 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==1)'];
        trial2_LED='trialTypes.optoGroup~=1 & (trialTypes.led_1forward==0 | trialTypes.led_2forward==0 | trialTypes.led_3forward==0 | trialTypes.led_4forward==0 | trialTypes.led_1back==0)';
    case 'cued failure'
        % did not touch pellet despite reaching in cued window
        trial1=['(trialTypes.optoGroup~=1 & trialTypes.cued_reach_1back==1 & trialTypes.isLongITI_1back==1 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touch_in_cued_window_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==0) | ' ...
                '(trialTypes.optoGroup~=1 & trialTypes.isLongITI_1forward==1                               & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touch_in_cued_window_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==0)']; 
        trial2='trialTypes.optoGroup~=1 & (trialTypes.led_1forward==1 | trialTypes.led_2forward==1 | trialTypes.led_3forward==1 | trialTypes.led_4forward==1 | trialTypes.led_1back==1)';
        trial1_LED=['(trialTypes.optoGroup~=1 & trialTypes.cued_reach_1back==1 & trialTypes.isLongITI_1back==1 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touch_in_cued_window_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==1) | ' ...
                    '(trialTypes.optoGroup~=1 & trialTypes.isLongITI_1forward==1                               & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touch_in_cued_window_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==1)']; 
        trial2_LED='trialTypes.optoGroup~=1 & (trialTypes.led_1forward==0 | trialTypes.led_2forward==0 | trialTypes.led_3forward==0 | trialTypes.led_4forward==0 | trialTypes.led_1back==0)';
%         trial1=['(trialTypes.optoGroup~=1 & trialTypes.touch_in_cued_window_1back==0 & trialTypes.noReach_1back==0 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touch_in_cued_window_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==0) | ' ...
%                 '(trialTypes.optoGroup~=1 & trialTypes.isLongITI_1forward==1                                       & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touch_in_cued_window_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==0)']; 
%         trial2='trialTypes.optoGroup~=1 & (trialTypes.led_1forward==1 | trialTypes.led_2forward==1 | trialTypes.led_3forward==1 | trialTypes.led_4forward==1 | trialTypes.led_1back==1)';
%         trial1_LED=['(trialTypes.optoGroup~=1 & trialTypes.touch_in_cued_window_1back==0 & trialTypes.noReach_1back==0 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touch_in_cued_window_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==1) | ' ...
%                     '(trialTypes.optoGroup~=1 & trialTypes.isLongITI_1forward==1                                       & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touch_in_cued_window_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==1)']; 
%         trial2_LED='trialTypes.optoGroup~=1 & (trialTypes.led_1forward==0 | trialTypes.led_2forward==0 | trialTypes.led_3forward==0 | trialTypes.led_4forward==0 | trialTypes.led_1back==0)';
    case 'delayed success'
        trial1=['(trialTypes.optoGroup~=1 & trialTypes.consumed_pellet_1back==1 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==0) | ' ...
                '(trialTypes.optoGroup~=1 & trialTypes.isLongITI_1forward==1    & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==0)']; 
        trial2='trialTypes.optoGroup~=1 & (trialTypes.led_1forward==1 | trialTypes.led_2forward==1 | trialTypes.led_3forward==1 | trialTypes.led_4forward==1 | trialTypes.led_1back==1)';
        trial1_LED=['(trialTypes.optoGroup~=1 & trialTypes.consumed_pellet_1back==1 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==1) | ' ...
                    '(trialTypes.optoGroup~=1 & trialTypes.isLongITI_1forward==1    & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==1)']; 
        trial2_LED='trialTypes.optoGroup~=1 & (trialTypes.led_1forward==0 | trialTypes.led_2forward==0 | trialTypes.led_3forward==0 | trialTypes.led_4forward==0 | trialTypes.led_1back==0)';
        % could put chewing_at_trial_start==0 for trial2 and trial2_LED but
        % it does NOT make a difference and reduces trials
    case 'uncued failure'
        % did not touch pellet despite reaching in delayed window after cue
        trial1=['(trialTypes.optoGroup~=1 & trialTypes.consumed_pellet_1back==0 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touched_pellet_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==0) |' ...
                '(trialTypes.optoGroup~=1 & trialTypes.consumed_pellet_1back==0 & trialTypes.reachedBeforeCue_1forward==1 & trialTypes.touched_pellet_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==0)'];
        trial2='trialTypes.optoGroup~=1 & (trialTypes.led_1forward==1 | trialTypes.led_2forward==1 | trialTypes.led_3forward==1 | trialTypes.led_4forward==1 | trialTypes.led_1back==1)';
        trial1_LED=['(trialTypes.optoGroup~=1 & trialTypes.consumed_pellet_1back==0 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touched_pellet_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==1) |' ...
                '(trialTypes.optoGroup~=1 & trialTypes.consumed_pellet_1back==0 & trialTypes.reachedBeforeCue_1forward==1 & trialTypes.touched_pellet_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==1)'];
        trial2_LED='trialTypes.optoGroup~=1 & (trialTypes.led_1forward==0 | trialTypes.led_2forward==0 | trialTypes.led_3forward==0 | trialTypes.led_4forward==0 | trialTypes.led_1back==0)';
%         trial1=['(trialTypes.optoGroup_2back~=1 & trialTypes.optoGroup~=1 & trialTypes.reachedInTimeWindow_1back==1 & trialTypes.touched_pellet_1back==0 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touched_pellet_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.isLongITI_1forward==1 & trialTypes.led_1forward==0)'];
%         trial2='trialTypes.optoGroup~=1 & (trialTypes.led_1forward==1 | trialTypes.led_2forward==1 | trialTypes.led_3forward==1 | trialTypes.led_4forward==1 | trialTypes.led_1back==1)';
%         trial1_LED=['(trialTypes.optoGroup_2back~=1 & trialTypes.optoGroup~=1 & trialTypes.reachedInTimeWindow_1back==1 & trialTypes.touched_pellet_1back==0 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touched_pellet_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.isLongITI_1forward==1 & trialTypes.led_1forward==1)'];
%         trial2_LED='trialTypes.optoGroup~=1 & (trialTypes.led_1forward==0 | trialTypes.led_2forward==0 | trialTypes.led_3forward==0 | trialTypes.led_4forward==0 | trialTypes.led_1back==0)';
end


