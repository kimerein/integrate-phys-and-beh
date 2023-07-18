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
trialTypes.isLongITI_2forward=[trialTypes.led(1+2:end); zeros(2,1)];
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
trialTypes.led_5forward=[trialTypes.led(1+5:end); zeros(5,1)];
trialTypes.led_6forward=[trialTypes.led(1+6:end); zeros(6,1)];
trialTypes.led_7forward=[trialTypes.led(1+7:end); zeros(7,1)];
trialTypes.led_8forward=[trialTypes.led(1+8:end); zeros(8,1)];
trialTypes.led_9forward=[trialTypes.led(1+9:end); zeros(9,1)];
trialTypes.led_10forward=[trialTypes.led(1+10:end); zeros(10,1)];
trialTypes.led_11forward=[trialTypes.led(1+11:end); zeros(11,1)];
trialTypes.led_12forward=[trialTypes.led(1+12:end); zeros(12,1)];

linkerForNoLED_accumulate=['((trialTypes.led_1forward==1 & trialTypes.led_2forward==1) | (trialTypes.led_2forward==1 & trialTypes.led_3forward==1) | (trialTypes.led_3forward==1 & trialTypes.led_4forward==1) | (trialTypes.led_4forward==1 & trialTypes.led_5forward==1) | (trialTypes.led_5forward==1 & trialTypes.led_6forward==1)' ...
                         ' | (trialTypes.led_6forward==1 & trialTypes.led_7forward==1) | (trialTypes.led_7forward==1 & trialTypes.led_8forward==1) | (trialTypes.led_8forward==1 & trialTypes.led_9forward==1) | (trialTypes.led_9forward==1 & trialTypes.led_10forward==1) | (trialTypes.led_10forward==1 & trialTypes.led_11forward==1) | (trialTypes.led_11forward==1 & trialTypes.led_12forward==1))'];
linkerForLED_accumulate=['((trialTypes.led_1forward==0 & trialTypes.led_2forward==0) | (trialTypes.led_2forward==0 & trialTypes.led_3forward==0) | (trialTypes.led_3forward==0 & trialTypes.led_4forward==0) | (trialTypes.led_4forward==0 & trialTypes.led_5forward==0) | (trialTypes.led_5forward==0 & trialTypes.led_6forward==0)' ...
                       ' | (trialTypes.led_6forward==0 & trialTypes.led_7forward==0) | (trialTypes.led_7forward==0 & trialTypes.led_8forward==0) | (trialTypes.led_8forward==0 & trialTypes.led_9forward==0) | (trialTypes.led_9forward==0 & trialTypes.led_10forward==0) | (trialTypes.led_10forward==0 & trialTypes.led_11forward==0) | (trialTypes.led_11forward==0 & trialTypes.led_12forward==0))'];
% FOR VARIED TIMING
linkerForVariedTimingForward=' & trialTypes.optoGroup_1forward==2';
linkerForVariedTimingSame=' & trialTypes.optoGroup==2';
linkerForVariedTimingUncued='';
% ELSE
% linkerForVariedTimingForward='';
% linkerForVariedTimingSame='';
% linkerForVariedTimingUncued='';

switch whichEventType
    case 'cued success'
        trial1=['(trialTypes.optoGroup~=1 & trialTypes.consumed_pellet_1back==1 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==0) | ' ...
                '(trialTypes.optoGroup~=1 & trialTypes.isLongITI_1forward==1 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==0)']; 
        trial2='trialTypes.optoGroup~=1 & (trialTypes.led_1forward==1 | trialTypes.led_2forward==1 | trialTypes.led_3forward==1 | trialTypes.led_4forward==1 | trialTypes.led_1back==1)';
        trial1_LED=['(trialTypes.optoGroup~=1 & trialTypes.consumed_pellet_1back==1 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==1' linkerForVariedTimingForward ') | ' ...
                    '(trialTypes.optoGroup~=1 & trialTypes.isLongITI_1forward==1 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==1' linkerForVariedTimingForward ')']; 
        trial2_LED='trialTypes.optoGroup~=1 & (trialTypes.led_1forward==0 | trialTypes.led_2forward==0 | trialTypes.led_3forward==0 | trialTypes.led_4forward==0 | trialTypes.led_1back==0)';
    case 'cued success accumulate'
        % need to use reachedInTimeWindow less than opto_duration and
        % reachBatch_success_reachStarts, reachedInTimeWindow=[0 0.4] works
        trial1=['trialTypes.optoGroup~=1 & (trialTypes.consumed_pellet_1back==1 | trialTypes.isLongITI_Xforward==1) SPLIT trialTypes.optoGroup~=1 & trialTypes.reachedInTimeWindow==1 & trialTypes.led==0'];
        trial2=['trialTypes.optoGroup~=1 & ' linkerForNoLED_accumulate];
        trial1_LED=['trialTypes.optoGroup~=1 & (trialTypes.consumed_pellet_1back==1 | trialTypes.isLongITI_Xforward==1) SPLIT trialTypes.optoGroup~=1 & trialTypes.reachedInTimeWindow==1 & trialTypes.led==1' linkerForVariedTimingSame];
        trial2_LED=['trialTypes.optoGroup~=1 & ' linkerForLED_accumulate];
    case 'backwards cued success'
        % need to use reachedInTimeWindow=[0 opto_duration] and reachBatch_success_reachStarts
        trial2=['(trialTypes.optoGroup~=1 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==0) & (trialTypes.led_2forward==1 | trialTypes.led_3forward==1 | trialTypes.led_4forward==1 | trialTypes.led_5forward==1)'];
        trial1=['trialTypes.optoGroup~=1'];
        trial2_LED=['trialTypes.optoGroup~=1 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==1 & (trialTypes.led_2forward==0 | trialTypes.led_3forward==0 | trialTypes.led_4forward==0 | trialTypes.led_5forward==0)' linkerForVariedTimingForward];
        trial1_LED=['trialTypes.optoGroup~=1'];
    case 'cued failure'
        % did not touch pellet despite reaching in cued window
        trial1=['(trialTypes.optoGroup~=1 & trialTypes.cued_reach_1back==1 & trialTypes.isLongITI_1back==1 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touch_in_cued_window_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==0) | ' ...
                '(trialTypes.optoGroup~=1 & trialTypes.isLongITI_1forward==1                               & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touch_in_cued_window_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==0)']; 
        trial2='trialTypes.optoGroup~=1 & (trialTypes.led_1forward==1 | trialTypes.led_2forward==1 | trialTypes.led_3forward==1 | trialTypes.led_4forward==1 | trialTypes.led_1back==1)';
        trial1_LED=['(trialTypes.optoGroup~=1 & trialTypes.cued_reach_1back==1 & trialTypes.isLongITI_1back==1 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touch_in_cued_window_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==1' linkerForVariedTimingForward ') | ' ...
                    '(trialTypes.optoGroup~=1 & trialTypes.isLongITI_1forward==1                               & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touch_in_cued_window_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==1' linkerForVariedTimingForward ')']; 
        trial2_LED='trialTypes.optoGroup~=1 & (trialTypes.led_1forward==0 | trialTypes.led_2forward==0 | trialTypes.led_3forward==0 | trialTypes.led_4forward==0 | trialTypes.led_1back==0)';
%         trial1=['(trialTypes.optoGroup~=1 & trialTypes.touch_in_cued_window_1back==0 & trialTypes.noReach_1back==0 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touch_in_cued_window_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==0) | ' ...
%                 '(trialTypes.optoGroup~=1 & trialTypes.isLongITI_1forward==1                                       & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touch_in_cued_window_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==0)']; 
%         trial2='trialTypes.optoGroup~=1 & (trialTypes.led_1forward==1 | trialTypes.led_2forward==1 | trialTypes.led_3forward==1 | trialTypes.led_4forward==1 | trialTypes.led_1back==1)';
%         trial1_LED=['(trialTypes.optoGroup~=1 & trialTypes.touch_in_cued_window_1back==0 & trialTypes.noReach_1back==0 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touch_in_cued_window_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==1) | ' ...
%                     '(trialTypes.optoGroup~=1 & trialTypes.isLongITI_1forward==1                                       & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touch_in_cued_window_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==1)']; 
%         trial2_LED='trialTypes.optoGroup~=1 & (trialTypes.led_1forward==0 | trialTypes.led_2forward==0 | trialTypes.led_3forward==0 | trialTypes.led_4forward==0 | trialTypes.led_1back==0)';
    case 'cued failure accumulate'
        trial1=['trialTypes.optoGroup~=1 & ((trialTypes.cued_reach_1back==1 & trialTypes.isLongITI_1back==1) | trialTypes.isLongITI_Xforward==1) SPLIT trialTypes.optoGroup~=1 & trialTypes.reachedInTimeWindow==1 & trialTypes.touch_in_cued_window==0 & trialTypes.led==0'];
        trial2=['trialTypes.optoGroup~=1 & ' linkerForNoLED_accumulate];
        trial1_LED=['trialTypes.optoGroup~=1 & ((trialTypes.cued_reach_1back==1 & trialTypes.isLongITI_1back==1) | trialTypes.isLongITI_Xforward==1) SPLIT trialTypes.optoGroup~=1 & trialTypes.reachedInTimeWindow==1 & trialTypes.touch_in_cued_window==0 & trialTypes.led==1' linkerForVariedTimingSame];
        trial2_LED=['trialTypes.optoGroup~=1 & ' linkerForLED_accumulate];
    case 'backwards cued failure' 
        trial2=['(trialTypes.optoGroup~=1 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touch_in_cued_window_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==0) & (trialTypes.led_2forward==1 | trialTypes.led_3forward==1 | trialTypes.led_4forward==1 | trialTypes.led_5forward==1)']; 
        trial1=['trialTypes.optoGroup~=1'];
        trial2_LED=['(trialTypes.optoGroup~=1 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touch_in_cued_window_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==1 & (trialTypes.led_2forward==0 | trialTypes.led_3forward==0 | trialTypes.led_4forward==0 | trialTypes.led_5forward==0)' linkerForVariedTimingForward ')']; 
        trial1_LED=['trialTypes.optoGroup~=1'];
    case 'delayed success'
        trial1=['(trialTypes.optoGroup~=1 & trialTypes.consumed_pellet_1back==1 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==0) | ' ...
                '(trialTypes.optoGroup~=1 & trialTypes.isLongITI_1forward==1    & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==0)' linkerForVariedTimingUncued]; 
        trial2=['trialTypes.optoGroup~=1 & (trialTypes.led_1forward==1 | trialTypes.led_2forward==1 | trialTypes.led_3forward==1 | trialTypes.led_4forward==1 | trialTypes.led_1back==1)' linkerForVariedTimingUncued];
        trial1_LED=['(trialTypes.optoGroup~=1 & trialTypes.consumed_pellet_1back==1 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==1) | ' ...
                    '(trialTypes.optoGroup~=1 & trialTypes.isLongITI_1forward==1    & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==1)' linkerForVariedTimingForward linkerForVariedTimingUncued]; 
        trial2_LED=['trialTypes.optoGroup~=1 & (trialTypes.led_1forward==0 | trialTypes.led_2forward==0 | trialTypes.led_3forward==0 | trialTypes.led_4forward==0 | trialTypes.led_1back==0)' linkerForVariedTimingUncued];
        % could put chewing_at_trial_start==0 for trial2 and trial2_LED but
        % it does NOT make a difference and reduces trials
    case 'backwards delayed success'
        trial2=['(trialTypes.optoGroup~=1 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==0) & (trialTypes.led_2forward==1 | trialTypes.led_3forward==1 | trialTypes.led_4forward==1 | trialTypes.led_5forward==1)' linkerForVariedTimingUncued]; 
        trial1=['trialTypes.optoGroup~=1' linkerForVariedTimingUncued];
        trial2_LED=['(trialTypes.optoGroup~=1 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==1) & (trialTypes.led_2forward==0 | trialTypes.led_3forward==0 | trialTypes.led_4forward==0 | trialTypes.led_5forward==0)' linkerForVariedTimingForward linkerForVariedTimingUncued]; 
        trial1_LED=['trialTypes.optoGroup~=1' linkerForVariedTimingUncued];
    case 'delayed success accumulate'
        trial1=['trialTypes.optoGroup~=1 & (trialTypes.consumed_pellet_1back==1 | trialTypes.isLongITI_Xforward==1) SPLIT trialTypes.optoGroup~=1 & trialTypes.reachedInTimeWindow==1 & trialTypes.led==0'];
        trial2=['trialTypes.optoGroup~=1 & ' linkerForNoLED_accumulate];
        trial1_LED=['trialTypes.optoGroup~=1 & (trialTypes.consumed_pellet_1back==1 | trialTypes.isLongITI_Xforward==1) SPLIT trialTypes.optoGroup~=1 & trialTypes.reachedInTimeWindow==1 & trialTypes.led==1' linkerForVariedTimingSame];
        trial2_LED=['trialTypes.optoGroup~=1 & ' linkerForLED_accumulate];
    case 'uncued failure'
        % did not touch pellet despite reaching in delayed window after cue
        % OR reached before cue and thus failed to get pellet
        trial1=['(trialTypes.optoGroup~=1 & trialTypes.consumed_pellet_1back==0 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touched_pellet_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==0) |' ...
                '(trialTypes.optoGroup~=1 & trialTypes.consumed_pellet_1back==0 & trialTypes.reachedBeforeCue_1forward==1 & trialTypes.touched_pellet_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==0)' linkerForVariedTimingUncued];
        trial2=['trialTypes.optoGroup~=1 & (trialTypes.led_1forward==1 | trialTypes.led_2forward==1 | trialTypes.led_3forward==1 | trialTypes.led_4forward==1 | trialTypes.led_1back==1)' linkerForVariedTimingUncued];
        trial1_LED=['(trialTypes.optoGroup~=1 & trialTypes.consumed_pellet_1back==0 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touched_pellet_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==1) |' ...
                '(trialTypes.optoGroup~=1 & trialTypes.consumed_pellet_1back==0 & trialTypes.reachedBeforeCue_1forward==1 & trialTypes.touched_pellet_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==1)' linkerForVariedTimingForward linkerForVariedTimingUncued];
        trial2_LED=['trialTypes.optoGroup~=1 & (trialTypes.led_1forward==0 | trialTypes.led_2forward==0 | trialTypes.led_3forward==0 | trialTypes.led_4forward==0 | trialTypes.led_1back==0)' linkerForVariedTimingUncued];
%         trial1=['(trialTypes.optoGroup_2back~=1 & trialTypes.optoGroup~=1 & trialTypes.reachedInTimeWindow_1back==1 & trialTypes.touched_pellet_1back==0 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touched_pellet_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.isLongITI_1forward==1 & trialTypes.led_1forward==0)'];
%         trial2='trialTypes.optoGroup~=1 & (trialTypes.led_1forward==1 | trialTypes.led_2forward==1 | trialTypes.led_3forward==1 | trialTypes.led_4forward==1 | trialTypes.led_1back==1)';
%         trial1_LED=['(trialTypes.optoGroup_2back~=1 & trialTypes.optoGroup~=1 & trialTypes.reachedInTimeWindow_1back==1 & trialTypes.touched_pellet_1back==0 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touched_pellet_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.isLongITI_1forward==1 & trialTypes.led_1forward==1)'];
%         trial2_LED='trialTypes.optoGroup~=1 & (trialTypes.led_1forward==0 | trialTypes.led_2forward==0 | trialTypes.led_3forward==0 | trialTypes.led_4forward==0 | trialTypes.led_1back==0)';
    case 'uncued failure accumulate'
        % note that reachratesettings.acrossSess_window3=[-2 0] picks up
        % more of the uncued decrease
        trial1=['trialTypes.optoGroup~=1 & trialTypes.consumed_pellet_1back==0 SPLIT trialTypes.optoGroup~=1 & (trialTypes.reachedInTimeWindow==1 | trialTypes.reachedBeforeCue==1) & trialTypes.touched_pellet==0 & trialTypes.led==0'];
        trial2=['trialTypes.optoGroup~=1 & ' linkerForNoLED_accumulate];
        trial1_LED=['trialTypes.optoGroup~=1 & trialTypes.consumed_pellet_1back==0 SPLIT trialTypes.optoGroup~=1 & (trialTypes.reachedInTimeWindow==1 | trialTypes.reachedBeforeCue==1) & trialTypes.touched_pellet==0' linkerForVariedTimingSame];
        trial2_LED=['trialTypes.optoGroup~=1 & ' linkerForLED_accumulate];
    case 'backwards uncued failure'
        trial2=['(trialTypes.optoGroup~=1 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touched_pellet_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==0) | ' ...
                '(trialTypes.optoGroup~=1 & trialTypes.reachedBeforeCue_1forward==1 & trialTypes.touched_pellet_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==0)' linkerForVariedTimingUncued ' & (trialTypes.led_2forward==1 | trialTypes.led_3forward==1 | trialTypes.led_4forward==1 | trialTypes.led_5forward==1)'];
        trial1=['trialTypes.optoGroup~=1' linkerForVariedTimingUncued];
        trial2_LED=['(trialTypes.optoGroup~=1 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touched_pellet_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==1) | ' ...
                    '(trialTypes.optoGroup~=1 & trialTypes.reachedBeforeCue_1forward==1 & trialTypes.touched_pellet_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==1)' linkerForVariedTimingUncued ' & (trialTypes.led_2forward==0 | trialTypes.led_3forward==0 | trialTypes.led_4forward==0 | trialTypes.led_5forward==0)'];
        trial1_LED=['trialTypes.optoGroup~=1' linkerForVariedTimingUncued];
    case 'no no opto'
        trial1='trialTypes.optoGroup~=1 & trialTypes.led==0';
        trial2='trialTypes.optoGroup~=1 & trialTypes.led==0';
        trial1_LED='trialTypes.optoGroup~=1 & trialTypes.led==0';
        trial2_LED=['trialTypes.optoGroup~=1 & trialTypes.led==1' linkerForVariedTimingSame];
end


