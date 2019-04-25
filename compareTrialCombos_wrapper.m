function out=compareTrialCombos_wrapper(alltbt,trialTypes,metadata)

% trial type 1
% templateSequence1{1}=trialTypes.chewing_at_trial_start==0 & trialTypes.touch_in_cued_window==1 & trialTypes.led==0 & trialTypes.nWheelTurns<3;
% templateSequence1{2}=trialTypes.chewing_at_trial_start==0 & trialTypes.paw_during_wheel==0;

% templateSequence1{1}=trialTypes.touch_in_cued_window==1 & trialTypes.cued_reach==1 & trialTypes.led==0 & trialTypes.nWheelTurns>2;
% templateSequence1{2}=trialTypes.chewing_at_trial_start==0 & trialTypes.paw_during_wheel==0;

% templateSequence1{1}=trialTypes.after_cue_drop==1 & trialTypes.led==0; 
% templateSequence1{2}=trialTypes.after_cue_drop==1 & trialTypes.led==0;
% templateSequence1{3}=trialTypes.chewing_at_trial_start==0 & trialTypes.paw_during_wheel==0 & trialTypes.dprimes>1;

% templateSequence1{1}=trialTypes.touched_pellet==1 & trialTypes.cued_reach==1 & trialTypes.chewing_at_trial_start==0 & trialTypes.led==0 & trialTypes.dprimes>1 & trialTypes.dprimes<4;
% templateSequence1{2}=trialTypes.chewing_at_trial_start==0 & trialTypes.paw_during_wheel==0 & trialTypes.dprimes>1 & trialTypes.dprimes<4;
% nNext1=1;


% trial type 2
% templateSequence2{1}=trialTypes.chewing_at_trial_start==0 & trialTypes.touch_in_cued_window==1 & trialTypes.led==1 & trialTypes.nWheelTurns<3;
% templateSequence2{2}=trialTypes.chewing_at_trial_start==0 & trialTypes.paw_during_wheel==0;

% templateSequence2{1}=trialTypes.touch_in_cued_window==1 & trialTypes.cued_reach==1 & trialTypes.led==1 & trialTypes.nWheelTurns>2;
% templateSequence2{2}=trialTypes.chewing_at_trial_start==0 & trialTypes.paw_during_wheel==0;

% templateSequence2{1}=trialTypes.after_cue_drop==1 & (trialTypes.led==1 | trialTypes.led_1forward==1); 
% templateSequence2{2}=trialTypes.after_cue_drop==1 & (trialTypes.led==1 | trialTypes.led_1back==1);
% templateSequence2{3}=trialTypes.after_cue_drop==0 & trialTypes.paw_during_wheel==0 & trialTypes.dprimes>1;

% templateSequence2{1}=trialTypes.touched_pellet==1 & trialTypes.cued_reach==1 & trialTypes.chewing_at_trial_start==0 & trialTypes.led==1 & trialTypes.dprimes>1 & trialTypes.dprimes<4;
% templateSequence2{2}=trialTypes.chewing_at_trial_start==0 & trialTypes.paw_during_wheel==0 & trialTypes.dprimes>1 & trialTypes.dprimes<4;
% nNext2=1;



% template sequence 1 
% templateSequence1{1}=trialTypes.consumed_pellet==1 & trialTypes.led==0;
% templateSequence1{2}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
% templateSequence1{2}=trialTypes.chewing_at_trial_start==0 & trialTypes.paw_during_wheel==0;
templateSequence1{1}=trialTypes.touched_pellet==1 & trialTypes.cued_reach==1 & trialTypes.led==0;
templateSequence1{2}=trialTypes.chewing_at_trial_start==0 & trialTypes.paw_during_wheel==0;
nNext1=1;

% template sequence 2 
% templateSequence2{1}=trialTypes.consumed_pellet==1 & trialTypes.led==1;
% templateSequence2{2}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
% templateSequence2{2}=trialTypes.chewing_at_trial_start==0 & trialTypes.paw_during_wheel==0;
templateSequence2{1}=trialTypes.touched_pellet==1 & trialTypes.cued_reach==1 & trialTypes.led==1;
templateSequence2{2}=trialTypes.chewing_at_trial_start==0 & trialTypes.paw_during_wheel==0;
nNext2=1;

% Format 3 (10 args); tbt, trialTypes, metadata, templateSequence1,
%   nNext1, templateSequence2, nNext2, useAsCue, whichReach, doPlot
[~,~,out]=compareTrialCombos(alltbt,trialTypes,metadata,templateSequence1,nNext1,templateSequence2,nNext2,'cueZone_onVoff','all_reachBatch',1);