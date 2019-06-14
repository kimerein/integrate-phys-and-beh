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
% templateSequence1{1}=trialTypes.touched_pellet==1 & trialTypes.cued_reach==1 & trialTypes.led==0;
% templateSequence1{2}=trialTypes.chewing_at_trial_start==0 & trialTypes.paw_during_wheel==0;
sumAcross=25;
sumsOfTouched=nan(1,length(trialTypes.touch_in_cued_window)-(sumAcross-1));
for i=1:length(sumsOfTouched)
    sumsOfTouched(i)=sum(trialTypes.touch_in_cued_window(i:i+(sumAcross-1)));
%     sumsOfTouched(i)=sum(trialTypes.touch_in_cued_window(i:i+(sumAcross-1))==0);
end
sumsOfTouched=[sumsOfTouched zeros(1,sumAcross-1)]';

sumAcrossLED=25;
sumsOfLED=nan(1,length(trialTypes.led)-(sumAcross-1));
for i=1:length(sumsOfLED)
    sumsOfLED(i)=sum(trialTypes.led(i:i+(sumAcrossLED-1)));
end
sumsOfLED=[sumsOfLED zeros(1,sumAcrossLED-1)]';

% templateSequence1{1}=trialTypes.paw_during_wheel==0;
% templateSequence1{2}=(sumsOfTouched>=1) & (sumsOfLED==0) & trialTypes.paw_during_wheel==0 & any(alltbt.isChewing(:,200:486),2)==0;
% % templateSequence1{2}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
% templateSequence1{3}=trialTypes.paw_during_wheel==0;
% templateSequence1{4}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
% templateSequence1{5}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
% templateSequence1{6}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
% templateSequence1{7}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
% templateSequence1{8}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
% templateSequence1{9}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
% templateSequence1{10}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
% templateSequence1{11}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
% templateSequence1{12}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;

templateSequence1{1}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
% templateSequence1{2}=(sumsOfTouched>=4) & (sumsOfLED<1);
% templateSequence1{2}=(sumsOfTouched>=4) & (sumsOfLED>=3);
templateSequence1{2}=(sumsOfTouched>=20) & (sumsOfLED>=12);
% templateSequence1{2}=(sumsOfTouched>=20) & (sumsOfLED<16);
% templateSequence1{2}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
templateSequence1{3}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
templateSequence1{4}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
templateSequence1{5}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
templateSequence1{6}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
templateSequence1{7}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
templateSequence1{8}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
templateSequence1{9}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
templateSequence1{10}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
templateSequence1{11}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
templateSequence1{12}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;

nNext1=1;

% template sequence 2 
% templateSequence2{1}=trialTypes.consumed_pellet==1 & trialTypes.led==1;
% templateSequence2{2}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
% templateSequence2{2}=trialTypes.chewing_at_trial_start==0 & trialTypes.paw_during_wheel==0;
% templateSequence2{1}=trialTypes.touched_pellet==1 & trialTypes.cued_reach==1 & trialTypes.led==1;
% templateSequence2{2}=trialTypes.chewing_at_trial_start==0 & trialTypes.paw_during_wheel==0;

% templateSequence2{1}=trialTypes.paw_during_wheel==0;
% templateSequence2{2}=(sumsOfTouched>=1) & (sumsOfLED>=1) & trialTypes.paw_during_wheel==0 & any(alltbt.isChewing(:,200:486),2)==0;
% % templateSequence2{2}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
% templateSequence2{3}=trialTypes.paw_during_wheel==0;
% templateSequence2{4}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
% templateSequence2{5}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
% templateSequence2{6}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
% templateSequence2{7}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
% templateSequence2{8}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
% templateSequence2{9}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
% templateSequence2{10}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
% templateSequence2{11}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
% templateSequence2{12}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;

templateSequence2{1}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
% templateSequence2{2}=(sumsOfTouched<=3);
templateSequence2{2}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
templateSequence2{3}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
templateSequence2{4}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
templateSequence2{5}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
templateSequence2{6}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
templateSequence2{7}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
templateSequence2{8}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
templateSequence2{9}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
templateSequence2{10}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
templateSequence2{11}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
templateSequence2{12}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;

nNext2=1;

% Format 3 (10 args); tbt, trialTypes, metadata, templateSequence1,
%   nNext1, templateSequence2, nNext2, useAsCue, whichReach, doPlot
[~,~,out]=compareTrialCombos(alltbt,trialTypes,metadata,templateSequence1,nNext1,templateSequence2,nNext2,'cueZone_onVoff','all_reachBatch',1);