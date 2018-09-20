function regressReactionTime(alltbt,out,metadata)

% Only take trials where mouse was not chewing at the trial start
% And mouse did not reach preemptively

% Terms
% reaction time 1 back
% reaction time 2 back
% reaction time 3 back
% reaction time 4 back
% led this trial
% led 1 back
% led 2 back
% led 3 back
% led 4 back
% consumed pellet 1 back
% touched pellet 1 back
% cued reach 1 back
% cued reach 2 back
% cued reach 3 back
% cued reach 4 back

% Response
% reaction time on current trial

% Get reaction times for all trials where mouse reached after cue onset
[reactionTimes,alltbt]=plotOnlyFirstReach(alltbt,1,'reachStarts_noPawOnWheel','cueZone_onVoff',out,'led',0);

% Only try to fit these trials
doUse=out.paw_during_wheel==0 & out.chewing_at_trial_start==0;

X=[reactionTimes(1:end-4); reactionTimes(2:end-3); reactionTimes(3:end-2); reactionTimes(4:end-1); ...
   out.led(5:end)'; out.led_1back(5:end)'; out.led_2back(5:end)'; out.led_3back(5:end)'; out.led_4back(5:end)'; ...
   out.consumed_pellet_1back(5:end)'; out.touched_pellet_1back(5:end)'; ...
   out.cued_reach_1back(5:end)'; out.cued_reach_2back(5:end)'; out.cued_reach_3back(5:end)'; out.cued_reach_4back(5:end)']';
response=reactionTimes(5:end)';
response(~doUse)=nan;
X(~doUse,:)=nan(sum(~doUse),size(X,2));
Y=ordinal(response,{'1','2','3','4','5','6','7','8','9','10'},...
                    [],[0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2]);
[B,dev,stats]=mnrfit(X,Y,'Model','ordinal');

disp('pause');
