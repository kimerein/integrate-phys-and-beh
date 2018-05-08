function [alltbt,out,rt_noLED,useTheseTrials_noLED,rt_LED,useTheseTrials_LED]=analyzeRTchanges(filedir)

alltbt=combineExptPieces(filedir,'cueZone_onVoff',0.25);

out=getSweepsFromBeh(alltbt,'cueZone_onVoff');

% out.reachedLastTimeNoLED_noPawOut=out.led_previousTrial==0 & out.reachedAfterCue_previousTrial==1 & out.pawOutDuringWheel_previousTrial==0  & out.touchedPellet_previousTrial==1;
out.reachedLastTimeNoLED_noPawOut=out.led_previousTrial==0  & out.reachedAfterCue_previousTrial==1 & out.touchedPellet_previousTrial==1 & out.rewarded_previousTrial==0 & out.pawOutDuringWheel_previousTrial==0;
[rt_noLED,useTheseTrials_noLED]=plotOnlyFirstReach(alltbt,1,'reachStarts','cueZone_onVoff',out,'reachedLastTimeNoLED_noPawOut',1);
useTheseTrials_noLED=useTheseTrials_noLED==1 & out.reachedLastTimeNoLED_noPawOut==1;
f=find(useTheseTrials_noLED==1 & out.pawOutDuringWheel==0);
% f=find(out.led_previousTrial==0 & out.reachedAfterCue_previousTrial==1 & out.touchedPellet_previousTrial==1);
changeInRT=rt_noLED(f-1)-rt_noLED(f);

% out.reachedLastTimeAndLED_noPawOut=out.led_previousTrial==1 & out.reachedAfterCue_previousTrial==1 & out.pawOutDuringWheel_previousTrial==0  & out.touchedPellet_previousTrial==1;
out.reachedLastTimeAndLED_noPawOut=out.led_previousTrial==1  & out.reachedAfterCue_previousTrial==1 & out.touchedPellet_previousTrial==1 & out.rewarded_previousTrial==0 & out.pawOutDuringWheel_previousTrial==0;
[rt_LED,useTheseTrials_LED]=plotOnlyFirstReach(alltbt,1,'reachStarts','cueZone_onVoff',out,'reachedLastTimeAndLED_noPawOut',1);
useTheseTrials_LED=useTheseTrials_LED==1 & out.reachedLastTimeAndLED_noPawOut==1;
f=find(useTheseTrials_LED==1 & out.pawOutDuringWheel==0);
% f=find(out.led_previousTrial==1 & out.reachedAfterCue_previousTrial==1 & out.touchedPellet_previousTrial==1);
changeInRT_LED=rt_LED(f-1)-rt_LED(f);


p=ranksum(changeInRT,changeInRT_LED);
disp('p-value of ranksum');
disp(p);

out.noReachLastTimeNoLED_noPawOut=out.led_previousTrial==0 & out.reachedAfterCue_previousTrial==0;
plotOnlyFirstReach(alltbt,1,'reachStarts','cueZone_onVoff',out,'noReachLastTimeNoLED_noPawOut',1);

out.noReachLastTimeAndLED_noPawOut=out.led_previousTrial==1 & out.reachedAfterCue_previousTrial==0;
plotOnlyFirstReach(alltbt,1,'reachStarts','cueZone_onVoff',out,'noReachLastTimeAndLED_noPawOut',1);
