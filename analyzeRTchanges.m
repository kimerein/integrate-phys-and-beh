function [alltbt,out,rt_noLED,useTheseTrials_noLED,rt_LED,useTheseTrials_LED]=analyzeRTchanges(filedir)

alltbt=combineExptPieces(filedir,'cueZone_onVoff',0.25);

out=getSweepsFromBeh(alltbt,'cueZone_onVoff');

% out.reachedLastTimeNoLED_noPawOut=out.led_previousTrial==0  & out.reachedAfterCue_previousTrial==1 & out.touchedPellet_previousTrial==1 & out.rewarded_previousTrial==1 & out.pawOutDuringWheel_previousTrial==0;
% out.reachedLastTimeNoLED_noPawOut=out.led_previousTrial==0  & out.reachedAfterCue_previousTrial==1 & out.touchedPellet_previousTrial==1 & out.pawOutDuringWheel_previousTrial==0;
out.reachedLastTimeNoLED_noPawOut=out.led==0;
[rt_noLED,useTheseTrials_noLED]=plotOnlyFirstReach(alltbt,1,'reachStarts','cueZone_onVoff',out,'reachedLastTimeNoLED_noPawOut',1);
useTheseTrials_noLED=useTheseTrials_noLED==1 & out.reachedLastTimeNoLED_noPawOut==1;
f=find(useTheseTrials_noLED==1 & out.pawOutDuringWheel==0);
% f=find(out.led_previousTrial==0 & out.reachedAfterCue_previousTrial==1 & out.touchedPellet_previousTrial==1);
changeInRT=rt_noLED(f-1)-rt_noLED(f);

% out.reachedLastTimeAndLED_noPawOut=out.led_previousTrial==1  & out.reachedAfterCue_previousTrial==1 & out.touchedPellet_previousTrial==1 & out.rewarded_previousTrial==1 & out.pawOutDuringWheel_previousTrial==0;
% out.reachedLastTimeAndLED_noPawOut=out.led_previousTrial==1  & out.reachedAfterCue_previousTrial==1 & out.touchedPellet_previousTrial==1 & out.pawOutDuringWheel_previousTrial==0;
out.reachedLastTimeAndLED_noPawOut=out.led==1;
[rt_LED,useTheseTrials_LED]=plotOnlyFirstReach(alltbt,1,'reachStarts','cueZone_onVoff',out,'reachedLastTimeAndLED_noPawOut',1);
useTheseTrials_LED=useTheseTrials_LED==1 & out.reachedLastTimeAndLED_noPawOut==1;
f=find(useTheseTrials_LED==1 & out.pawOutDuringWheel==0);
% f=find(out.led_previousTrial==1 & out.reachedAfterCue_previousTrial==1 & out.touchedPellet_previousTrial==1);
changeInRT_LED=rt_LED(f-1)-rt_LED(f);

% figure();
% hold on; f=find(useTheseTrials_noLED==1 & out.pawOutDuringWheel==0);
% figure(); scatter(rt_noLED(f-1)+rand(size(rt_noLED(f-1))).*0.05,rt_noLED(f)+rand(size(rt_noLED(f))).*0.05,[],'b','filled');
% f=find(useTheseTrials_LED==1 & out.pawOutDuringWheel==0);
% hold on; scatter(rt_LED(f-1)+rand(size(rt_LED(f-1))).*0.05,rt_LED(f)+rand(size(rt_LED(f))).*0.05,[],'c','filled');
% line([0 10],[0 10]);
% 
% p=ranksum(changeInRT,changeInRT_LED);
% disp('p-value of ranksum');
% disp(p);

out.noReachLastTimeNoLED_noPawOut=out.led_previousTrial==0 & out.reachedAfterCue_previousTrial==0;
plotOnlyFirstReach(alltbt,1,'reachStarts','cueZone_onVoff',out,'noReachLastTimeNoLED_noPawOut',1);

out.noReachLastTimeAndLED_noPawOut=out.led_previousTrial==1 & out.reachedAfterCue_previousTrial==0;
plotOnlyFirstReach(alltbt,1,'reachStarts','cueZone_onVoff',out,'noReachLastTimeAndLED_noPawOut',1);
