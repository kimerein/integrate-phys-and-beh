function out=getSweepsFromBeh(tbt,cueName)

lowThresh=0.05;
durationOfWheelTurn=1; % in seconds
durationOfWheelTurn_inds=floor(durationOfWheelTurn/mode(diff(nanmean(tbt.times,1))));
includePawOnWheelDuringCue=1;
cueDuration=0.2; % in seconds
cueDuration_inds=floor(cueDuration/mode(diff(nanmean(tbt.times,1))));
reachAfterCueWindow=1; % in seconds
reachAfterCueWindow_inds=floor(reachAfterCueWindow/mode(diff(nanmean(tbt.times,1))));

% Get trial durations from tbt
figure();
plot(nanmean(tbt.(cueName),1));
xlabel('indices');
ylabel('av');
title('average of movie cue across trials');
ind=input('Enter index showing minimum trial length -- note that any opto within a trial should occur before this index. ');
if ~isnumeric(ind)
    error('Please enter a number.');
end

% Get trials with opto
led=any(tbt.optoOn(:,1:ind)>lowThresh,2);
led=single(led);

% Save to structure for output
out.led=led;

% Opto status from previous trial
out.led_previousTrial=[nan; led(1:end-1)];
out.led_2back=[nan; nan; led(1:end-2)];
out.led_3back=[nan; nan; nan; led(1:end-3)];
out.led_4back=[nan; nan; nan; nan; led(1:end-4)];

% Get trials where reach occurred before cue
% cueInd=find(nanmean(tbt.cueZone_onVoff,1)>0,1,'first');
pawOutDuringWheel=zeros(size(tbt.cueZone_onVoff,1),1);
out.rewarded=nan(size(tbt.cueZone_onVoff,1),1);
out.reachedAfterCue=nan(size(tbt.cueZone_onVoff,1),1);
for i=1:size(tbt.cueZone_onVoff,1)
    presentInd=find(tbt.pelletPresented(i,:)>lowThresh,1,'first');
    temp=tbt.cueZone_onVoff;
    cueInd=find(temp(i,:)>lowThresh,1,'first');
    if includePawOnWheelDuringCue==1
        cueInd=cueInd+cueDuration_inds;
    end
    if cueInd-durationOfWheelTurn_inds<1
        if any(tbt.reach_ongoing(i,1:cueInd)>lowThresh) || any(tbt.isHold(i,1:cueInd)>lowThresh) || any(tbt.reachStarts(i,1:cueInd)>lowThresh)
            pawOutDuringWheel(i)=1;
        else
            pawOutDuringWheel(i)=0;
        end
%     elseif any(tbt.reach_ongoing(i,cueInd-durationOfWheelTurn_inds:cueInd)>lowThresh) || any(tbt.reachStarts(i,cueInd-durationOfWheelTurn_inds:cueInd)>lowThresh)
    elseif any(tbt.reach_ongoing(i,cueInd-durationOfWheelTurn_inds:cueInd)>lowThresh) || any(tbt.isHold(i,cueInd-durationOfWheelTurn_inds:cueInd)>lowThresh) || any(tbt.reachStarts(i,cueInd-durationOfWheelTurn_inds:cueInd)>lowThresh)
        pawOutDuringWheel(i)=1;
    else
        pawOutDuringWheel(i)=0;
    end
    
    % If mouse got a pellet successfully in this trial
    if any(tbt.success_reachStarts(i,cueInd:ind)>lowThresh) || any(tbt.success_reachStarts_pawOnWheel(i,cueInd:ind)>lowThresh)
        out.rewarded(i)=1;
    else
        out.rewarded(i)=0;
    end
    
    % If mouse reached within a few seconds of cue
    if any(tbt.reachStarts(i,cueInd:cueInd+reachAfterCueWindow_inds)>lowThresh)
        out.reachedAfterCue(i)=1;
    else
        out.reachedAfterCue(i)=0;
    end
end
out.pawOutDuringWheel=pawOutDuringWheel;

% Rewarded status and behavior from previous trial
out.rewarded_previousTrial=[nan; out.rewarded(1:end-1)];
out.reachedAfterCue_previousTrial=[nan; out.reachedAfterCue(1:end-1)];
out.rewarded_2back=[nan; nan; out.rewarded(1:end-2)];
out.rewarded_3back=[nan; nan; nan; out.rewarded(1:end-3)];
out.rewarded_4back=[nan; nan; nan; nan; out.rewarded(1:end-4)];
out.reachedAfterCue_2back=[nan; nan; out.reachedAfterCue(1:end-2)];
out.reachedAfterCue_3back=[nan; nan; nan; out.reachedAfterCue(1:end-3)];
out.reachedAfterCue_4back=[nan; nan; nan; nan; out.reachedAfterCue(1:end-4)];

% Classify trial types
% out.trialTypes=classifyTrialTypes(tbt,led);

end
