function settings=RTanalysis_settings()

settings.lowThresh=0.05; % threshold for detecting events in trial-by-trial behavior data structure
settings.durationOfWheelTurn=1.029; % duration of pellet presenter wheel turn, in seconds
settings.wheelStopsThisManySecsBeforeCue=0.1; % time window between pellet presenter wheel stopping and cue onset, in seconds
settings.cueDuration=0.25; % duration of cue, in seconds
settings.maxTrialDuration=35; % maximum duration of a trial, in seconds
settings.timeSlop=0.1; % max possible error in timing due to frame rate in movie and/or alignment approach, in seconds
settings.reachAfterCueWindow_start=0.15; % define start of time window in which to say mouse "reached after cue", seconds from cue onset
settings.reachAfterCueWindow_end=3.5; % define end of time window in which to say mouse "reached after cue", seconds from cue onset 
settings.preCueWindow_start=0; % define start of time window from trial onset, in seconds
settings.preCueWindow_end=1.5; % define end of time window from trial onset, in seconds

settings.returnPreCueRTs=0; % also return "reaction times" pretending that time moves backwards from cue
settings.divideByBaseRate=0; % note that reaction times pick first reach after cue -- higher baseline reaching rate will mean more RTs near 0 ... do you want to normalize reaction time distribution by base rate? if yes, then 1
settings.subtractBaselineReaching=0; % if is 1, will subtract off the baseline reaching rate, effectively zeroing out time bins where distribution falls at or below baseline
settings.throwOutAllAfter1stBaseline=0; % if 1, will throw out all reaction times after reaction time distribution first returns to baseline reach rate
settings.noRTlessThan=0; % no reaction times can be less than this (in seconds); if 0, this is the classic definition of reaction time

% settings.preCueWindow_start=3; % define start of time window from trial onset, in seconds
% settings.preCueWindow_end=4.5; % define end of time window from trial onset, in seconds
% settings.preCueWindow_start={[0 1.5],[3 14]}; % define start of time window from trial onset, in seconds
% settings.preCueWindow_end=1.5; % doesn't matter
settings.excludePawOnWheelDuringCue=0; % 1 if want to exclude trials where mouse reached during cue, else 0
settings.longRT_ifNoReach=1; % if mouse does not reach in this trial, fill in reaction time as longer than trial length IF this is set to 1, else throw out trial

% % Settings on 4/7/2019 before changing
% settings.lowThresh=0.05; % threshold for detecting events in trial-by-trial behavior data structure
% settings.durationOfWheelTurn=1.029; % duration of pellet presenter wheel turn, in seconds
% settings.wheelStopsThisManySecsBeforeCue=0.1; % time window between pellet presenter wheel stopping and cue onset, in seconds
% settings.excludePawOnWheelDuringCue=1; % 1 if want to exclude trials where mouse reached during cue, else 0
% settings.cueDuration=0.25; % duration of cue, in seconds
% settings.reachAfterCueWindow_start=0; % define start of time window in which to say mouse "reached after cue", seconds from cue onset
% settings.reachAfterCueWindow_end=1.5; % define end of time window in which to say mouse "reached after cue", seconds from cue onset 
% settings.maxTrialDuration=35; % maximum duration of a trial, in seconds
% settings.timeSlop=0.1; % max possible error in timing due to frame rate in movie and/or alignment approach, in seconds
% settings.preCueWindow_start=0; % define start of time window from trial onset, in seconds
% settings.preCueWindow_end=1.5; % define end of time window from trial onset, in seconds
% settings.longRT_ifNoReach=0; % if mouse does not reach in this trial, fill in reaction time as longer than trial length IF this is set to 1, else throw out trial
% settings.returnPreCueRTs=1; % also return "reaction times" pretending that time moves backwards from cue
% settings.divideByBaseRate=0; % note that reaction times pick first reach after cue -- higher baseline reaching rate will mean more RTs near 0 ... do you want to normalize reaction time distribution by base rate? if yes, then 1
% settings.subtractBaselineReaching=1; % if is 1, will subtract off the baseline reaching rate, effectively zeroing out time bins where distribution falls at or below baseline
% settings.throwOutAllAfter1stBaseline=0; % if 1, will throw out all reaction times after reaction time distribution first returns to baseline reach rate
% settings.noRTlessThan=0; % no reaction times can be less than this (in seconds); if 0, this is the classic definition of reaction time
% % settings.preCueWindow_start=3; % define start of time window from trial onset, in seconds
% % settings.preCueWindow_end=4.5; % define end of time window from trial onset, in seconds
% % settings.preCueWindow_start={[0 1.5],[3 14]}; % define start of time window from trial onset, in seconds
% % settings.preCueWindow_end=1.5; % doesn't matter