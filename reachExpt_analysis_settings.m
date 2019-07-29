function settings=reachExpt_analysis_settings()

settings.lowThresh=0.05; % threshold for detecting events in trial-by-trial behavior data structure
settings.durationOfWheelTurn=1.029; % duration of pellet presenter wheel turn, in seconds
settings.wheelStopsThisManySecsBeforeCue=0.1; % time window between pellet presenter wheel stopping and cue onset, in seconds
settings.cueDuration=0.25; % duration of cue, in seconds
settings.maxTrialDuration=35; % maximum duration of a trial, in seconds
settings.timeSlop=0.1; % max possible error in timing due to frame rate in movie and/or alignment approach, in seconds
settings.reachAfterCueWindow_start=0.05; % define start of time window in which to say mouse "reached after cue", seconds from cue onset
settings.reachAfterCueWindow_end=5; % define end of time window in which to say mouse "reached after cue", seconds from cue onset 
settings.preCueWindow_start=0; % define start of time window from trial onset, in seconds
settings.preCueWindow_end=1.5; % define end of time window from trial onset, in seconds