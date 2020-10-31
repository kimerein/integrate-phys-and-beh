function settings=trialTypeSettings(RTsettings)

% This function defines types of behavioral trials
% Modify this function to change trial-type classification scheme

if isempty(RTsettings)
    RTsettings=reachExpt_analysis_settings(); % go to this file to change experiment-specific parameters
end


% Threshold for event detection
settings.lowThresh=RTsettings.lowThresh; % threshold for detecting events in trial-by-trial behavior data structure


% Name of cue in data structure
settings.nameOfCue=RTsettings.nameOfCue;


% Get information about trial structure from RTsettings
trialStructure.durationOfWheelTurn=RTsettings.durationOfWheelTurn; % duration of pellet presenter wheel turn, in seconds 
trialStructure.wheelStopsThisManySecsBeforeCue=RTsettings.wheelStopsThisManySecsBeforeCue; % time window between pellet presenter wheel stopping and cue onset, in seconds
trialStructure.cueDuration=RTsettings.cueDuration; % duration of cue, in seconds
trialStructure.maxTrialDuration=RTsettings.maxTrialDuration; % maximum duration of a trial, in seconds
trialStructure.maxDelayUntilOpto=RTsettings.maxDelayUntilOpto; % maximum time from trial onset until opto
trialStructure.timeSlop=RTsettings.timeSlop; % max possible error in timing due to frame rate in movie and/or alignment approach, in seconds


% Define time windows within a trial RELATIVE TO CUE ONSET
% Mouse "reached after cue" if mouse reached within cued_reach_window:
timeWindow(1).name='cued_reach_window';
timeWindow(1).start=RTsettings.reachAfterCueWindow_start; % define start of time window in which to say mouse "reached after cue", seconds from cue onset
timeWindow(1).end=RTsettings.reachAfterCueWindow_end;  % define end of time window in which to say mouse "reached after cue", seconds from cue onset 
% wheel_turning is time window relative to cue onset when wheel is turning
timeWindow(2).name='wheel_turning';
timeWindow(2).start=-trialStructure.durationOfWheelTurn-trialStructure.wheelStopsThisManySecsBeforeCue;
timeWindow(2).end=-trialStructure.wheelStopsThisManySecsBeforeCue;
% after_cue is time window from cue onset to end of trial
timeWindow(3).name='after_cue';
timeWindow(3).start=0-trialStructure.timeSlop; % in seconds from cue onset
timeWindow(3).end=trialStructure.maxTrialDuration; 
% all_times refers to entire trial
timeWindow(4).name='all_times';
timeWindow(4).start=0-trialStructure.durationOfWheelTurn-trialStructure.wheelStopsThisManySecsBeforeCue; % trial begins when pellet presenter wheel begins to turn
timeWindow(4).end=trialStructure.maxTrialDuration; 
% trial_start refers to beginning of trial
timeWindow(5).name='trial_start';
timeWindow(5).start=0-3;
timeWindow(5).end=0-2.3;


% Result of first "reach batch" after cue
% Note that "reach batch" refers to the outcome of the last reach in set of
% contiguous reaches, per definition above
i=1;
bool_test(i).testwhat='reach batch';
bool_test(i).fieldname='success_reachStarts';
bool_test(i).test='first';
bool_test(i).inEventSet='all_reachBatch';
bool_test(i).thresh=settings.lowThresh;
bool_test(i).comparator='>';
bool_test(i).window='after_cue';

i=2;
bool_test(i).testwhat='reach batch';
bool_test(i).fieldname='drop_reachStarts';
bool_test(i).test='first';
bool_test(i).inEventSet='all_reachBatch';
bool_test(i).thresh=settings.lowThresh;
bool_test(i).comparator='>';
bool_test(i).window='after_cue';

i=3;
bool_test(i).testwhat='reach batch';
bool_test(i).fieldname='miss_reachStarts';
bool_test(i).test='first';
bool_test(i).inEventSet='all_reachBatch';
bool_test(i).thresh=settings.lowThresh;
bool_test(i).comparator='>';
bool_test(i).window='after_cue';

i=4;
bool_test(i).testwhat='single reach';
bool_test(i).fieldname='pelletmissingreach_reachStarts';
bool_test(i).test='first';
bool_test(i).inEventSet='all_reachBatch';
bool_test(i).thresh=settings.lowThresh;
bool_test(i).comparator='>';
bool_test(i).window='after_cue';

% Did mouse consume pellet on this trial?
% Build up boolean tests to assess this
i=5;
bool_test(i).testwhat='single reach';
bool_test(i).fieldname='success_reachStarts_pawOnWheel';
bool_test(i).test='any';
bool_test(i).thresh=settings.lowThresh;
bool_test(i).comparator='>';
bool_test(i).window='all_times';

i=6;
bool_test(i).testwhat='single reach';
bool_test(i).fieldname='success_reachStarts';
bool_test(i).test='any';
bool_test(i).thresh=settings.lowThresh;
bool_test(i).comparator='>';
bool_test(i).window='all_times';

% Did mouse touch pellet on this trial?
% Build up boolean tests to assess this
i=7;
bool_test(i).testwhat='single reach';
bool_test(i).fieldname='drop_reachStarts_pawOnWheel';
bool_test(i).test='any';
bool_test(i).thresh=settings.lowThresh;
bool_test(i).comparator='>';
bool_test(i).window='all_times';

i=8;
bool_test(i).testwhat='single reach';
bool_test(i).fieldname='drop_reachStarts';
bool_test(i).test='any';
bool_test(i).thresh=settings.lowThresh;
bool_test(i).comparator='>';
bool_test(i).window='all_times';

% Did mouse perform a cued reach?
% Build up boolean tests to assess this
i=9;
bool_test(i).testwhat='single reach';
bool_test(i).fieldname='all_reachBatch';
bool_test(i).test='any';
bool_test(i).thresh=settings.lowThresh;
bool_test(i).comparator='>';
bool_test(i).window='cued_reach_window';

% Test whether paw was on wheel while wheel turning
i=10;
bool_test(i).testwhat='single reach';
bool_test(i).fieldname='reach_ongoing';
bool_test(i).test='any';
bool_test(i).thresh=settings.lowThresh;
bool_test(i).comparator='>';
bool_test(i).window='wheel_turning';

i=11;
bool_test(i).testwhat='single reach';
bool_test(i).fieldname='isHold';
bool_test(i).test='any';
bool_test(i).thresh=settings.lowThresh;
bool_test(i).comparator='>';
bool_test(i).window='wheel_turning';

i=12;
bool_test(i).testwhat='single reach';
bool_test(i).fieldname='reachStarts';
bool_test(i).test='any';
bool_test(i).thresh=settings.lowThresh;
bool_test(i).comparator='>';
bool_test(i).window='wheel_turning';

% Was mouse still chewing at the beginning of this trial?
i=13;
bool_test(i).testwhat='single reach';
bool_test(i).fieldname='isChewing';
bool_test(i).test='any';
bool_test(i).thresh=settings.lowThresh;
bool_test(i).comparator='>';
bool_test(i).window='trial_start';

% Any reach batch successes or drops for touching pellet check
i=14;
bool_test(i).testwhat='reach batch';
bool_test(i).fieldname='success_reachStarts';
bool_test(i).test='any';
bool_test(i).thresh=settings.lowThresh;
bool_test(i).comparator='>';
bool_test(i).window='all_times';

i=15;
bool_test(i).testwhat='reach batch';
bool_test(i).fieldname='drop_reachStarts';
bool_test(i).test='any';
bool_test(i).thresh=settings.lowThresh;
bool_test(i).comparator='>';
bool_test(i).window='all_times';

% Any success in cued window
i=16;
bool_test(i).testwhat='reach batch';
bool_test(i).fieldname='success_reachStarts';
bool_test(i).test='any';
bool_test(i).thresh=settings.lowThresh;
bool_test(i).comparator='>';
bool_test(i).window='cued_reach_window';

% Any drop in cued window
i=17;
bool_test(i).testwhat='reach batch';
bool_test(i).fieldname='drop_reachStarts';
bool_test(i).test='any';
bool_test(i).thresh=settings.lowThresh;
bool_test(i).comparator='>';
bool_test(i).window='cued_reach_window';

% Put together bool_test results to get trial classifications
% Each index into trialtype(i).outcomes vector of length n refers to outcome of bool_tests 1:n
% If index i is nan, bool_tests(i) can be either 0 or 1 (false or true)
trialtype(1).outcomes=[1 0 0 0 nan nan nan nan nan nan nan nan nan nan nan nan nan];
% trialtype(1).outcomes=[1 0 0 0 nan nan nan nan nan 0 0 0 nan nan nan];
trialtype(1).name='after_cue_success';
trialtype(1).bool='&';

trialtype(2).outcomes=[0 1 0 0 nan nan nan nan nan nan nan nan nan nan nan nan nan];
% trialtype(2).outcomes=[0 1 0 0 nan nan nan nan nan 0 0 0 nan nan nan];
trialtype(2).name='after_cue_drop';
trialtype(2).bool='&';

trialtype(3).outcomes=[0 0 1 0 nan nan nan nan nan nan nan nan nan nan nan nan nan];
% trialtype(3).outcomes=[0 0 1 0 nan nan nan nan nan 0 0 0 nan nan nan];
trialtype(3).name='after_cue_miss';
trialtype(3).bool='&';

trialtype(4).outcomes=[0 0 0 1 nan nan nan nan nan nan nan nan nan nan nan nan nan];
% trialtype(4).outcomes=[0 0 0 1 nan nan nan nan nan 0 0 0 nan nan nan];
trialtype(4).name='after_cue_no_pellet';
trialtype(4).bool='&';

trialtype(5).outcomes=[nan nan nan nan 1 1 nan nan nan nan nan nan nan 1 nan nan nan];
trialtype(5).name='consumed_pellet';
trialtype(5).bool='|';

trialtype(6).outcomes=[nan nan nan nan 1 1 1 1 nan nan nan nan nan 1 1 nan nan];
trialtype(6).name='touched_pellet';
trialtype(6).bool='|';

trialtype(7).outcomes=[nan nan nan nan nan nan nan nan 1 nan nan nan nan nan nan nan nan];
trialtype(7).name='cued_reach';
trialtype(7).bool='&';

trialtype(8).outcomes=[nan nan nan nan nan nan nan nan nan 1 1 1 nan nan nan nan nan];
trialtype(8).name='paw_during_wheel';
trialtype(8).bool='|';

trialtype(9).outcomes=[nan nan nan nan nan nan nan nan nan nan nan nan 1 nan nan nan nan];
trialtype(9).name='chewing_at_trial_start';
trialtype(9).bool='&';

trialtype(10).outcomes=[nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan 1 nan];
trialtype(10).name='success_in_cued_window';
trialtype(10).bool='&';

trialtype(11).outcomes=[nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan 1 1];
trialtype(11).name='touch_in_cued_window';
trialtype(11).bool='|';

% Output settings
settings.trialStructure=trialStructure;
settings.timeWindow=timeWindow;
settings.bool_test=bool_test;
settings.trialtype=trialtype;

end