function settings=trialTypeSettings()

% This function defines types of behavioral trials
% Modify this function to change trial-type classification scheme

% "Reach batch" definition
reach_batch.window=0.3; % in seconds, reaches must occur within this many seconds of each other to be in same batch
reach_batch.firstreach_type={'miss_reachStarts'}; % first reach must be one of these types
reach_batch.secondreach_type={'success_reachStarts_pawOnWheel','drop_reachStarts_pawOnWheel','miss_reachStarts_pawOnWheel'}; % second reach must be one of these types
reach_batch.take_first_or_second_type=2; % number indicates whether to convert batch to the type of the first or second reach

% Test whether paw was on wheel while wheel turning
bool_test(1).testwhat='single reach';
bool_test(1).fieldname='pawOnWheel';
bool_test(1).test='any';
bool_test(1).thresh=0.5;
bool_test(1).comparator='>';
bool_test(1).window='wheel_turning';

% Result of first "reach batch" after cue
% Note that "reach batch" refers to the outcome of the last reach in set of
% contiguous reaches
bool_test(2).testwhat='reach batch';
bool_test(2).fieldname='success_reachStarts';
bool_test(2).test='first';
bool_test(2).thresh=0.5;
bool_test(2).comparator='>';
bool_test(2).window='after cue';

bool_test(3).testwhat='reach batch';
bool_test(3).fieldname='drop_reachStarts';
bool_test(3).test='first';
bool_test(3).thresh=0.5;
bool_test(3).comparator='>';
bool_test(3).window='after cue';

bool_test(4).testwhat='reach batch';
bool_test(4).fieldname='miss_reachStarts';
bool_test(4).test='first';
bool_test(4).thresh=0.5;
bool_test(4).comparator='>';
bool_test(4).window='after cue';

bool_test(5).testwhat='reach batch';
bool_test(5).fieldname='pelletmissingreach_reachStarts';
bool_test(5).test='first';
bool_test(5).thresh=0.5;
bool_test(5).comparator='>';
bool_test(5).window='after cue';

% Put together bool_test results to get trial classifications
% Each index into trialtype(i).outcomes vector of length n refers to outcome of bool_tests 1:n
trialtype(1).outcomes=[0 1 0 0 0];
trialtype(1).name='success';
trialtype(2).outcomes=[0 0 1 0 0];
trialtype(2).name='drop';
trialtype(3).outcomes=[0 0 0 1 0];
trialtype(3).name='miss';
trialtype(4).outcomes=[0 0 0 0 1];
trialtype(4).name='no pellet';
trialtype(5).outcomes=[1 1 0 0 0];
trialtype(5).name='on wheel success';
trialtype(6).outcomes=[1 0 1 0 0];
trialtype(6).name='on wheel drop';
trialtype(7).outcomes=[1 0 0 1 0];
trialtype(7).name='on wheel miss';
trialtype(8).outcomes=[1 0 0 0 1];
trialtype(8).name='on wheel no pellet';

settings.reach_batch=reach_batch;
settings.bool_test=bool_test;
settings.trialtype=trialtype;

end