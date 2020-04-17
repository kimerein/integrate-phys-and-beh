% runReactionTimeAnalysis.m

% script for running a frequently used subset of analyses

%% load in data

exptDataDir='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\alltbt'; % directory containing experimental data

alltbt=loadStructFieldByField([exptDataDir '\alltbt']); % load alltbt
out=loadStructFieldByField([exptDataDir '\out']); % load out
metadata=loadStructFieldByField([exptDataDir '\metadata']); % load metadata
a=load([exptDataDir '\reachExptAnalysis_settings.mat']); % load reach expt analysis settings 
reachExptSettings=a.settings;

%% choose additional settings for reaction time analysis

settings=RTanalysis_settings('display settings');

%% build relevant data sets

% buildReachingRTModel.m
[dataset,correctedDistributions]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir);

% save datasets

%% fit model terms

% putTogetherTermsForRTModel.m

%% plot

% plotBehaviorEventFx.m

% fitRTMtermToOtherReaches.m

% fitLocalCorrelationStructureOfReaches.m

% fitReachingRTModel.m

% pelletTouchAlignedReaching.m

