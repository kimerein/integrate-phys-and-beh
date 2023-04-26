% script for combining reach data across videos and days

% note that this script takes about X mins to run for 55 days of data

%% combine experiments

temp=datestr(datetime('now'));
temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))=''; 
temp=temp(~isspace(temp));
saveDir=['\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\' 'alltbt' temp]; % save combined data to this directory
mkdir(saveDir);

% settings for how to combine experiments
settings=reachExpt_analysis_settings('display settings'); % modify settings in reachExpt_analysis_settings.m to change experiment-specific settings

% combine data
disp('Running combineExptPieces');
[alltbt,metadata]=combineExptPieces(continuingAnalysisDir,settings.nameOfCue,settings.cueDuration,settings.doRealign,settings);
% get rid of reaches during chewing
disp('Running clearReachesAfterSuccess');
alltbt=clearReachesAfterSuccess(alltbt,'all_reachBatch','reachBatch_success_reachStarts');
% make all_reachBatch field have only reaches NOT starting on wheel
disp('Running clearReachesFromPawOnWheel');
alltbt=clearReachesFromPawOnWheel(alltbt);
% add times with respect to session start
disp('Running getTimesWrtSessionStart');
alltbt=getTimesWrtSessionStart(alltbt,metadata);

% define and classify trial types
% modify settings in trialTypeSettings.m to change trial types
[out,alltbt]=getSweepsFromBeh(alltbt,settings);
metadata=howFarThroughSession(metadata,false,out);
[metadata,alltbt,out]=add_dprimes_to_tbt(alltbt,out,metadata,[],'all_reachBatch','cueZone_onVoff',settings);

% save field by field
mkdir([saveDir '\alltbt\']); mkdir([saveDir '\out\']); mkdir([saveDir '\metadata\']);
saveStructFieldByField(alltbt,[saveDir '\alltbt\']); % save alltbt
saveStructFieldByField(out,[saveDir '\out\']); % save out
saveStructFieldByField(metadata,[saveDir '\metadata\']); % save metadata
save([saveDir '\reachExptAnalysis_settings.mat'],'settings'); % save reach expt analysis settings