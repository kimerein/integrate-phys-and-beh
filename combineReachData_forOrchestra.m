% script for combining reach data across videos and days -- FOR ORCHESTRA

%% combine experiments

load('settings.mat');
load('continuingAnalysisDir');

temp=datestr(datetime('now'));
temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))=''; 
temp=temp(~isspace(temp));
saveDir=['alltbt' temp]; % save combined data to this directory
mkdir(saveDir);

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
[out,alltbt]=getSweepsFromBeh(alltbt);
metadata=howFarThroughSession(metadata);
[metadata,alltbt,out]=add_dprimes_to_tbt(alltbt,out,metadata,'all_reachBatch','cueZone_onVoff',settings);

% save field by field
saveStructFieldByField(alltbt,saveDir); % save alltbt
saveStructFieldByField(out,saveDir); % save out
saveStructFieldByField(metadata,saveDir); % save metadata
save([saveDir '\reachExptAnalysis_settings.mat'],'settings'); % save reach expt analysis settings