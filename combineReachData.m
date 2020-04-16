% script for combining reach data across videos and days

%% set up filter (i.e., choose which experiments to include)

% load in mouse database
a=load([DB_savedir DB_filename 'wProcData.mat']);
db=a.db;

% choose filter parameters
filtStruct='db_bymouse';
filtParam='mouse_name';
useMice={'March_A',...
         'March_B',...
         'March_C'};
excludeControls=true;
     
 % filter mouse database    
db=filterMouseDatabase(db,filtStruct,filtParam,useMice,excludeControls);

%% get experiments

if isfield(db,'dbs_procData')
    % take processed_data folders from db
    continuingAnalysisDir=db.dbs_procData.dataFolders;
else
    % manually add all processed_data folders to one directory
    % continuingAnalysisDir='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\20191106 batch opto update\interleaved and batch opto'; % this is the directory containing all processed_data folders for further analysis
    % OR
    % get all processed_data folders from .AVI locations
    [continuingAnalysisDir,db.dbs_procData]=getProcessedData(db);
end

%% combine experiments

% settings for how to combine experiments
settings.check_for_human=1; % 1 if want to exclude data without humanChecked.txt
settings.discardPreemptive=true; % true if want to discard data sets with preemptive reaching
useAsCue='cueZone_onVoff'; % name of cue in tbt
cueDuration=0.25; % in seconds
doRealign=1; % 1 if want to realign data based on cue starts (i.e., align all cue onsets)

[alltbt,allmetadata]=combineExptPieces(continuingAnalysisDir,useAsCue,cueDuration,doRealign,settings);
% define and classify trial types
% modify settings in trialTypeSettings.m to change trial types
% modify settings in reachExpt_analysis_settings.m to change
% experiment-specific settings
[out,alltbt]=getSweepsFromBeh(alltbt);
metadata=howFarThroughSession(metadata);

% save field by field