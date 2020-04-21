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

temp=datestr(datetime('now'));
temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))=''; 
temp=temp(~isspace(temp));
saveDir=['\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\' 'alltbt' temp]; % save combined data to this directory

% settings for how to combine experiments
settings=reachExpt_analysis_settings('display settings'); % modify settings in reachExpt_analysis_settings.m to change experiment-specific settings

% combine data
[alltbt,metadata]=combineExptPieces(continuingAnalysisDir,settings.nameOfCue,settings.cueDuration,settings.doRealign,settings);
% define and classify trial types
% modify settings in trialTypeSettings.m to change trial types
[out,alltbt]=getSweepsFromBeh(alltbt);
metadata=howFarThroughSession(metadata);
[metadata,alltbt]=add_dprimes_to_tbt(alltbt,out,metadata,'all_reachBatch','cueZone_onVoff',settings);

% save field by field
saveStructFieldByField(alltbt,saveDir); % save alltbt
saveStructFieldByField(out,saveDir); % save out
saveStructFieldByField(metadata,saveDir); % save metadata
save([saveDir '\reachExptAnalysis_settings.mat'],'settings'); % save reach expt analysis settings
save([saveDir '\mouse_database.mat'],'db'); % save filtered mouse database