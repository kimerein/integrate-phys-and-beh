% script for combining reach data across videos and days

%% set up filter (i.e., choose which experiments to include)

% load in mouse database
a=load([DB_savedir DB_filename '.mat']);
db=a.db;

% filter mouse database


% get processed_data folders in database
[continuingAnalysisDir,dbs_procData]=getProcessedData(db);


%% get experiments



% manually add all processed_data folders to one directory
% continuingAnalysisDir='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\20191106 batch opto update\interleaved and batch opto'; % this is the directory containing all processed_data folders for further analysis
% OR
% get all processed_data folders from .AVI locations
[continuingAnalysisDir,dbs_procData]=getProcessedData(db);


%% combine experiments



%% save combined data set and filter

[alltbt,metadata]=combineExptPieces('\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\mouse summary data\bi_both','cueZone_onVoff',0.25,1);

[out,alltbt]=getSweepsFromBeh(alltbt);

metadata=howFarThroughSession(metadata);