% script for preparing directory structure before combining data across videos and days

%% Make mouse database
% Step 1: Add all mice to makeMouseDatabase.m
% Step 2: Run makeMouseDatabase.m
% Notes:
%   Mouse database gathers knowledge about the directory structure that contains
%   all .AVI video files
%
%   DB_savedir and DB_filename are required
%
%   behLogFile is optional -- points to a .tsv file
%   
%   behLogFile is from the Google Sheet daily behavior log
%
%   each row of the sheet behLogFile MUST begin with a date (i.e., column 1
%   is experiment dates)
% 
%   existingDB is optional -- points to a previously constructed mouse
%   database
% 
%   pass in existingDB to maintain existing mouse >> mouse ID mapping;
%   otherwise, a new mapping is generated every time you run
%   makeMouseDatabase.m

DB_savedir='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\';
DB_filename='mouse_database';
DB_filename(~ismember(DB_filename,['A':'Z' 'a':'z' '0':'9']))='';
DB_filename=DB_filename(~isspace(DB_filename));
behLogFile='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\KER Behavior\By date\Behavior Log Starting 20190910 - Sheet1.tsv'; % make empty [] if unused
existingDB='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\20191106 batch opto update\mouse_database_20191106.mat'; % make empty [] if unused

db=makeMouseDatabase(behLogFile,existingDB);
save([DB_savedir DB_filename '.mat'],'db');

%% Add expt details to processed_data folders
% Add mouse ID, session numbers, optogenetic stimulation threshold, behavior 
% type (REGULAR or CONTROL), and whether mouse exhibited preemptive reaching
% to folders for further analysis

% load in mouse database
a=load([DB_savedir DB_filename '.mat']);
db=a.db;
% AND/OR
% combine mouse databases
% be sure that you now have the biggest, most inclusive mouse database, before assigning mouse IDs and session numbers 

% manually add all processed_data folders to one directory
% continuingAnalysisDir='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\20191106 batch opto update\interleaved and batch opto'; % this is the directory containing all processed_data folders for further analysis
% OR
% get all processed_data folders from .AVI locations
[continuingAnalysisDir,dbs_procData]=getProcessedData(db);

getMouseAndSessionNumbers(continuingAnalysisDir,db); % assign mouse ID and session numbers based on mouse database
getOptoThreshForExpts(continuingAnalysisDir,false,[]); % manually set optogenetic stimulation threshold for experiments
checkForESPreaching(continuingAnalysisDir); % manually check for preemptive reaching
checkForControl(continuingAnalysisDir); % use the parsedOutput.mat file to check for missing cues, suggesting a control experiment
if isfield(db,'behLog')
    getReachExptType(continuingAnalysisDir,dbs_procData,db.behLog); % use the behavior log to check for control experiments
end