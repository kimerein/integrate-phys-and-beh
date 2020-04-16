function [dataFolders,dbs_matchDataFolders]=getProcessedData(db)

% using the mouse database, find all corresponding processed_data folders
% Assumption: corresponding processed_data folder should be in the same
% directory as the .AVI file, for example,
% \\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\KER Behavior\3F_white\20190306\O2 output\2019-03-06 16-16-03-C_processed_data
% and \\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\KER
% Behavior\3F_white\20190306\O2 output\2019-03-06 16-16-03-C.AVI

dataFolders={};
usedWhichVids=[];
for i=1:length(db.dbs.vids_to_match_mouseIDs)
    currVid=db.dbs.vids_to_match_mouseIDs{i};
    endOfDir=regexp(currVid,'\');
    currDir=currVid(1:endOfDir(end));
    justVid=currVid(endOfDir(end)+1:end-4); % drop file extension
    % look for processed_data folder
    if exist([currDir justVid '_processed_data'],'dir')
        dataFolders{length(dataFolders)+1}=[currDir justVid '_processed_data'];
        usedWhichVids(length(usedWhichVids)+1)=i;
    end
end

f=fieldnames(db.dbs);
for i=1:length(f)
    temp=db.dbs.(f{i});
    dbs_matchDataFolders.(f{i})=temp(usedWhichVids);
end
dbs_matchDataFolders.usedWhichVids=usedWhichVids;

% Get rid of duplicates
% Take the first instance of any duplicated processed_data folders
[dataFolders,si]=unique(dataFolders);
f=fieldnames(dbs_matchDataFolders);
for i=1:length(f)
    temp=dbs_matchDataFolders.(f{i});
    dbs_matchDataFolders.(f{i})=temp(si);
end

dbs_matchDataFolders.dataFolders=dataFolders;
end