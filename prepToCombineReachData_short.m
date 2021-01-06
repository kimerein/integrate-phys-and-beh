% script for preparing directory structure before combining data across videos and days

% assign mouse and session numbers
expt_dir='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\Behavior Final Data Sets\GtACR2 batch opto';
optoOnHere=1; % 1 if there was opto, else 0
startMouseNumberingAt=1; % if want to start mouse numbers at some offset
% organiziation of data in this directory must be 
% top folder: by mouse name
% next folder: by session
% each folder within should contain a .txt with the date of the session,
% e.g., "20200203.txt"
% get a cell array containing the locations of processed_data folders
continuingAnalysisDir={};
k=1;
ls=dir(expt_dir);
for i=1:length(ls) % iterate through mice
    thisname=ls(i).name;
    thisisdir=ls(i).isdir;
    dates_db={};
    disp(thisname);
    if thisisdir==1 % must be the name of a mouse
        sub_ls=dir([expt_dir '\' thisname]);
        for j=1:length(sub_ls) % iterate through sessions
            sub_thisname=sub_ls(j).name;
            sub_thisisdir=sub_ls(j).isdir;
            if ~isempty(regexp(sub_thisname,'processed_data','ONCE')) && sub_thisisdir==1
                % add to list
                continuingAnalysisDir{k}=[expt_dir '\' thisname '\' sub_thisname];
                % set mouse id
                mouse_id=i-2+(startMouseNumberingAt-1);
                save([expt_dir '\' thisname '\' sub_thisname '\mouse_id.mat'],'mouse_id');
                % set session id
                sub_sub_ls=dir([expt_dir '\' thisname '\' sub_thisname]);
                for l=1:length(sub_sub_ls)
                    % look for date .txt file
                    if ismember(1,regexp(sub_sub_ls(l).name,'20'))
                        % assign from database
                        if ismember(sub_sub_ls(l).name,dates_db)
                            f=find(ismember(dates_db,sub_sub_ls(l).name));
                            nth_session=f;
                        else % or add to database
                            dates_db{length(dates_db)+1}=sub_sub_ls(l).name;
                            nth_session=length(dates_db);
                        end
                    end
                end
                save([expt_dir '\' thisname '\' sub_thisname '\nth_session.mat'],'nth_session');
                % save optoOnHere or not
                save([expt_dir '\' thisname '\' sub_thisname '\optoOnHere.mat'],'optoOnHere');
                k=k+1;
            end
        end
    end
end

checkForESPreaching(continuingAnalysisDir); % manually check for preemptive reaching