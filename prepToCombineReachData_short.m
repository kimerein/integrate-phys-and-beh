% script for preparing directory structure before combining data across videos and days

% assign mouse and session numbers
expt_dir='Z:\Kim\Behavior Final Data Sets\Controls for learning curves';
optoOnHere=1; % 1 if there was opto, else 0
onlySaveOptoIfDoesNotExist=true;
startMouseNumberingAt=1; % if want to start mouse numbers at some offset
% organization of data in this directory must be 
% top folder: by mouse name
% next folder: by session
% each folder within should contain a .txt with the date of the session,
% e.g., "20200203.txt"
% get a cell array containing the locations of processed_data folders
continuingAnalysisDir={};
k=1;
ls=dir(expt_dir);
mouse_database.names={};
for i=1:length(ls) % iterate through mice
    thisname=ls(i).name;
    thisisdir=ls(i).isdir;
    dates_db={};
    disp(thisname);
    if thisisdir==1 % must be the name of a mouse
        sub_ls=dir([expt_dir '\' thisname]);
        if ~any(strcmp(mouse_database.names,thisname))
            % add this mouse to database
            currmousenameind=length(mouse_database.names)+1;
            mouse_database.names{currmousenameind}=thisname;
            notyetsavedmouseid=true;
        end
        for j=1:length(sub_ls) % iterate through sessions
            sub_thisname=sub_ls(j).name;
            sub_thisisdir=sub_ls(j).isdir;
            if ~isempty(regexp(sub_thisname,'processed_data','ONCE')) && sub_thisisdir==1
                % add to list
                continuingAnalysisDir{k}=[expt_dir '\' thisname '\' sub_thisname];
                % set mouse id
                mouse_id=i-2+(startMouseNumberingAt-1);
                save([expt_dir '\' thisname '\' sub_thisname '\mouse_id.mat'],'mouse_id');
                if notyetsavedmouseid==true
                    mouse_database.mouse_id(currmousenameind)=mouse_id;
                    notyetsavedmouseid=false;
                end
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
                if onlySaveOptoIfDoesNotExist==true
                    if isfile([expt_dir '\' thisname '\' sub_thisname '\optoOnHere.mat'])
                    else
                        save([expt_dir '\' thisname '\' sub_thisname '\optoOnHere.mat'],'optoOnHere');
                    end
                else
                    save([expt_dir '\' thisname '\' sub_thisname '\optoOnHere.mat'],'optoOnHere');
                end
                k=k+1;
            end
        end
    end
end

checkForESPreaching(continuingAnalysisDir); % manually check for preemptive reaching