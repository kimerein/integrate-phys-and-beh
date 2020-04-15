function getMouseAndSessionNumbers(expt_dir,mouseDatabase)

% get mouse ID database
if isempty(mouseDatabase)
    mouseDatabase=makeMouseDatabase;
elseif isstr(mouseDatabase)
    a=load(mouseDatabase);
    mouseDatabase=a.mouseDatabase;
else
    % passed in mouseDatabase
end
    
% this should be a structure with the following fields
% 'mouse_name'          mouse's name
% 'mouse_id'            mouse's ID number
% 'description'         description of mouse genotype, injections and
%                       experimental goals
% 'trainers'            who trained this mouse
% 'cue start date'      first date of cue training on rig
% 'sac date'            date mouse was sac'ced
% 'data_directories'    a list of all directories containing behavior data
%                       related to this mouse
% 'low_speed_videos'    a list of all low-speed videos (.AVI) with behavior
%                       data for this mouse

% use this list of data_directories and low_speed_videos associated with
% each mouse to figure out mouse's id for analyzed data

for i=1:length(mouseDatabase.dbs.vids_to_match_mouseIDs)
    temp=mouseDatabase.dbs.vids_to_match_mouseIDs{i};
    rout=regexp(temp,'\');
    vidnames{i}=temp(rout(end)+1:end);
end

if ~iscell(expt_dir)
    ls=dir(expt_dir);
else
    ls=expt_dir;
end
for i=1:length(ls)
    if ~iscell(expt_dir)
        thisname=ls(i).name;
        thisisdir=ls(i).isdir;
    else
        currdir=ls{i};
        temp=regexp(currdir,'\');
        thisname=currdir(temp(end)+1:end);
        thisisdir=isempty(regexp(thisname,'\.','ONCE'));
    end
    r=regexp(thisname,'processed_data','ONCE');
    if ~isempty(r) && thisisdir==1
        avi_name=thisname(1:r-2);
        % look for which mouse has this AVI video
        % and in which session
        f=find(ismember(vidnames,[avi_name '.AVI']));
        mouse_id=mouseDatabase.dbs.mouseIDs_to_match_vids(f);
        nth_session=mouseDatabase.dbs.sessIDs_to_match_vids(f);
        if length(mouse_id)>1
            if any(diff(mouse_id))~=0
                error('Problem! Multiple mice associated with this video file');
            end
            mouse_id=mouse_id(1);
        end
        if length(nth_session)>1
            if any(diff(nth_session))~=0
                error('Problem! Multiple dates associated with this video file');
            end
            nth_session=nth_session(1);
        end
        save([expt_dir '\' thisname '\mouse_id.mat'],'mouse_id');
        save([expt_dir '\' thisname '\nth_session.mat'],'nth_session');
    end
end