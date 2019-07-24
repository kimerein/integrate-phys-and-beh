function outdbs=makeMouseDatabase()

% Assumption about directory structure:
% AVIs from the same session should be in the same directory!

% Mouse Database should be a structure with the following fields
% 'mouse_name'          mouse's name
% 'mouse_id'            mouse's ID number
% 'description'         description of mouse genotype, injections and
%                       experimental goals
% 'trainers'            who trained this mouse
% 'cue start date'      first date of cue training on rig
% 'data_directories'    a list of all directories containing behavior data
%                       related to this mouse
% 'low_speed_videos'    a list of all low-speed videos (.AVI) with behavior
%                       data for this mouse
% 'directory_per_video' the directory in which each low-speed video file
%                       is stored

db=[];

% March_A
mn='March_A';
db=addNameAndDescription(db,mn,'NkxCre X ReaChR, batch opto during learning');
db=addTrainer(db,mn,'Kim');
db=addTrainer(db,mn,'Marci');
db=addDirectoriesRecursive(db,mn,'\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\KER Behavior',false);

% March_B
mn='March_B';
db=addNameAndDescription(db,mn,'NkxCre X ReaChR, batch opto during learning');
db=addTrainer(db,mn,'Kim');
db=addTrainer(db,mn,'Marci');
db=addDirectoriesRecursive(db,mn,'\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\KER Behavior',false);

% March_C
mn='March_C';
db=addNameAndDescription(db,mn,'NkxCre X ReaChR, batch opto during learning');
db=addTrainer(db,mn,'Kim');
db=addTrainer(db,mn,'Marci');
db=addDirectoriesRecursive(db,mn,'\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\KER Behavior',false);

% March_D
mn='March_D';
db=addNameAndDescription(db,mn,'NkxCre X ReaChR, batch opto during learning');
db=addTrainer(db,mn,'Kim');
db=addTrainer(db,mn,'Marci');
db=addDirectoriesRecursive(db,mn,'\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\KER Behavior',false);

% Feb_1
mn='Feb_1';
db=addNameAndDescription(db,mn,'NkxCre X ReaChR, continuous opto during learning');
db=addTrainer(db,mn,'Kim');
db=addTrainer(db,mn,'Marci');
db=addDirectoriesRecursive(db,mn,'\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\KER Behavior',false);

% Feb_2
mn='Feb_2';
db=addNameAndDescription(db,mn,'NkxCre X ReaChR, continuous opto during learning');
db=addTrainer(db,mn,'Kim');
db=addTrainer(db,mn,'Marci');
db=addDirectoriesRecursive(db,mn,'\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\KER Behavior',false);

% Feb_3
mn='Feb_3';
db=addNameAndDescription(db,mn,'control genotype, no str silencing, continuous opto during learning');
db=addTrainer(db,mn,'Kim');
db=addTrainer(db,mn,'Marci');
db=addDirectoriesRecursive(db,mn,'\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\KER Behavior',false);

% Feb_4
mn='Feb_4';
db=addNameAndDescription(db,mn,'control genotype, no str silencing, continuous opto during learning');
db=addTrainer(db,mn,'Kim');
db=addTrainer(db,mn,'Marci');
db=addDirectoriesRecursive(db,mn,'\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\KER Behavior',false);

% Jan2_gray

% Jan3_gray
mn='Jan3_gray';
db=addNameAndDescription(db,mn,'NkxCre X ReaChR, interleaved opto during learning');
db=addTrainer(db,mn,'Kim');
db=addTrainer(db,mn,'Marci');
db=addDirectoriesRecursive(db,mn,'\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\KER Behavior',false);

% Jan3_white

% 3F_white
mn='3F_white';
db=addNameAndDescription(db,mn,'NkxCre X ReaChR, interleaved opto during learning');
db=addTrainer(db,mn,'Kim');
db=addTrainer(db,mn,'JW');
db=addDirectoriesRecursive(db,mn,'\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\KER Behavior',false);
db=addDirectoriesRecursive(db,mn,'\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\JW\Behavior',false);

% 3F_spot
mn='3F_spot';
db=addNameAndDescription(db,mn,'NkxCre X ReaChR, continuous opto during learning');
db=addTrainer(db,mn,'Kim');
db=addTrainer(db,mn,'JW');
db=addDirectoriesRecursive(db,mn,'\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\KER Behavior',false);
db=addDirectoriesRecursive(db,mn,'\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\JW\Behavior',false);

% 2F_white
mn='2F_white';
db=addNameAndDescription(db,mn,'NkxCre X ReaChR, interleaved opto during learning');
db=addTrainer(db,mn,'Kim');
db=addTrainer(db,mn,'JW');
db=addDirectoriesRecursive(db,mn,'\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\KER Behavior',false);
db=addDirectoriesRecursive(db,mn,'\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\JW\Behavior',false);

% bi_agouti
mn='bi_agouti';
db=addNameAndDescription(db,mn,'NkxCre X ReaChR, interleaved opto during learning');
db=addTrainer(db,mn,'Julia G');
db=addDirectoriesRecursive(db,mn,'\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\JuliaG\Behavior Expts',false);

% bi_bl
mn='bi_bl';
db=addNameAndDescription(db,mn,'NkxCre X ReaChR, interleaved opto during learning');
db=addTrainer(db,mn,'Julia G');
db=addDirectoriesRecursive(db,mn,'\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\JuliaG\Behavior Expts',false);

% bi_both
mn='bi_both';
db=addNameAndDescription(db,mn,'NkxCre X ReaChR, interleaved opto during learning');
db=addTrainer(db,mn,'Julia G');
db=addDirectoriesRecursive(db,mn,'\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\JuliaG\Behavior Expts',false);

% str_agouti
mn='str_agouti';
db=addNameAndDescription(db,mn,'NkxCre X ReaChR, interleaved opto after learning');
db=addTrainer(db,mn,'Julia G');
db=addTrainer(db,mn,'Whitney');
db=addDirectoriesRecursive(db,mn,'\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\JuliaG\Behavior Expts',false);
db=addDirectoriesRecursive(db,mn,'\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\Whitney\Behavior Experiments',false);

% str_bl
mn='str_bl';
db=addNameAndDescription(db,mn,'NkxCre X ReaChR, interleaved opto after learning');
db=addTrainer(db,mn,'Julia G');
db=addTrainer(db,mn,'Whitney');
db=addDirectoriesRecursive(db,mn,'\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\JuliaG\Behavior Expts',false);
db=addDirectoriesRecursive(db,mn,'\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\Whitney\Behavior Experiments',false);

% str_wh
mn='str_wh';
db=addNameAndDescription(db,mn,'NkxCre X ReaChR, interleaved opto after learning');
db=addTrainer(db,mn,'Julia G');
db=addTrainer(db,mn,'Whitney');
db=addDirectoriesRecursive(db,mn,'\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\JuliaG\Behavior Expts',false);
db=addDirectoriesRecursive(db,mn,'\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\Whitney\Behavior Experiments',false);

% assign mouse IDs and others
[db,dbs]=assignIDnums(db);

outdbs.db_bymouse=db;
outdbs.dbs=dbs;

end

function db=addNameAndDescription(db,name,description)

if isempty(db)
    i=1;
else
    i=length(db)+1;
end

db(i).mouse_name=name;
db(i).description=description;

disp(['Adding ' name ' to database']);

end

function [db,dbs]=assignIDnums(db)

notdefined=true;
for i=1:length(db)
    % for each mouse, assign a new ID
    db(i).mouse_id=i;
    % for each mouse, get the unique sessions
    % AVIs from the same session should be in the same directory!
    % [~,~,sessid]=unique(db(i).directory_per_video);
    % db(i).sessids=sessid;
    
    % for each mouse, sort sessions according to date
    [~,ui]=unique(db(i).date_per_video);
    [~,si]=sort(db(i).date_per_video);
    clear unique_sessids
    unique_sessids(si)=1:length(db(i).date_per_video);
    for j=1:length(ui)
        currdatename=db(i).date_per_video(ui(j));
        for k=1:length(db(i).date_per_video)
            if strcmp(db(i).date_per_video(k),currdatename)==1
                % this is a duplicate
                unique_sessids(k)=unique_sessids(ui(j));
            end
        end
    end 
    ordered_sessids=nan(size(unique_sessids));
    k=1;
    for j=1:length(unique_sessids)
        if any(unique_sessids==j)
            ordered_sessids(unique_sessids==j)=k;
            k=k+1;
        end
    end 
    % rename sessids for this mouse to match order in which sessions
    % occurred
    db(i).sessids=ordered_sessids;
    if notdefined==true
        dbs.vids_to_match_mouseIDs(1:1+length(db(i).low_speed_videos)-1)=db(i).low_speed_videos;
        dbs.mouseIDs_to_match_vids(1:1+length(db(i).low_speed_videos)-1)=db(i).mouse_id;
        dbs.sessIDs_to_match_vids(1:1+length(db(i).low_speed_videos)-1)=db(i).sessids;
        dbs.sessDates_to_match_vids(1:1+length(db(i).low_speed_videos)-1)=db(i).date_per_video;
        notdefined=false;
    else
        dbs.vids_to_match_mouseIDs(length(dbs.vids_to_match_mouseIDs)+1:length(dbs.vids_to_match_mouseIDs)+1+length(db(i).low_speed_videos)-1)=db(i).low_speed_videos;
        dbs.mouseIDs_to_match_vids(length(dbs.mouseIDs_to_match_vids)+1:length(dbs.mouseIDs_to_match_vids)+1+length(db(i).low_speed_videos)-1)=db(i).mouse_id;
        dbs.sessIDs_to_match_vids(length(dbs.sessIDs_to_match_vids)+1:length(dbs.sessIDs_to_match_vids)+1+length(db(i).low_speed_videos)-1)=db(i).sessids;
        dbs.sessDates_to_match_vids(length(dbs.sessDates_to_match_vids)+1:length(dbs.sessDates_to_match_vids)+1+length(db(i).low_speed_videos)-1)=db(i).date_per_video;
    end
end

end

function db=addTrainer(db,mouse_name,addTrainer)

% list of trainers
t={ 'Kim'; ...
    'Marci'; ...
    'Mitchell'; ...
    'Julia G'; ...
    'JW'; ...
    'Whitney'; ...
    'Ian'};

if ~any(ismember(t,addTrainer))
    disp('Problem: do not recognize trainer');
end

% find which index into db
for i=1:length(db)
    mouse_names{i}=db(i).mouse_name;
end
f=find(ismember(mouse_names,mouse_name));

if isfield(db(f),'trainers')
    db(f).trainers{length(db(f).trainers)+1}=addTrainer;
else
    db(f).trainers={addTrainer};
end

end

function [db,mouse_name,addDir,inRightMouse]=addDirectoriesRecursive(db,mouse_name,addDir,inRightMouse)

% search for all directories in this directory with name matching
% mouse_name
ls=dir(addDir);
for i=1:length(ls)
    thisname=ls(i).name;
    thisisdir=ls(i).isdir;
    if thisisdir==1 && ~isempty(regexp(thisname,mouse_name,'ONCE')) % matches current mouse
        % either this is the bottom directory that contains the .AVI videos   
        % or this directory contains dates matching this mouse
        [db,mouse_name]=addDirectoriesRecursive(db,mouse_name,[addDir '\' thisname],true); 
    elseif thisisdir==1 % but does not match current mouse
        if strcmp(thisname,'.') || strcmp(thisname,'..')
            continue
        end
        [db,mouse_name]=addDirectoriesRecursive(db,mouse_name,[addDir '\' thisname],inRightMouse);
    end
end

if inRightMouse==true
    db=addDirectoryData(db,mouse_name,addDir); % add .AVI files from current directory
end

end

function db=addDirectoryData(db,mouse_name,addDir)

% add all .AVI files in addDir
% and add addDir

% find which index into db
for i=1:length(db)
    mouse_names{i}=db(i).mouse_name;
end
f=find(ismember(mouse_names,mouse_name));

if ~isfield(db(f),'data_directories')
    db(f).data_directories=[];
end
if ~isfield(db(f),'low_speed_videos')
    db(f).low_speed_videos=[];
end
if ~isfield(db(f),'directory_per_video')
    db(f).directory_per_video=[];
end
if ~isfield(db(f),'date_per_video')
    db(f).date_per_video=[];
end

data_directories=db(f).data_directories;
low_speed_videos=db(f).low_speed_videos;
directory_per_video=db(f).directory_per_video;
date_per_video=db(f).date_per_video;

data_directories{length(data_directories)+1}=addDir;

ls=dir(addDir);
for i=1:length(ls)
    thisname=ls(i).name;
    if ~isempty(regexp(thisname,'.AVI','ONCE')) % is a low-speed video file
        temp=getDateFromFilename(addDir);
        if isempty(temp) % no date associated with this .AVI file
            continue
        end
        currind=length(low_speed_videos)+1;
        date_per_video{currind}=temp;
        % add to list
        low_speed_videos{currind}=[addDir '\' thisname];
        directory_per_video{currind}=addDir;
    end
end

db(f).data_directories=data_directories;
db(f).low_speed_videos=low_speed_videos;
db(f).directory_per_video=directory_per_video;
db(f).date_per_video=date_per_video;

end

function date=getDateFromFilename(fn)

date=[];

r=regexp(fn,'2019','ONCE');
if ~isempty(r)
    date=fn(r:r+7);
    return
end
r=regexp(fn,'2018','ONCE');
if ~isempty(r)
    date=fn(r:r+7);
    return
end
r=regexp(fn,'2017','ONCE');
if ~isempty(r)
    date=fn(r:r+7);
    return
end
r=regexp(fn,'2016','ONCE');
if ~isempty(r)
    date=fn(r:r+7);
    return
end

end