function db=makeMouseDatabase()

% Mouse Database should be a structure with the following fields
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

function assignIDnums

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

data_directories=db(f).data_directories;
low_speed_videos=db(f).low_speed_videos;

data_directories{length(data_directories)+1}=addDir;

ls=dir(addDir);
for i=1:length(ls)
    thisname=ls(i).name;
    if ~isempty(regexp(thisname,'.AVI','ONCE')) % is a low-speed video file
        % add to list
        low_speed_videos{length(low_speed_videos)+1}=[addDir '\' thisname];
    end
end

db(f).data_directories=data_directories;
db(f).low_speed_videos=low_speed_videos;

end