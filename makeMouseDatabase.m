function outdbs=makeMouseDatabase(varargin)

% Assumption about directory structure:
% AVIs from the same session should be in the same directory!

% Assumption about the behavior log file .tsv
% Date, formatted as 6/4/1987 indicates the beginning of a row

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

if ~isempty(varargin)
    behLogFile=varargin{1};
    if length(varargin)>1
        % pass in existing mouse database, if want to keep existing mouse
        % to ID mappings
        temp=varargin{2};
        T_mouseToID=extractMouseNametoIDmapping(temp.db_bymouse);
    else
        T_mouseToID=[];
    end
else
    behLogFile=[];
    T_mouseToID=[];
end

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

% April_white
mn='April_white';
db=addNameAndDescription(db,mn,'A2a-Cre x D1-Cre, batch opto during learning');
db=addTrainer(db,mn,'Marci');
db=addDirectoriesRecursive(db,mn,'\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\KER Behavior',false);

% May
mn='May';
db=addNameAndDescription(db,mn,'A2a-Cre x D1-Cre, batch opto during learning');
db=addTrainer(db,mn,'Marci');
db=addDirectoriesRecursive(db,mn,'\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\KER Behavior',false);

% April_short
mn='April_short';
db=addNameAndDescription(db,mn,'NkxCre X ReaChR, batch opto during learning');
db=addTrainer(db,mn,'Marci');
db=addDirectoriesRecursive(db,mn,'\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\KER Behavior',false);

% April_long
mn='April_long';
db=addNameAndDescription(db,mn,'NkxCre X ReaChR, batch opto during learning');
db=addTrainer(db,mn,'Marci');
db=addDirectoriesRecursive(db,mn,'\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\KER Behavior',false);

% assign mouse IDs and others
[db,dbs]=assignIDnums(db,T_mouseToID);

% Optional: get more details about experiments from behavior log file
if ~isempty(behLogFile)
    behLog=readBehLogFile(behLogFile);
    [behLog,dbs]=behLogByVid(behLog,db,dbs);
end

outdbs.db_bymouse=db;
outdbs.dbs=dbs;
if ~isempty(behLogFile)
    outdbs.behLog=behLog;
end

end

function T=extractMouseNametoIDmapping(db_bymouse)

mapping=cell(length(db_bymouse),2); % first col is mouse name, second col is mouse ID
for i=1:length(db_bymouse)
    mapping{i,1}=db_bymouse(i).mouse_name;
    mapping{i,2}=db_bymouse(i).mouse_id;
end

T=cell2table(mapping,'VariableNames',{'name','id'});

end

function [behLog,dbs]=behLogByVid(behLog,db_bymouse,dbs)

% for each row in behLog, add a pointer to the relevant indices into
% vids_to_match_mouseIDs, and vice versa

for i=1:length(db.dbs.sessDates_to_match_vids)
    currstr=db.dbs.sessDates_to_match_vids{i};
    % assume yyyymmdd format in string
    dat2(i)=datetime(str2double(currstr(1:4)),str2double(currstr(5:6)),str2double(currstr(7:8)));
end

% convert mouse names in behLog to mouseIDs consistent with dbs
T_mouseToID=extractMouseNametoIDmapping(db_bymouse);
behLog.ID=nan(1,length(behLog.Mouse));
namesToIgnore={};
for i=1:length(behLog.Mouse)
    indIntoID=find(strcmp(T_mouseToID.name,behLog.Mouse(i)));
    if isempty(indIntoID)
        % did not recognize mouse name in behLog
        % check if this is a mouse to ignore
        if any(strcmp(namesToIgnore,behLog.Mouse(i)))
            continue
        end
        % otherwise, ask user
        tit='Problem: Unrecognized mouse name';
        promp={['Mouse name in behavior log is ' behLog.Mouse(i) '.'],'Please select the matching mouse name from this list.'};
        temp=T_mouseToID.name;
        temp{length(temp)+1}='None of the above';
        indx=listdlg('ListString',temp,'PromptString',promp,'Name',tit,'SelectionMode','single');
        if isempty(indx)
            % user canceled
            error('Unrecognized mouse name in behavior log');
        end
        if indx==length(temp)
            % none of the above
            answer=questdlg('Do you want to completely ignore this mouse in further analysis?');
            switch answer
                case 'Yes'
                    namesToIgnore{length(namesToIgnore)+1}=behLog.Mouse(i);
                case 'No'
                    error('Please add this mouse in makeMouseDatabase.m.');
                otherwise
                    error('Please add this mouse in makeMouseDatabase.m.');
            end
        else
            % add this mouse name as a synonym
            T_mouseToID(size(T_mouseToID,1)+1,:)={behLog.Mouse(i),T_mouseToID.id(indx)};
            % get mouse ID
            behLog.ID(i)=T_mouseToID.id(indx);
        end
    else
        % get mouse ID
        behLog.ID(i)=T_mouseToID.id(indIntoID);
    end
end

% Match mouse name and date for each video
for i=1:size(behLog,1)
    dat1=behLog.Date(i);
    id1=behLog.ID(i);
    
    indsInto_dbs=find(dbs.mouseIDs_to_match_vids==id1 & ismember(dat2,dat1));
    behLog.indsInto_vids(i)={indsInto_dbs};
end

dbs.indsInto_behLog=nan(size(dbs.mouseIDs_to_match_vids));
for i=1:length(dbs.mouseIDs_to_match_vids)
    for j=1:length(behLog.indsInto_vids)
        temp=behLog.indsInto_vids(j);
        temp=temp{1};
        if any(ismember(temp,i))
            dbs.indsInto_behLog(i)=j;
            break
        end
    end
end

end

function T=readBehLogFile(behLogFile)

approxColumns=26; % approximate number of columns for table pre-allocation

% Check that is .tsv
fi=regexp(behLogFile,'\.');
if ~strcmp(behLogFile(fi(end):end),'.tsv')
    error('Expected behLogFile to be .tsv format');
end
fid=fopen(behLogFile);
str=textscan(fid,'%s','Delimiter','\t');
s=str{1};
expression = ['(?<month>\d+)/(?<day>\d+)/(?<year>\d+)|','(?<day>\d+)-(?<month>\d+)-(?<year>\d+)'];
k=1;
stillInHeader=true;
for i=1:length(s)
    % find dates -- should be row beginnings
    tokenNames=regexp(s{i},expression,'names');
    % assume that everything before first date is part of the header
    if ~isempty(tokenNames)
        % this is a date
        if k==1 
            % this is the first date
            % everything before this is header
            header=s(1:i-1);
            % make table
            sz=[ceil(length(s)/approxColumns)+100 length(header)];
            C=cell(sz);
            varTypes{1}='datetime';
            for j=2:length(s)+10
                varTypes{j}='string';
            end
            for m=1:length(header)
                if isempty(header{m})
                    header{m}=['Empty Col ' num2str(m)];
                end
            end
            stillInHeader=false;
            % make first row
            t=datetime(str2num(tokenNames.year),str2num(tokenNames.month),str2num(tokenNames.day));
            l=1; % l is column
            C{k,l}=t;
            l=l+1;
            k=k+1; % k is row
        else
            % start a new row
            t=datetime(str2num(tokenNames.year),str2num(tokenNames.month),str2num(tokenNames.day));
            l=1; % l is column
            C{k,l}=t;
            l=l+1;
            k=k+1; % k is row
        end     
    else
        if stillInHeader==true
        else
            C{k,l}=s{i};
            l=l+1;
        end
    end
end

if length(header)<size(C,2)
    % add to header
    for m=length(header)+1:size(C,2)
        header{m}=['Empty Col ' num2str(m)];
    end
end
for i=1:length(header)
    temp=header{i};
    temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))=''; 
    header{i}=temp(~isspace(temp));
end
T=cell2table(C,'VariableNames',header);

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

function [db,dbs]=assignIDnums(db,existIDs)

notdefined=true;
currAddedMouse=1;
for i=1:length(db)
    % for each mouse, assign a new ID
    if ~isempty(existIDs)
        indx=find(strcmp(existIDs.name,db(i).mouse_name));
        if ~isempty(indx)
            db(i).mouse_id=existIDs.id(indx);
        else
            db(i).mouse_id=nanmax(existIDs.id)+currAddedMouse;
            currAddedMouse=currAddedMouse+1;
        end
    else
        db(i).mouse_id=i;
    end
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

r=regexp(fn,'2020','ONCE');
if ~isempty(r)
    date=fn(r:r+7);
    return
end
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