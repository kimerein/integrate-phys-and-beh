function getNthSession(dataTable,mouse_database,metadata,alsoFixOptoOnHere,excludeTrainingRig)

% Use TAB delimiter when downloading table
    
% fix nth_session in metadata
% for days/mice missing from this table, leave nth_session alone

% Structure of mouse_database needs to be:
% mouse_database.names{currmousenameind} is name of mouse
% mouse_database.mouse_id(currmousenameind) is mouseid

% Fields of metadata must include dateFromTextFile, which contains a
% datetime

% Fields of Behavior Log are:
% Date
% Mouse
% Rig
% Order
% Start
% End
% Opto_attach
% Perc_pellet
% Perc_opto
% Code

% ingest this from .csv file
%dataTable='/Users/kim/Downloads/Combo Behavior Log - Slimmed down table.csv';
% BETTER TO USE TAB DELIMITER
% data_array=table2cell(readtable(dataTable,'Format','%s%s%s%s%s%s%s%s%s%s','Delimiter', '\t', 'HeaderLines', 0, 'ReadVariableNames', true));
data_array=table2cell(readtable(dataTable,'Format','%s%s%s%s%s%s%s%s%s%s','Delimiter', ',', 'HeaderLines', 0, 'ReadVariableNames', true));

%data_array=clearEmptyRows(data_array);

% save datetimes corresponding to dates in table
warning('off');
clear datetime_array
for i=1:size(data_array,1)
    temp=data_array{i,1};
    r=regexp(data_array{i,1},'/');
    if ~isempty(r)
        try
            datetime_array(i)=datetime(str2num(temp(r(end)+1:end)),str2num(temp(1:r(1)-1)),str2num(temp(r(1)+1:r(2)-1)));
        catch
%             disp(data_array{i,1});
        end
    else
        try
            datetime_array(i)=datetime(data_array{i,1},'InputFormat','uuuuMMdd');
        catch
%             disp(data_array{i,1});
        end
    end
end

% For each mouse, find all dates
% If mouse name doesn't match a mouse name in mouse_database.names, display
% Then fix nth_session in metadata for this mouse, based on which dates are
% in table
table_mousenames=unique(data_array(:,2));
for i=1:length(table_mousenames)
    currmousename=table_mousenames{i,2};
    [isInDB,correspondingMouseID]=findMouseNameInDatabase(currmousename,mouse_database);
    if ~isInDB
        disp(['This mouse name from table is not in database: ' currmousename]);
    else
        % Exclude days on training rig!!!
        datesForThisMouse=whichDatesForThisMouse(datetime_array,data_array,currmousename,excludeTrainingRig);
    end
    % Fix nth_session in metadata
    for j=1:length(datesForThisMouse)
        metadata.nth_session(eq(metadata.dateFromTextFile,datesForThisMouse(j)) & metadata.mouseid==correspondingMouseID)=j;
    end
end

if alsoFixOptoOnHere==true
    % if optoOnHere field does not match field Opto_attach, fix it
end

end

function datesForThisMouse=whichDatesForThisMouse(datetime_array,data_array,currmousename,excludeTrainingRig)

datesForThisMouse=[];
for i=1:size(data_array,1)
    if strcmp(data_array{i,1},currmousename)
        if ~ismember(datetime_array(i),datesForThisMouse)
            datesForThisMouse=[datesForThisMouse; datetime_array(i)]; 
        end
    end
end

end

function [isInDB,correspondingMouseID]=findMouseNameInDatabase(currmousename,mouse_database)

correspondingMouseID=nan;
isInDB=false;

% Structure of mouse_database needs to be:
% mouse_database.names{currmousenameind} is name of mouse
% mouse_database.mouse_id(currmousenameind) is mouseid

for i=1:length(mouse_database.names)
    if strcmp(mouse_database.names{i},currmousename)
        isInDB=true;
        correspondingMouseID=mouse_database.mouse_id(i);
    end
end

end

function data_array=clearEmptyRows(data_array)

% Dump any empty rows at the end of the table
for i=size(data_array,1):-1:1
    allempt=true;
    for j=1:size(data_array,2)
        if ~isempty(data_array{i,j})
            allempt=false;
            break
        end
    end
    if allempt==false
        break
    else
        % drop this empty row
        data_array=data_array(~ismember(1:size(data_array,1),i),:);
    end
end

end
