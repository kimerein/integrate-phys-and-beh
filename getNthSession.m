function metadata=getNthSession(dataTable,mouse_database,metadata,alsoFixOptoOnHere,excludeTrainingRig)

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
% Download google spreadsheet as .tsv, then change extension to .csv
data_array=table2cell(readtable(dataTable,'Format','%s%s%s%s%s%s%s%s%s%s','Delimiter', '\t', 'HeaderLines', 0, 'ReadVariableNames', true));
% data_array=table2cell(readtable(dataTable,'Format','%s%s%s%s%s%s%s%s%s%s','Delimiter', ',', 'HeaderLines', 0, 'ReadVariableNames', true));

%data_array=clearEmptyRows(data_array);

if ~isfield(metadata,'pelletPresentFromTable')
    metadata.pelletPresentFromTable=nan(size(metadata.nth_session));
end

% save datetimes corresponding to dates in table
warning('off');
clear datetime_array
for i=1:size(data_array,1)
    temp=data_array{i,1};
    r=regexp(data_array{i,1},'/');
    if ~isempty(r)
        try
            year=temp(r(end)+1:end);
            if length(year)==2
                year=['20' year];
            end
            year=str2num(year);
            datetime_array(i)=datetime(year,str2num(temp(1:r(1)-1)),str2num(temp(r(1)+1:r(2)-1)));
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
    currmousename=table_mousenames{i};
    [isInDB,correspondingMouseID]=findMouseNameInDatabase(currmousename,mouse_database);
    if ~isInDB
        disp(['This mouse name from table is not in database: ' currmousename]);
        continue
    else
        % Exclude days on training rig!!!
        [datesForThisMouse,optoAttachForThisMouse,percOptoForThisMouse,percPelletForThisMouse]=whichDatesForThisMouse(datetime_array,data_array,currmousename,excludeTrainingRig);
    end
    % Fix nth_session in metadata
    for j=1:length(datesForThisMouse)
        metadata.nth_session(eq(metadata.dateFromTextFile,datesForThisMouse(j)) & metadata.mouseid==correspondingMouseID)=j;
        if alsoFixOptoOnHere==true
            % if optoOnHere field does not match field Opto_attach, fix it
            switch optoAttachForThisMouse{j}
                case 'yes'
                    if str2num(percOptoForThisMouse{j})~=0
                        metadata.optoOnHere(eq(metadata.dateFromTextFile,datesForThisMouse(j)) & metadata.mouseid==correspondingMouseID)=1;
                    else
                        metadata.optoOnHere(eq(metadata.dateFromTextFile,datesForThisMouse(j)) & metadata.mouseid==correspondingMouseID)=0;
                    end
                case 'no'
                    metadata.optoOnHere(eq(metadata.dateFromTextFile,datesForThisMouse(j)) & metadata.mouseid==correspondingMouseID)=0;
                otherwise
                    metadata.optoOnHere(eq(metadata.dateFromTextFile,datesForThisMouse(j)) & metadata.mouseid==correspondingMouseID)=0;
            end
        end
        try
            metadata.pelletPresentFromTable(eq(metadata.dateFromTextFile,datesForThisMouse(j)) & metadata.mouseid==correspondingMouseID)=str2double(percPelletForThisMouse{j});
        catch
            disp('could not convert str2num in getNthSession.m');
        end
    end
end

end

function [datesForThisMouse,optoAttachForThisMouse,percOptoForThisMouse,percPelletForThisMouse]=whichDatesForThisMouse(datetime_array,data_array,currmousename,excludeTrainingRig)

datesForThisMouse(1)=datetime(1900,3,3); % just a placeholder to define
optoOnForThisMouse={}; percOptoForThisMouse={}; percPelletForThisMouse={};
k=1;

for i=1:size(data_array,1)
    if strcmp(data_array{i,2},currmousename)
        if ~ismember(datetime_array(i),datesForThisMouse)
            if excludeTrainingRig==true
                if ~strcmp(data_array{i,3},'Training')
                    datesForThisMouse=[datesForThisMouse; datetime_array(i)];
                    optoAttachForThisMouse{k}=data_array{i,7}; percOptoForThisMouse{k}=data_array{i,9}; percPelletForThisMouse{k}=data_array{i,8}; k=k+1;
                end
            else
                datesForThisMouse=[datesForThisMouse; datetime_array(i)]; 
                optoAttachForThisMouse{k}=data_array{i,7}; percOptoForThisMouse{k}=data_array{i,9}; percPelletForThisMouse{k}=data_array{i,8}; k=k+1;
            end
        end
    end
end
% drop placeholder
datesForThisMouse=datesForThisMouse(2:end);

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
