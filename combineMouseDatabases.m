function db_out=combineMouseDatabases(db1,db2)

% Make sure that mouse IDs match for db1 and db2
T_mouseToID=extractMouseNametoIDmapping(db1.db_bymouse);
[db2.db_bymouse,db2.dbs]=assignIDnums(db2.db_bymouse,T_mouseToID);

% Concatenate behLogs
if isfield(db1,'behLog') && isfield(db2,'behLog')
    % Look for duplicate days in behLog
    % A duplicate is the same day, same mouse
    % Take the first instance of any duplicate
    % Assume no duplication within each beginning behLog
    toRemove=false(1,length(db2.behLog.Date));
    for i=1:length(db1.behLog.Date)
        for j=1:length(db2.behLog.Date)
            if db1.behLog.Date(i)==db2.behLog.Date(j) % first col of behLog is datetime
                % check if this is the same mouse
                if strcmp(db1.behLog.Mouse(i),db2.behLog.Mouse(j))
                    % is a duplicate
                    % remove from db2.behLog
                    toRemove(j)=true;
                end
            end
        end
    end
    db2.behLog=db2.behLog(~toRemove,:);
    % Concatenate behLogs
    f1=fieldnames(db1.behLog);
    f2=fieldnames(db2.behLog);
    f=unique(cat(1,f1,f2));
    for i=1:length(f)
        if isfield(db1.behLog,f{i}) && isfield(db2.behLog,f{i}) % in both behLogs, duplicates have already been removed
            db_out.behLog.(f{i})=cat(1,db1.behLog.(f{i}),db2.behLog.(f{i}));
        elseif isfield(db1.behLog,f{i}) && ~isfield(db2.behLog,f{i}) % only in db1
            if iscell(db1.behLog.(f{i}))
                db_out.behLog.(f{i})=cat(1,db1.behLog.(f{i}),cell(length(db2.behLog.Date),1));
            else
                db_out.behLog.(f{i})=cat(1,db1.behLog.(f{i}),nan(length(db2.behLog.Date),1));
            end
        elseif ~isfield(db1.behLog,f{i}) && isfield(db2.behLog,f{i}) % only in db2
            if iscell(db1.behLog.(f{i}))
                db_out.behLog.(f{i})=cat(1,cell(length(db1.behLog.Date),1),db2.behLog.(f{i}));
            else
                db_out.behLog.(f{i})=cat(1,nan(length(db1.behLog.Date),1),db2.behLog.(f{i}));
            end
        end
    end
end

% Concatenate mouse by mouse
% Duplicates are OK here
added_mice_in_db2=false(1,length(db2.db_bymouse));
added_mice_in_db1=false(1,length(db1.db_bymouse));
l=1;
for i=1:length(db1.db_bymouse)
    for j=1:length(db2.db_bymouse)
        if strcmp(db1.db_bymouse(i).mouse_name,db2.db_bymouse(j).mouse_name)
            % same mouse
            added_mice_in_db2(j)=true;
            added_mice_in_db1(i)=true;
            % check for same mouse_id
            if db1.db_bymouse(i).mouse_id==db2.db_bymouse(j).mouse_id
                % OK, continue
            else
                % problem
                error('Mouse IDs do not match in db1 and db2');
            end
            % concatenate other fields
            f=fieldnames(db1.db_bymouse(i)); % assume the same fields in db1.db_bymouse and db2.db_bymouse
            for k=1:length(f)
                if length(db1.db_bymouse(i).(f{i}))>1
                    db_out.db_bymouse(l).(f{i})=cat(2,db1.db_bymouse(i).(f{i}),db2.db_bymouse(j).(f{i}));
                else
                    db_out.db_bymouse(l).(f{i})=db1.db_bymouse(i).(f{i});
                end
            end
            l=l+1; % increment for next mouse to add
        end
    end
end
% add any mice that were only in db1 or db2
fi=find(added_mice_in_db1==false);
for i=1:length(fi)
    % mice that were in db1 but not db2
    f=fieldnames(db1.db_bymouse(i));
    for k=1:length(f)
        db_out.db_bymouse(l).(f{i})=db1.db_bymouse(i).(f{i});
    end
    l=l+1; % increment for next mouse to add   
end
fi=find(added_mice_in_db2==false);
for i=1:length(fi)
    % mice that were in db2 but not db1
    f=fieldnames(db2.db_bymouse(i));
    for k=1:length(f)
        db_out.db_bymouse(l).(f{i})=db2.db_bymouse(i).(f{i});
    end
    l=l+1; % increment for next mouse to add   
end

% Re-run assign mouse IDs and sess IDs to be sure that session IDs are
% correctly ordered in output database
% This will also populate db_out.dbs
T_mouseToID=extractMouseNametoIDmapping(db_out.db_bymouse);
[db_out.db_bymouse,db_out.dbs]=assignIDnums(db_out.db_bymouse,T_mouseToID);
% Match behLog and dbs by vids
[db_out.behLog,db_out.dbs]=matchBehLogAndDbs_by_vids(db_out.behLog,db_out.dbs);

end

function [behLog,dbs]=matchBehLogAndDbs_by_vids(behLog,dbs)

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

function T=extractMouseNametoIDmapping(db_bymouse)

mapping=cell(length(db_bymouse),2); % first col is mouse name, second col is mouse ID
for i=1:length(db_bymouse)
    mapping{i,1}=db_bymouse(i).mouse_name;
    mapping{i,2}=db_bymouse(i).mouse_id;
end

T=cell2table(mapping,'VariableNames',{'name','id'});

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
