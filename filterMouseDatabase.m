function db_out=filterMouseDatabase(db,firstField,secondField,values,excludeControls)

% set up reference pointers
indIntoDbsByVids=1;
for i=1:length(db.db_bymouse)
    db.db_bymouse(i).indInDbsByVids=indIntoDbsByVids:indIntoDbsByVids+length(db.db_bymouse(i).low_speed_videos)-1;
    db.dbs.musbymus_indsInVids(indIntoDbsByVids:indIntoDbsByVids+length(db.db_bymouse(i).low_speed_videos)-1)=1:length(db.db_bymouse(i).low_speed_videos);
    indIntoDbsByVids=indIntoDbsByVids+length(db.db_bymouse(i).low_speed_videos);
end
% reference fields for behLog are
% behLog.indsInto_vids
% dbs.indsInto_behLog

% filter by which (mouse by mouse or video or behavior log database)
filtfirst=db.(firstField);
k=1;
switch firstField
    case 'db_bymouse' 
        useMouse=nan(1,length(filtfirst));
        takeIndsInDbs=zeros(1,length(db.dbs.vids_to_match_mouseIDs));
        for i=1:length(filtfirst)
            % will either include all data from this mouse or throw out all
            % data from this mouse
            if any(ismember(filtfirst(i).(secondField),values))
                % include this mouse
                useMouse(i)=1;
                db_out.db_bymouse(k)=filtfirst(i);
                takeIndsInDbs(db.db_bymouse(i).indInDbsByVids)=1;
            else
                % exclude this mouse
                useMouse(i)=0;
            end
        end
        if isfield(dbs,'indsInto_behLog')
            takeRowsInBehLog=dbs.indsInto_behLog(takeIndsInDbs==1);
        end
        % filt dbs by mice
        f=fieldnames(db.dbs);
        for j=1:length(f)
            temp=db.dbs.(f{j});
            db_out.dbs.(f{j})=temp(takeIndsInDbs==1);
        end       
        % filt behLog by mice
        if isfield(db,'behLog')
            db_out.behLog=db.behLog(takeRowsInBehLog,:);
        end  
    case 'dbs'
        takeIndsInDbs=ismember(filtfirst.(secondField),values);
        f=fieldnames(db.dbs);
        for j=1:length(f)
            temp=db.dbs(f{j});
            db_out.dbs.(f{j})=temp(takeIndsInDbs==1);
        end
        % filt mouse by mouse
        for j=1:length(db.db_bymouse)
            takeIndsInMouse=db.dbs.musbymus_indsInVids(takeIndsInDbs==1 && ismember(db.dbs.mouseIDs_to_match_vids,db.db_bymouse(j).mouse_id));
            expectedLength=length(db.db_bymouse(j).sessids);
            f=fieldnames(db.db_bymouse(j));
            for l=1:length(f)
                temp=db.db_bymouse(j).(f{l});
                if length(temp)==expectedLength
                    db_out.db_bymouse(j).(f{l})=temp(takeIndsInMouse);
                else
                    db_out.db_bymouse(j).(f{l})=temp;
                end
            end
        end
        % filt behLog
        if isfield(dbs,'indsInto_behLog')
            takeRowsInBehLog=dbs.indsInto_behLog(takeIndsInDbs==1);
        end
        if isfield(db,'behLog')
            db_out.behLog=db.begLog(takeRowsInBehLog,:);
        end
    case 'behLog'
        
    otherwise
        error('firstField is not a field of mouse database');
end
