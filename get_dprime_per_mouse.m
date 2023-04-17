function [alltbt_backup,trialTypes_backup,metadata_backup,isreachout_permouse,u]=get_dprime_per_mouse(alltbt,trialTypes,metadata)

% Assign unique sessids
u=unique(metadata.mouseid);
j=0;
for i=1:length(u)
    metadata.sessid(metadata.mouseid==u(i))=metadata.nth_session(metadata.mouseid==u(i))+j;
    j=j+nanmax(metadata.sessid(metadata.mouseid==u(i)));
end

alltbt_backup=alltbt;
trialTypes_backup=trialTypes;
metadata_backup=metadata;

% Per mouse
u=unique(metadata.mouseid); % u is the mouseid
isreachout_permouse=cell(length(u),1);
for i=1:length(u)
    currMouseID=u(i);
    alltbt=alltbt_backup; 
    trialTypes=trialTypes_backup; 
    metadata=metadata_backup;
    
    temp=datestr(datetime('now'));
    temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))='';
    temp=temp(~isspace(temp));
    saveDir=['\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\RT pairs data sets\' temp]; % where to save details of alltbt filtering and RT pairs data set
    
    % filter settings
    tbt_filter.sortField='mouseid';
    tbt_filter.range_values=[currMouseID-0.5 currMouseID+0.5];
    tbt_filter.name=[tbt_filter.sortField num2str(tbt_filter.range_values(1)) 'to' num2str(tbt_filter.range_values(2))];
    temp=tbt_filter.name;
    temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))='';
    temp=temp(~isspace(temp));
    tbt_filter.name=temp;
    tbt_filter.clock_progress=true;
    
    % filter alltbt
    [alltbt,trialTypes,metadata,tookThese]=filtTbt(alltbt,trialTypes,tbt_filter.sortField,tbt_filter.range_values,metadata,tbt_filter.clock_progress);
    
    % get dprimes for this mouse
    settingsForDp=settingsForDprimes(alltbt,'cueZone_onVoff',false);
    [isreaching_out,dprimes]=getCuedResponseVsSuppression(alltbt,metadata,trialTypes,'cueZone_onVoff','all_reachBatch',[],1,settingsForDp.reachAfterCueWindow_start,settingsForDp.reachAfterCueWindow_end,false,0);
    isreachout_permouse{i}=isreaching_out;
    [metadata,alltbt,trialTypes]=add_dprimes_to_tbt(alltbt,trialTypes,metadata,dprimes);
    
    % add back to multi-mouse alltbt
    if ~isfield(alltbt_backup,'dprimes')
        alltbt_backup.dprimes=nan(size(alltbt.cue,1),1);
    end
    f={'dprimes'}; % just add dprimes, other fields OK
    for j=1:length(f)
        temp=alltbt_backup.(f{j});
        if length(size(temp))>1
            temp(tookThese,:)=alltbt.(f{j});
        else
            temp(tookThese)=alltbt.(f{j});
        end
        alltbt_backup.(f{j})=temp;
    end
    f={'dprimes'}; % just add dprimes, other fields OK
    for j=1:length(f)
        temp=metadata_backup.(f{j});
        temp(tookThese)=metadata.(f{j});
        metadata_backup.(f{j})=temp;
    end
    f={'dprimes'}; % just add dprimes, other fields OK
    for j=1:length(f)
        temp=trialTypes_backup.(f{j});
        temp(tookThese)=trialTypes.(f{j});
        trialTypes_backup.(f{j})=temp;
    end
end