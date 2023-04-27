function [alltbt_backup,trialTypes_backup,metadata_backup,isreachout_permouse,u]=get_dprime_per_mouse(varargin)

if length(varargin)==4
    settingsForDp=[];
    alltbt=varargin{1};
    trialTypes=varargin{2};
    metadata=varargin{3};
    getRatesInstead=varargin{4};
else
    alltbt=varargin{1};
    trialTypes=varargin{2};
    metadata=varargin{3};
    getRatesInstead=varargin{4};
    settingsForDp=varargin{5};
end

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
    [isreaching_out,dprimes]=getCuedResponseVsSuppression(alltbt,metadata,trialTypes,'cueZone_onVoff','all_reachBatch',[],1,settingsForDp.reachAfterCueWindow_start,settingsForDp.reachAfterCueWindow_end,false,getRatesInstead,settingsForDp);
    if isfield(isreaching_out,'cued_reach_rate')
        RRcued=isreaching_out.cued_reach_rate;
        RRuncued=isreaching_out.noncued_reach_rate;
    else
        RRcued=[];
        RRuncued=[];
    end
    metadata=addRR(metadata,RRcued,'reachrate_cued'); 
    metadata=addRR(metadata,RRuncued,'reachrate_uncued');
    isreachout_permouse{i}=isreaching_out;

    if ~isempty(dprimes)
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

    if ~isempty(RRcued)
        % add to multi-mouse metadata
        if ~isfield(metadata_backup,'reachrate_cued')
            metadata_backup.reachrate_cued=nan(size(metadata_backup.sessid));
        end
        f={'reachrate_cued'}; 
        for j=1:length(f)
            temp=metadata_backup.(f{j});
            temp(tookThese)=metadata.(f{j});
            metadata_backup.(f{j})=temp;
        end
    end
    if ~isempty(RRuncued)
        % add to multi-mouse metadata
        if ~isfield(metadata_backup,'reachrate_uncued')
            metadata_backup.reachrate_uncued=nan(size(metadata_backup.sessid));
        end
        f={'reachrate_uncued'}; 
        for j=1:length(f)
            temp=metadata_backup.(f{j});
            temp(tookThese)=metadata.(f{j});
            metadata_backup.(f{j})=temp;
        end
    end
end

end

function metadata=addRR(metadata,rr,fname)

if isempty(rr)
    return
end

% assume ordered by sess
u=unique(metadata.sessid);
if ~isfield(metadata,fname)
    metadata.(fname)=nan(size(metadata.sessid));
end
for i=1:length(u)
    temp=metadata.(fname);
    temp(metadata.sessid==u(i))=rr(i);
    metadata.(fname)=temp;
end

end