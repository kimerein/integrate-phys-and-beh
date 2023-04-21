function [alltbt_backup,trialTypes_backup,metadata_backup]=get_DistractorDprime_per_mouse(alltbt,trialTypes,metadata)

metadata.distract_dprimes=nan(size(metadata.dprimes));
trialTypes.distract_dprimes=nan(size(trialTypes.dprimes));

% Realign alltbt to distractor instead of cue
alltbt=realignToADistractor(alltbt,'movie_distractor');
% Change cueZone_onVoff field to be movie_distractor
alltbt.cueZone_onVoff=alltbt.movie_distractor;

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
u=unique(metadata.mouseid);
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
    [metadata,alltbt,trialTypes]=add_dprimes_to_tbt(alltbt,trialTypes,metadata,dprimes);
    
    % add back to multi-mouse alltbt
    if ~isfield(alltbt_backup,'distract_dprimes')
        alltbt_backup.distract_dprimes=nan(size(alltbt.cue,1),1);
    end
    f={'distract_dprimes'}; % just add dprimes, other fields OK
    for j=1:length(f)
        temp=alltbt_backup.(f{j});
        if length(size(temp))>1
            temp(tookThese,:)=alltbt.(f{j});
        else
            temp(tookThese)=alltbt.(f{j});
        end
        alltbt_backup.(f{j})=temp;
    end
    f={'distract_dprimes'}; % just add dprimes, other fields OK
    for j=1:length(f)
        temp=metadata_backup.(f{j});
        temp(tookThese)=metadata.(f{j});
        metadata_backup.(f{j})=temp;
    end
    f={'distract_dprimes'}; % just add dprimes, other fields OK
    for j=1:length(f)
        temp=trialTypes_backup.(f{j});
        temp(tookThese)=trialTypes.(f{j});
        trialTypes_backup.(f{j})=temp;
    end
end

end

function alltbt=realignToADistractor(alltbt,distractName)

randomDistractor=false;

% for each trial, randomly choose one of the distractors and align trial to
% this instead of cue
% where is cue currently
[~,ma]=nanmax(nanmean(alltbt.cueZone_onVoff,1),[],2);
temp=alltbt.(distractName);
temp(temp>=0.5)=1; temp(temp<0.5)=0;
alltbt.(distractName)=temp;
figure(); plot(nanmean(alltbt.(distractName),1),'Color','k');
shiftBy=nan(size(temp,1),1);
if randomDistractor==true
    for i=1:size(temp,1)
        % find distractor onsets
        f=find(diff(temp(i,:))==1);
        if isempty(f)
            continue
        end
        f=f(randperm(length(f)));
        f=f(1);
        % realign to this
        shiftBy(i)=f-ma; % if positive, will shift backwards, else will shift forward in time
    end
else
    % The distractor after the cue
    for i=1:size(temp,1)
        % find distractor onsets
        tempie=temp(i,:);
        tempie(1:ma)=0;
        f=find(diff(tempie)==1);
        if isempty(f)
            continue
        end
        f=f(randperm(length(f)));
        f=f(1);
        % realign to this
        shiftBy(i)=f-ma; % if positive, will shift backwards, else will shift forward in time
    end
end
% shift all relevant fields
f={'all_reachBatch','cue','isChewing','isHold','movie_distractor','optoOn','optoZone','pawOnWheel','pelletPresent','pelletmissingreach_reachStarts',...
   'reachBatch_all_pawOnWheel','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','reachBatch_success_reachStarts','reachStarts',...
   'reachStarts_pelletPresent'};
for i=1:length(f)
    temp=alltbt.(f{i});
    for j=1:size(temp,1)
        if ~isnan(shiftBy(j))
            if shiftBy(j)>0
                temp(j,:)=circshift(temp(j,:),[0 shiftBy(j)]);
            end
        end
    end
    alltbt.(f{i})=temp;
end

hold on; plot(nanmean(alltbt.(distractName),1),'Color','b');
title('black before, blue after shift');

end

function [metadata,alltbt,out]=add_dprimes_to_tbt(varargin)

if length(varargin)==4
    alltbt=varargin{1};
    out=varargin{2};
    metadata=varargin{3};
    dprimes=varargin{4};
    disp('input 1: alltbt, input 2: trialTypes, input 3: metadata, input 4: dprimes');
elseif length(varargin)==7
    alltbt=varargin{1};
    out=varargin{2};
    metadata=varargin{3};
    dprimes=varargin{4};
    reachName=varargin{5};
    cueName=varargin{6};
    settings=varargin{7};
end

if ~isempty(dprimes)
    disp('adding dprimes passed in');
    % assume dprimes are ordered by session
    u=unique(metadata.sessid);
    for i=1:length(u)
        alltbt.distract_dprimes(metadata.sessid==u(i))=dprimes(i);
    end
    metadata.distract_dprimes=alltbt.distract_dprimes;
    out.distract_dprimes=alltbt.distract_dprimes;
else
    [dprimes]=get_dprime_per_session(alltbt,out,metadata,reachName,cueName,settings);
    
    metadata.distract_dprimes=nan(size(metadata.sessid));
    u=unique(metadata.sessid);
    for i=1:length(u)
        curru=u(i);
        metadata.distract_dprimes(ismember(metadata.sessid,curru))=dprimes(i);
    end
    
    alltbt.distract_dprimes=metadata.distract_dprimes;
    out.distract_dprimes=metadata.distract_dprimes;
end

end