function get_dprime_per_mouse()
backup.alltbt=alltbt;
backup.trialTypes=trialTypes;
backup.metadata=metadata;

%% Optional: Fix sessids to match nth_sessions
alltbt=backup.alltbt; trialTypes=backup.trialTypes; metadata=backup.metadata;
u=unique(metadata.mouseid);
j=0;
for i=1:length(u)
    metadata.sessid(metadata.mouseid==u(i))=metadata.nth_session(metadata.mouseid==u(i))+j;
    j=j+nanmax(metadata.sessid(metadata.mouseid==u(i)));
end

%% perform any filtering on alltbt
% for example, filter by d-prime

temp=datestr(datetime('now'));
temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))=''; 
temp=temp(~isspace(temp));
saveDir=['\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\RT pairs data sets\' temp]; % where to save details of alltbt filtering and RT pairs data set

% filter settings
tbt_filter.sortField='mouseid';
tbt_filter.range_values=[-100 100];
tbt_filter.name=[tbt_filter.sortField num2str(tbt_filter.range_values(1)) 'to' num2str(tbt_filter.range_values(2))];
temp=tbt_filter.name;
temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))=''; 
temp=temp(~isspace(temp));
tbt_filter.name=temp;
tbt_filter.clock_progress=true;

% filter alltbt
[alltbt,trialTypes,metadata]=filtTbt(alltbt,trialTypes,tbt_filter.sortField,tbt_filter.range_values,metadata,tbt_filter.clock_progress);

[isreaching_out,dprimes]=getCuedResponseVsSuppression(alltbt,metadata,trialTypes,'cueZone_onVoff','all_reachBatch',[],1,0,1.5,1);


end