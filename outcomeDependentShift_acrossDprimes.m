function outcomeDependentShift_acrossDprimes(backup_alltbt,backup_trialTypes,backup_metadata)

dprimesBins={[-100 0]; [0 0.5]; [0.5 1.8]; [1.8 3]}; % success
f1=[];
f2=[];
for i=1:length(dprimesBins)
    alltbt=backup_alltbt;
    trialTypes=backup_trialTypes;
    metadata=backup_metadata;
    saveDir=[]; 
    % filter settings
    tbt_filter.sortField='dprimes';
    tbt_filter.range_values=dprimesBins{i};
    tbt_filter.name=[tbt_filter.sortField num2str(tbt_filter.range_values(1)) 'to' num2str(tbt_filter.range_values(2))];
    temp=tbt_filter.name;
    temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))='';
    temp=temp(~isspace(temp));
    tbt_filter.name=temp;
    tbt_filter.clock_progress=true;
    % filter alltbt
    [alltbt,trialTypes,metadata]=filtTbt(alltbt,trialTypes,tbt_filter.sortField,tbt_filter.range_values,metadata,tbt_filter.clock_progress);
    [f1,f2]=outcomeDependentShift_wrapper(alltbt,trialTypes,metadata,saveDir,f1,f2);  
end

