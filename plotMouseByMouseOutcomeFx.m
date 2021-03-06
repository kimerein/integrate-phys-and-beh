function plotMouseByMouseOutcomeFx(alltbt,metadata,trialTypes)

% select for these dprimes
dprimerange=[];

u=unique(metadata.mouseid);
j=0;
for i=1:length(u)
    metadata.sessid(metadata.mouseid==u(i))=metadata.nth_session(metadata.mouseid==u(i))+j;
    j=j+nanmax(metadata.sessid(metadata.mouseid==u(i)));
end

% for each mouse, filter alltbt
u=unique(metadata.mouseid);
touch_noLED_vectors=nan(length(u),2); % uncued then cued
touch_LED_vectors=nan(length(u),2);
noTouch_noLED_vectors=nan(length(u),2);
noTouch_LED_vectors=nan(length(u),2);
backup.alltbt=alltbt; 
backup.metadata=metadata;
backup.trialTypes=trialTypes;
temp=datestr(datetime('now'));
temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))='';
temp=temp(~isspace(temp));
saveDir=['\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\RT pairs data sets\' temp]; % where to save details of alltbt filtering and RT pairs data set
for i=1:length(u)
    alltbt=backup.alltbt;
    metadata=backup.metadata;
    trialTypes=backup.trialTypes;
    
    close all;
    
    currMouseID=u(i);
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
    [alltbt,trialTypes,metadata]=filtTbt(alltbt,trialTypes,tbt_filter.sortField,tbt_filter.range_values,metadata,tbt_filter.clock_progress);
    
    if ~isempty(dprimerange)
        [isreaching_out,dprimes]=getCuedResponseVsSuppression(alltbt,metadata,trialTypes,'cueZone_onVoff','all_reachBatch',trialTypes.led==0,1,0,1.5,1);
        close all;
        [metadata,alltbt,trialTypes]=add_dprimes_to_tbt(alltbt,trialTypes,metadata,dprimes);
        % filter settings
        tbt_filter.sortField='dprimes';
        tbt_filter.range_values=dprimerange;
        tbt_filter.name=[tbt_filter.sortField num2str(tbt_filter.range_values(1)) 'to' num2str(tbt_filter.range_values(2))];
        temp=tbt_filter.name;
        temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))='';
        temp=temp(~isspace(temp));
        tbt_filter.name=temp;
        tbt_filter.clock_progress=true;
        
        % filter alltbt
        [alltbt,trialTypes,metadata]=filtTbt(alltbt,trialTypes,tbt_filter.sortField,tbt_filter.range_values,metadata,tbt_filter.clock_progress);
    end
    
    % TOUCHED PELLET NO LED
    test.nInSequence=[2]; % defines trial pairs, e.g., 2 means will compare each trial with its subsequent trial, 3 means will compare each trial with the trial after next, etc.
    trial1='any(alltbt.all_reachBatch(:,94:end)>0.5,2) & trialTypes.touched_pellet==1  & trialTypes.isLongITI==1'; % e.g., mouse reached, touched pellet, no LED
    test.trial1=trial1;
    test.templateSequence2_cond=eval(trial1);
%     trial2='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1';
    trial2='trialTypes.led==0';
    test.trial2=trial2;
    test.templateSequence2_end=eval(trial2);
    test.fillInBetweenWithAnything=false; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
    test.event_name=['alltrials' tbt_filter.name 'inBetweenAnything' num2str(test.fillInBetweenWithAnything)];
    saveDir2=[saveDir '\' test.event_name];
    mkdir(saveDir2);
    % save settings for paired RT data set
    save([saveDir2 '\test_settings.mat'],'test');
    % build paired RT data set
    fakeCueInd=50; % in indices
    skipCorrected=true;
    [dataset,correctedDistributions]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); % actually, this function just builds the RT pairs dataset
    rrs_touchpellet_noLED=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',false);

    % TOUCHED PELLET LED
    test.nInSequence=[2]; % defines trial pairs, e.g., 2 means will compare each trial with its subsequent trial, 3 means will compare each trial with the trial after next, etc.
    trial1=['any(alltbt.all_reachBatch(:,94:end)>0.5,2) & trialTypes.touched_pellet==1 & trialTypes.isLongITI==1']; % e.g., mouse reached, touched pellet, no LED
    test.trial1=trial1;
    test.templateSequence2_cond=eval(trial1);
%     trial2='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1';
    trial2='trialTypes.led==0';
    test.trial2=trial2;
    test.templateSequence2_end=eval(trial2);
    test.fillInBetweenWithAnything=false; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
    test.event_name=['alltrials' tbt_filter.name 'inBetweenAnything' num2str(test.fillInBetweenWithAnything)];
    saveDir2=[saveDir '\' test.event_name];
    mkdir(saveDir2);
    % save settings for paired RT data set
    save([saveDir2 '\test_settings.mat'],'test');
    % build paired RT data set
    fakeCueInd=50; % in indices
    skipCorrected=true;
    [dataset,correctedDistributions]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); % actually, this function just builds the RT pairs dataset
    rrs_touchpellet_LED=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',false);

    % MISSED PELLET NO LED
    test.nInSequence=[2]; % defines trial pairs, e.g., 2 means will compare each trial with its subsequent trial, 3 means will compare each trial with the trial after next, etc.
    trial1='any(alltbt.all_reachBatch(:,94:end)>0.5,2) & trialTypes.touched_pellet==0 & trialTypes.consumed_pellet==0'; 
    test.trial1=trial1;
    test.templateSequence2_cond=eval(trial1);
%     trial2='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1';
    trial2='trialTypes.led==0';
    test.trial2=trial2;
    test.templateSequence2_end=eval(trial2);
    test.fillInBetweenWithAnything=false; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
    test.event_name=['alltrials' tbt_filter.name 'inBetweenAnything' num2str(test.fillInBetweenWithAnything)];
    saveDir2=[saveDir '\' test.event_name];
    mkdir(saveDir2);
    % save settings for paired RT data set
    save([saveDir2 '\test_settings.mat'],'test');
    % build paired RT data set
    fakeCueInd=50; % in indices
    skipCorrected=true;
    [dataset,correctedDistributions]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); % actually, this function just builds the RT pairs dataset
    rrs_misspellet_noLED=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',false);

    % MISSED PELLET LED
    test.nInSequence=[2]; % defines trial pairs, e.g., 2 means will compare each trial with its subsequent trial, 3 means will compare each trial with the trial after next, etc.
    trial1=['any(alltbt.all_reachBatch(:,94:end)>0.5,2) & trialTypes.touched_pellet==0 & trialTypes.consumed_pellet==0']; 
    test.trial1=trial1;
    test.templateSequence2_cond=eval(trial1);
%     trial2='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1';
    trial2='trialTypes.led==0';
    test.trial2=trial2;
    test.templateSequence2_end=eval(trial2);
    test.fillInBetweenWithAnything=false; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
    test.event_name=['alltrials' tbt_filter.name 'inBetweenAnything' num2str(test.fillInBetweenWithAnything)];
    saveDir2=[saveDir '\' test.event_name];
    mkdir(saveDir2);
    % save settings for paired RT data set
    save([saveDir2 '\test_settings.mat'],'test');
    % build paired RT data set
    fakeCueInd=50; % in indices
    skipCorrected=true;
    [dataset,correctedDistributions]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); % actually, this function just builds the RT pairs dataset
    rrs_misspellet_LED=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',false);

    % ALL OUTCOMES
    test.nInSequence=[2]; % defines trial pairs, e.g., 2 means will compare each trial with its subsequent trial, 3 means will compare each trial with the trial after next, etc.
    trial1='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1'; 
%     trial1='any(alltbt.all_reachBatch(:,94:end)>0.5,2) & trialTypes.led==0';
%     trial1='trialTypes.led==0';
%     trial1='trialTypes.led==0';
    test.trial1=trial1;
    test.templateSequence2_cond=eval(trial1);
    trial2='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1';
%     trial2='trialTypes.led==0';
    test.trial2=trial2;
    test.templateSequence2_end=eval(trial2);
    test.fillInBetweenWithAnything=false; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
    test.event_name=['alltrials' tbt_filter.name 'inBetweenAnything' num2str(test.fillInBetweenWithAnything)];
    saveDir2=[saveDir '\' test.event_name];
    mkdir(saveDir2);
    % save settings for paired RT data set
    save([saveDir2 '\test_settings.mat'],'test');
    % build paired RT data set
    fakeCueInd=50; % in indices
    skipCorrected=true;
    [dataset,correctedDistributions]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); % actually, this function just builds the RT pairs dataset
    rrs_alloutcomes=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',false);

    [uncued,cued]=plotShiftForEachSectionOfSession(rrs_touchpellet_noLED,rrs_alloutcomes,false); % data 2 is ref
    touch_noLED_vectors(i,:)=[uncued cued]; % uncued then cued
    [uncued,cued]=plotShiftForEachSectionOfSession(rrs_touchpellet_LED,rrs_alloutcomes,false); % data 2 is ref
    touch_LED_vectors(i,:)=[uncued cued]; % uncued then cued
    [uncued,cued]=plotShiftForEachSectionOfSession(rrs_misspellet_noLED,rrs_alloutcomes,false); % data 2 is ref
    noTouch_noLED_vectors(i,:)=[uncued cued]; % uncued then cued
    [uncued,cued]=plotShiftForEachSectionOfSession(rrs_misspellet_LED,rrs_alloutcomes,false); % data 2 is ref
    noTouch_LED_vectors(i,:)=[uncued cued]; % uncued then cued
end

figure();
for i=1:size(touch_noLED_vectors,1)
    if any(isnan([touch_noLED_vectors(i,1) touch_noLED_vectors(i,2)]))
    else
        quiver(0,0,touch_noLED_vectors(i,1),touch_noLED_vectors(i,2),'Color','k');
    end
    hold on;
    if any(isnan([touch_LED_vectors(i,1) touch_LED_vectors(i,2)]))
    else
        quiver(0,0,touch_LED_vectors(i,1),touch_LED_vectors(i,2),'Color','r');
    end
end
title('After touch pellet');

figure();
for i=1:size(noTouch_noLED_vectors,1)
    if any(isnan([noTouch_noLED_vectors(i,1) noTouch_noLED_vectors(i,2)]))
    else
        quiver(0,0,noTouch_noLED_vectors(i,1),noTouch_noLED_vectors(i,2),'Color','k');
    end
    hold on;
    if any(isnan([noTouch_LED_vectors(i,1) noTouch_LED_vectors(i,2)]))
    else
        quiver(0,0,noTouch_LED_vectors(i,1),noTouch_LED_vectors(i,2),'Color','r');
    end
end
title('After miss pellet');

end