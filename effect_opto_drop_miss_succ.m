% script for no effect of str silencing on drop, miss, success

trialTypes.optoGroup_1back=[1; trialTypes.optoGroup(1:end-1)];
trialTypes.optoGroup_1forward=[trialTypes.optoGroup(2:end); 0];
trialTypes.optoGroup_2forward=[trialTypes.optoGroup(3:end); zeros(2,1)];
trialTypes.optoGroup_2back=[ones(2,1); trialTypes.optoGroup(1:end-2)];

% trial2='trialTypes.led==0 & trialTypes.optoGroup_1forward==2';
% trial2='trialTypes.led==0 & trialTypes.led_1forward==1';
% trial2='trialTypes.led==1';
% trial2='trialTypes.optoGroup==2';
% trial2='trialTypes.optoGroup==2 & (trialTypes.led_1forward==0 | trialTypes.led_1back==0 | trialTypes.led_2back==0 | trialTypes.led_2forward==0)';
trial2='trialTypes.led==1 & (trialTypes.led_1forward==0 | trialTypes.led_1back==0 | trialTypes.led_2back==0 | trialTypes.led_2forward==0)';
% trial2='trialTypes.led==0 & (trialTypes.optoGroup_1forward==2 | trialTypes.optoGroup_1back==2 | trialTypes.optoGroup_2back==2 | trialTypes.optoGroup_2forward==2)';
% trial2='trialTypes.led==0 & (trialTypes.led_1forward==1 | trialTypes.led_1back==1 | trialTypes.led_2back==1 | trialTypes.led_2forward==1)';

backup.alltbt=alltbt;
backup.trialTypes=trialTypes;
backup.metadata=metadata;

% Filter by dprime
temp=datestr(datetime('now'));
temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))=''; 
temp=temp(~isspace(temp));
saveDir=['/Volumes/Neurobio/MICROSCOPE/Kim/RT pairs data sets/' temp]; % where to save details of alltbt filtering and RT pairs data set
% BEGINNER
alltbt.mouseid=metadata.mouseid;
alltbt.sessid=metadata.sessid;
trialTypes.sessid=metadata.sessid;
tbt_filter.sortField='dprimes';
tbt_filter.range_values=[-100 0.25]; % beginner: d<0.25, intermediate: 0.25<=d<0.75, expert: d>=0.75
tbt_filter.name=[tbt_filter.sortField num2str(tbt_filter.range_values(1)) 'to' num2str(tbt_filter.range_values(2))];
temp=tbt_filter.name;
temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))=''; 
temp=temp(~isspace(temp));
tbt_filter.name=temp;
tbt_filter.clock_progress=true;
% filter alltbt
[alltbt,trialTypes,metadata]=filtTbt(alltbt,trialTypes,tbt_filter.sortField,tbt_filter.range_values,metadata,tbt_filter.clock_progress);

% SUCCESS
alltbt.all_reachBatch=alltbt.reachBatch_success_reachStarts;
test.nInSequence=[2]; % defines trial pairs, e.g., 2 means will compare each trial with its subsequent trial, 3 means will compare each trial with the trial after next, etc.
trial1='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1';
test.templateSequence2_cond=eval(trial1);
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=true; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name=['alltrials' tbt_filter.name 'inBetweenAnything' num2str(test.fillInBetweenWithAnything)];
test.onlyConsiderReachType=[];
if ~isempty(test.onlyConsiderReachType)
    alltbt.backup_all_reachBatch=alltbt.all_reachBatch;
    alltbt.all_reachBatch=alltbt.(test.onlyConsiderReachType);
end
saveDir2=[saveDir '\' test.event_name];
mkdir(saveDir2);
save([saveDir2 '\test_settings.mat'],'test');
fakeCueInd=50;
skipCorrected=true;
% this function builds the dataset using the trial type sequences specified above
[dataset,correctedDistributions]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
if ~isempty(test.onlyConsiderReachType)
    alltbt.all_reachBatch=alltbt.backup_all_reachBatch;
end
[returnThis,returnThisRef]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rawReaching');
savePSTH.beginner.success=returnThis;

% DROP
alltbt.all_reachBatch=alltbt.reachBatch_drop_reachStarts;
test.nInSequence=[2]; % defines trial pairs, e.g., 2 means will compare each trial with its subsequent trial, 3 means will compare each trial with the trial after next, etc.
trial1='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1';
test.templateSequence2_cond=eval(trial1);
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=true; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name=['alltrials' tbt_filter.name 'inBetweenAnything' num2str(test.fillInBetweenWithAnything)];
test.onlyConsiderReachType=[];
if ~isempty(test.onlyConsiderReachType)
    alltbt.backup_all_reachBatch=alltbt.all_reachBatch;
    alltbt.all_reachBatch=alltbt.(test.onlyConsiderReachType);
end
saveDir2=[saveDir '\' test.event_name];
mkdir(saveDir2);
save([saveDir2 '\test_settings.mat'],'test');
fakeCueInd=50;
skipCorrected=true;
% this function builds the dataset using the trial type sequences specified above
[dataset,correctedDistributions]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
if ~isempty(test.onlyConsiderReachType)
    alltbt.all_reachBatch=alltbt.backup_all_reachBatch;
end
[returnThis,returnThisRef]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rawReaching');
savePSTH.beginner.drop=returnThis;

% MISS
alltbt.all_reachBatch=alltbt.reachBatch_miss_reachStarts;
test.nInSequence=[2]; % defines trial pairs, e.g., 2 means will compare each trial with its subsequent trial, 3 means will compare each trial with the trial after next, etc.
trial1='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1';
test.templateSequence2_cond=eval(trial1);
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=true; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name=['alltrials' tbt_filter.name 'inBetweenAnything' num2str(test.fillInBetweenWithAnything)];
test.onlyConsiderReachType=[];
if ~isempty(test.onlyConsiderReachType)
    alltbt.backup_all_reachBatch=alltbt.all_reachBatch;
    alltbt.all_reachBatch=alltbt.(test.onlyConsiderReachType);
end
saveDir2=[saveDir '\' test.event_name];
mkdir(saveDir2);
save([saveDir2 '\test_settings.mat'],'test');
fakeCueInd=50;
skipCorrected=true;
% this function builds the dataset using the trial type sequences specified above
[dataset,correctedDistributions]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
if ~isempty(test.onlyConsiderReachType)
    alltbt.all_reachBatch=alltbt.backup_all_reachBatch;
end
[returnThis,returnThisRef]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rawReaching');
savePSTH.beginner.miss=returnThis;





%%%%%%%%%%%%
alltbt=backup.alltbt;
trialTypes=backup.trialTypes;
metadata=backup.metadata;

% Filter by dprime
temp=datestr(datetime('now'));
temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))=''; 
temp=temp(~isspace(temp));
saveDir=['/Volumes/Neurobio/MICROSCOPE/Kim/RT pairs data sets/' temp]; % where to save details of alltbt filtering and RT pairs data set
% INTERMEDIATE
alltbt.mouseid=metadata.mouseid;
alltbt.sessid=metadata.sessid;
trialTypes.sessid=metadata.sessid;
tbt_filter.sortField='dprimes';
tbt_filter.range_values=[0.25 0.75]; % beginner: d<0.25, intermediate: 0.25<=d<0.75, expert: d>=0.75
tbt_filter.name=[tbt_filter.sortField num2str(tbt_filter.range_values(1)) 'to' num2str(tbt_filter.range_values(2))];
temp=tbt_filter.name;
temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))=''; 
temp=temp(~isspace(temp));
tbt_filter.name=temp;
tbt_filter.clock_progress=true;
% filter alltbt
[alltbt,trialTypes,metadata]=filtTbt(alltbt,trialTypes,tbt_filter.sortField,tbt_filter.range_values,metadata,tbt_filter.clock_progress);

% SUCCESS
alltbt.all_reachBatch=alltbt.reachBatch_success_reachStarts;
test.nInSequence=[2]; % defines trial pairs, e.g., 2 means will compare each trial with its subsequent trial, 3 means will compare each trial with the trial after next, etc.
trial1='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1';
test.templateSequence2_cond=eval(trial1);
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=true; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name=['alltrials' tbt_filter.name 'inBetweenAnything' num2str(test.fillInBetweenWithAnything)];
test.onlyConsiderReachType=[];
if ~isempty(test.onlyConsiderReachType)
    alltbt.backup_all_reachBatch=alltbt.all_reachBatch;
    alltbt.all_reachBatch=alltbt.(test.onlyConsiderReachType);
end
saveDir2=[saveDir '\' test.event_name];
mkdir(saveDir2);
save([saveDir2 '\test_settings.mat'],'test');
fakeCueInd=50;
skipCorrected=true;
% this function builds the dataset using the trial type sequences specified above
[dataset,correctedDistributions]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
if ~isempty(test.onlyConsiderReachType)
    alltbt.all_reachBatch=alltbt.backup_all_reachBatch;
end
[returnThis,returnThisRef]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rawReaching');
savePSTH.intermediate.success=returnThis;

% DROP
alltbt.all_reachBatch=alltbt.reachBatch_drop_reachStarts;
test.nInSequence=[2]; % defines trial pairs, e.g., 2 means will compare each trial with its subsequent trial, 3 means will compare each trial with the trial after next, etc.
trial1='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1';
test.templateSequence2_cond=eval(trial1);
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=true; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name=['alltrials' tbt_filter.name 'inBetweenAnything' num2str(test.fillInBetweenWithAnything)];
test.onlyConsiderReachType=[];
if ~isempty(test.onlyConsiderReachType)
    alltbt.backup_all_reachBatch=alltbt.all_reachBatch;
    alltbt.all_reachBatch=alltbt.(test.onlyConsiderReachType);
end
saveDir2=[saveDir '\' test.event_name];
mkdir(saveDir2);
save([saveDir2 '\test_settings.mat'],'test');
fakeCueInd=50;
skipCorrected=true;
% this function builds the dataset using the trial type sequences specified above
[dataset,correctedDistributions]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
if ~isempty(test.onlyConsiderReachType)
    alltbt.all_reachBatch=alltbt.backup_all_reachBatch;
end
[returnThis,returnThisRef]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rawReaching');
savePSTH.intermediate.drop=returnThis;

% MISS
alltbt.all_reachBatch=alltbt.reachBatch_miss_reachStarts;
test.nInSequence=[2]; % defines trial pairs, e.g., 2 means will compare each trial with its subsequent trial, 3 means will compare each trial with the trial after next, etc.
trial1='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1';
test.templateSequence2_cond=eval(trial1);
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=true; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name=['alltrials' tbt_filter.name 'inBetweenAnything' num2str(test.fillInBetweenWithAnything)];
test.onlyConsiderReachType=[];
if ~isempty(test.onlyConsiderReachType)
    alltbt.backup_all_reachBatch=alltbt.all_reachBatch;
    alltbt.all_reachBatch=alltbt.(test.onlyConsiderReachType);
end
saveDir2=[saveDir '\' test.event_name];
mkdir(saveDir2);
save([saveDir2 '\test_settings.mat'],'test');
fakeCueInd=50;
skipCorrected=true;
% this function builds the dataset using the trial type sequences specified above
[dataset,correctedDistributions]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
if ~isempty(test.onlyConsiderReachType)
    alltbt.all_reachBatch=alltbt.backup_all_reachBatch;
end
[returnThis,returnThisRef]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rawReaching');
savePSTH.intermediate.miss=returnThis;





%%%%%%%%%%%%
alltbt=backup.alltbt;
trialTypes=backup.trialTypes;
metadata=backup.metadata;

% Filter by dprime
temp=datestr(datetime('now'));
temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))=''; 
temp=temp(~isspace(temp));
saveDir=['/Volumes/Neurobio/MICROSCOPE/Kim/RT pairs data sets/' temp]; % where to save details of alltbt filtering and RT pairs data set
% EXPERT
alltbt.mouseid=metadata.mouseid;
alltbt.sessid=metadata.sessid;
trialTypes.sessid=metadata.sessid;
tbt_filter.sortField='dprimes';
tbt_filter.range_values=[0.75 100]; % beginner: d<0.25, intermediate: 0.25<=d<0.75, expert: d>=0.75
tbt_filter.name=[tbt_filter.sortField num2str(tbt_filter.range_values(1)) 'to' num2str(tbt_filter.range_values(2))];
temp=tbt_filter.name;
temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))=''; 
temp=temp(~isspace(temp));
tbt_filter.name=temp;
tbt_filter.clock_progress=true;
% filter alltbt
[alltbt,trialTypes,metadata]=filtTbt(alltbt,trialTypes,tbt_filter.sortField,tbt_filter.range_values,metadata,tbt_filter.clock_progress);

% SUCCESS
alltbt.all_reachBatch=alltbt.reachBatch_success_reachStarts;
test.nInSequence=[2]; % defines trial pairs, e.g., 2 means will compare each trial with its subsequent trial, 3 means will compare each trial with the trial after next, etc.
trial1='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1';
test.templateSequence2_cond=eval(trial1);
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=true; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name=['alltrials' tbt_filter.name 'inBetweenAnything' num2str(test.fillInBetweenWithAnything)];
test.onlyConsiderReachType=[];
if ~isempty(test.onlyConsiderReachType)
    alltbt.backup_all_reachBatch=alltbt.all_reachBatch;
    alltbt.all_reachBatch=alltbt.(test.onlyConsiderReachType);
end
saveDir2=[saveDir '\' test.event_name];
mkdir(saveDir2);
save([saveDir2 '\test_settings.mat'],'test');
fakeCueInd=50;
skipCorrected=true;
% this function builds the dataset using the trial type sequences specified above
[dataset,correctedDistributions]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
if ~isempty(test.onlyConsiderReachType)
    alltbt.all_reachBatch=alltbt.backup_all_reachBatch;
end
[returnThis,returnThisRef]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rawReaching');
savePSTH.expert.success=returnThis;

% DROP
alltbt.all_reachBatch=alltbt.reachBatch_drop_reachStarts;
test.nInSequence=[2]; % defines trial pairs, e.g., 2 means will compare each trial with its subsequent trial, 3 means will compare each trial with the trial after next, etc.
trial1='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1';
test.templateSequence2_cond=eval(trial1);
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=true; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name=['alltrials' tbt_filter.name 'inBetweenAnything' num2str(test.fillInBetweenWithAnything)];
test.onlyConsiderReachType=[];
if ~isempty(test.onlyConsiderReachType)
    alltbt.backup_all_reachBatch=alltbt.all_reachBatch;
    alltbt.all_reachBatch=alltbt.(test.onlyConsiderReachType);
end
saveDir2=[saveDir '\' test.event_name];
mkdir(saveDir2);
save([saveDir2 '\test_settings.mat'],'test');
fakeCueInd=50;
skipCorrected=true;
% this function builds the dataset using the trial type sequences specified above
[dataset,correctedDistributions]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
if ~isempty(test.onlyConsiderReachType)
    alltbt.all_reachBatch=alltbt.backup_all_reachBatch;
end
[returnThis,returnThisRef]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rawReaching');
savePSTH.expert.drop=returnThis;

% MISS
alltbt.all_reachBatch=alltbt.reachBatch_miss_reachStarts;
test.nInSequence=[2]; % defines trial pairs, e.g., 2 means will compare each trial with its subsequent trial, 3 means will compare each trial with the trial after next, etc.
trial1='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1';
test.templateSequence2_cond=eval(trial1);
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=true; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name=['alltrials' tbt_filter.name 'inBetweenAnything' num2str(test.fillInBetweenWithAnything)];
test.onlyConsiderReachType=[];
if ~isempty(test.onlyConsiderReachType)
    alltbt.backup_all_reachBatch=alltbt.all_reachBatch;
    alltbt.all_reachBatch=alltbt.(test.onlyConsiderReachType);
end
saveDir2=[saveDir '\' test.event_name];
mkdir(saveDir2);
save([saveDir2 '\test_settings.mat'],'test');
fakeCueInd=50;
skipCorrected=true;
% this function builds the dataset using the trial type sequences specified above
[dataset,correctedDistributions]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
if ~isempty(test.onlyConsiderReachType)
    alltbt.all_reachBatch=alltbt.backup_all_reachBatch;
end
[returnThis,returnThisRef]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rawReaching');
savePSTH.expert.miss=returnThis;