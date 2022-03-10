function memoryEffect(alltbt,metadata,trialTypes,nInSeq,useFractionThroughSession)


temp=datestr(datetime('now'));
temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))=''; 
temp=temp(~isspace(temp));
saveDir=['/Volumes/Neurobio/MICROSCOPE/Kim/RT pairs data sets/' temp]; % where to save details of alltbt filtering and RT pairs data set

%%%%%%%%% NO LED FIRST

% settings for paired RT data set
test.nInSequence=[nInSeq]; % defines trial pairs, e.g., 2 means will compare each trial with its subsequent trial, 3 means will compare each trial with the trial after next, etc.
% requirement for first trial in pair
trialTypes.isLongITI_1forward=[trialTypes.isLongITI(2:end); 0];
trialTypes.optoGroup_1forward=[trialTypes.optoGroup(2:end); 0];
% memory
trial1='trialTypes.led~=1'; 
test.trial1=trial1;
test.templateSequence2_cond=eval(trial1);
% memory
trial2='trialTypes.led~=1';
test.trial2=trial2;
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=false; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
tbt_filter.name='temp';
test.event_name=['alltrials' tbt_filter.name 'inBetweenAnything' num2str(test.fillInBetweenWithAnything)];
% test.onlyConsiderReachType='reachBatch_success_reachStarts';
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
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
if ~isempty(test.onlyConsiderReachType)
    alltbt.all_reachBatch=alltbt.backup_all_reachBatch;
end

% plot trial to trial change in reach CDF
% plotChangeInReachCDF(dataset.realDistributions,alltbt); title('No LED');

% measure how reach rate changes over course of session
shuffleTrialOrder=false; % if want to randomly permute trial order to test for ordering effects
reachratesettings.epsilon_cue=0; % in seconds
reachratesettings.epsilon_uncue=2; % in seconds
reachratesettings.epsilon_beforecue=1; % in seconds
reachratesettings.percentOfReachesFromSess_forInitCond=20; % use this fraction of trials from beginning of session to get initial windows
reachratesettings.percentOfReachesFromSess_forInitRate=20; % use this fraction of trials from beginning of session to get initial reach rates
reachratesettings.maxTrialLength=9.5; % in sec, wrt cue
reachratesettings.minTrialLength=-2; % wrt cue, in sec
reachratesettings.suppressPlots=true;
 % sec wrt cue onset
reachratesettings.acrossSess_window1=[0.05 1]; % cued window [0.05 1]
% reachratesettings.acrossSess_window1=[4 7];
% note that after mouse gets a pellet, reaching is suppressed
reachratesettings.acrossSess_window2=[7 reachratesettings.maxTrialLength]; % beware reach suppression after a success
reachratesettings.acrossSess_window3=[reachratesettings.minTrialLength -1];
reachratesettings.scatterPointSize=50; % size for points in scatter plot
reachratesettings.addSatietyLines=false; % whether to add proportionality lines to figure
% reachratesettings.stopPlottingTrialsAfterN=500; % will stop plotting after this nth trial in session, also only use this many trials for regression fit -- see next line, also controls colormap
reachratesettings.stopPlottingTrialsAfterN=175; % will stop plotting
% after this nth trial in session, also only use this many trials for
% regression fit -- see next line, also controls colormap
reachratesettings.showFitLine=true; % whether to show linear fit to change across trials
reachratesettings.useWindowsForUncued=[3]; % to use window2 or window3 or both for the uncued reach rate
reachratesettings.initWindows=[]; % empty if want to calculate from dataset
reachratesettings.addSessionLines=false; % for no averaging across sessions plot, whether to connect trial bins within same session with lines
reachratesettings.binTrialsForAvAcrossSess=true; % whether to bin multiple trials for first figure, will bin into binThisManyTrials
reachratesettings.binThisManyTrials=6; % how many trials to bin within each session
reachratesettings.nBinsForZones=40; % will be nBinsForZones squared total bins, this is # bins for each x and y axis
reachratesettings.useRateMethod=3; % 1, 2 or 3 (see explanation below)
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 

% get average and s.e. reach rate for useFractionThroughSession
uncued_inThisPartOfSess=reachrates.alltrials_uncued;
cued_inThisPartOfSess=reachrates.alltrials_cued;
fracs_inThisPartOfSess_noLED=reachrates.fracsThroughSess;
for i=1:size(reachrates.alltrials_uncued,1)
    for j=1:size(reachrates.alltrials_uncued,2)
        if reachrates.fracsThroughSess(i,j)>=useFractionThroughSession(1) && reachrates.fracsThroughSess(i,j)<=useFractionThroughSession(2) 
        else
            uncued_inThisPartOfSess(i,j)=nan;
            cued_inThisPartOfSess(i,j)=nan;
            fracs_inThisPartOfSess_noLED(i,j)=nan;
        end
    end
end
mean_uncued=nanmean(uncued_inThisPartOfSess(1:end));
se_uncued=nanstd(uncued_inThisPartOfSess(1:end),[],2)./sqrt(nansum(~isnan(uncued_inThisPartOfSess(1:end))));
mean_cued=nanmean(cued_inThisPartOfSess(1:end));
se_cued=nanstd(cued_inThisPartOfSess(1:end),[],2)./sqrt(nansum(~isnan(cued_inThisPartOfSess(1:end))));




%%%%%%%%% LED

% settings for paired RT data set
test.nInSequence=[nInSeq]; % defines trial pairs, e.g., 2 means will compare each trial with its subsequent trial, 3 means will compare each trial with the trial after next, etc.
% requirement for first trial in pair
trialTypes.isLongITI_1forward=[trialTypes.isLongITI(2:end); 0];
trialTypes.optoGroup_1forward=[trialTypes.optoGroup(2:end); 0];
% memory
trial1='trialTypes.led~=1'; 
test.trial1=trial1;
test.templateSequence2_cond=eval(trial1);
% memory
trial2='trialTypes.led==1 & trialTypes.optoGroup~=1 & trialTypes.optoGroup~=3';
test.trial2=trial2;
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=false; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name=['alltrials' tbt_filter.name 'inBetweenAnything' num2str(test.fillInBetweenWithAnything)];
% test.onlyConsiderReachType='reachBatch_success_reachStarts';
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
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
if ~isempty(test.onlyConsiderReachType)
    alltbt.all_reachBatch=alltbt.backup_all_reachBatch;
end

% plot trial to trial change in reach CDF
% plotChangeInReachCDF(dataset.realDistributions,alltbt); title('No LED');

% measure how reach rate changes over course of session
shuffleTrialOrder=false; % if want to randomly permute trial order to test for ordering effects
reachratesettings.epsilon_cue=0; % in seconds
reachratesettings.epsilon_uncue=2; % in seconds
reachratesettings.epsilon_beforecue=1; % in seconds
reachratesettings.percentOfReachesFromSess_forInitCond=20; % use this fraction of trials from beginning of session to get initial windows
reachratesettings.percentOfReachesFromSess_forInitRate=20; % use this fraction of trials from beginning of session to get initial reach rates
reachratesettings.maxTrialLength=9.5; % in sec, wrt cue
reachratesettings.minTrialLength=-2; % wrt cue, in sec
reachratesettings.suppressPlots=true;
 % sec wrt cue onset
reachratesettings.acrossSess_window1=[0.05 1]; % cued window [0.05 1]
% reachratesettings.acrossSess_window1=[4 7];
% note that after mouse gets a pellet, reaching is suppressed
reachratesettings.acrossSess_window2=[7 reachratesettings.maxTrialLength]; % beware reach suppression after a success
reachratesettings.acrossSess_window3=[reachratesettings.minTrialLength -1];
reachratesettings.scatterPointSize=50; % size for points in scatter plot
reachratesettings.addSatietyLines=false; % whether to add proportionality lines to figure
% reachratesettings.stopPlottingTrialsAfterN=500; % will stop plotting after this nth trial in session, also only use this many trials for regression fit -- see next line, also controls colormap
reachratesettings.stopPlottingTrialsAfterN=175; % will stop plotting
% after this nth trial in session, also only use this many trials for
% regression fit -- see next line, also controls colormap
reachratesettings.showFitLine=true; % whether to show linear fit to change across trials
reachratesettings.useWindowsForUncued=[3]; % to use window2 or window3 or both for the uncued reach rate
reachratesettings.initWindows=[]; % empty if want to calculate from dataset
reachratesettings.addSessionLines=false; % for no averaging across sessions plot, whether to connect trial bins within same session with lines
reachratesettings.binTrialsForAvAcrossSess=true; % whether to bin multiple trials for first figure, will bin into binThisManyTrials
reachratesettings.binThisManyTrials=6; % how many trials to bin within each session
reachratesettings.nBinsForZones=40; % will be nBinsForZones squared total bins, this is # bins for each x and y axis
reachratesettings.useRateMethod=3; % 1, 2 or 3 (see explanation below)
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 

% get average and s.e. reach rate for useFractionThroughSession
uncued_inThisPartOfSess=reachrates.alltrials_uncued;
cued_inThisPartOfSess=reachrates.alltrials_cued;
fracs_inThisPartOfSess_LED=reachrates.fracsThroughSess;
for i=1:size(reachrates.alltrials_uncued,1)
    for j=1:size(reachrates.alltrials_uncued,2)
        if reachrates.fracsThroughSess(i,j)>=useFractionThroughSession(1) && reachrates.fracsThroughSess(i,j)<=useFractionThroughSession(2) 
        else
            uncued_inThisPartOfSess(i,j)=nan;
            cued_inThisPartOfSess(i,j)=nan;
            fracs_inThisPartOfSess_LED(i,j)=nan;
        end
    end
end
mean_uncued_LED=nanmean(uncued_inThisPartOfSess(1:end));
se_uncued_LED=nanstd(uncued_inThisPartOfSess(1:end),[],2)./sqrt(nansum(~isnan(uncued_inThisPartOfSess(1:end))));
mean_cued_LED=nanmean(cued_inThisPartOfSess(1:end));
se_cued_LED=nanstd(cued_inThisPartOfSess(1:end),[],2)./sqrt(nansum(~isnan(cued_inThisPartOfSess(1:end))));

figure();
[n,x]=hist(fracs_inThisPartOfSess_noLED(~isnan(fracs_inThisPartOfSess_noLED)),0:0.05:1);
[new_n,new_x]=cityscape_hist(n,x);
plot(new_x,new_n,'Color','k');
hold on;
[n,x]=hist(fracs_inThisPartOfSess_LED(~isnan(fracs_inThisPartOfSess_LED)),0:0.05:1);
[new_n,new_x]=cityscape_hist(n,x);
plot(new_x,new_n,'Color','r');

% plot output
figure();
scatter(mean_uncued,mean_cued,[],'k','filled');
hold on; 
scatter(mean_uncued_LED,mean_cued_LED,[],'r','filled');
line([mean_uncued-se_uncued mean_uncued+se_uncued],[mean_cued mean_cued],'Color','k');
line([mean_uncued mean_uncued],[mean_cued-se_uncued mean_cued+se_uncued],'Color','k');
line([mean_uncued_LED-se_uncued_LED mean_uncued_LED+se_uncued_LED],[mean_cued_LED mean_cued_LED],'Color','r');
line([mean_uncued_LED mean_uncued_LED],[mean_cued_LED-se_cued_LED mean_cued_LED+se_cued_LED],'Color','r');




end