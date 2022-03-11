function memoryEffect(alltbt,metadata,trialTypes,nInSeq,useFractionThroughSession)

nInSequence=2;
flankingTrials='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1'; % take every trial

trialTypes.isLongITI_1forward=[trialTypes.isLongITI(2:end); 0];
trialTypes.optoGroup_1forward=[trialTypes.optoGroup(2:end); 0];
temp=datestr(datetime('now'));
temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))=''; 
temp=temp(~isspace(temp));
saveDir=['/Volumes/Neurobio/MICROSCOPE/Kim/RT pairs data sets/' temp]; % where to save details of alltbt filtering and RT pairs data set
alltbt.sessid=metadata.sessid;
trialTypes.sessid=metadata.sessid;
shuffleTrialOrder=false; % if want to randomly permute trial order to test for ordering effects
reachratesettings=getReachRateSettings();

% FIRST PLOT CHANGE OVER ALL TRIALS WITHIN SESSION
trial1=flankingTrials; 
test.trial1=trial1;
trial2=flankingTrials;
test.trial2=trial2;
reachratesettings.
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir);
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected);
reachratesettings.suppressPlots=false;
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 

% THEN FILTER TBT 

backup_alltbt=alltbt; backup_metadata=metadata; backup_trialTypes=trialTypes;
% NO LED FIRST
test.nInSequence=[nInSeq]; % defines trial pairs, e.g., 2 means will compare each trial with its subsequent trial, 3 means will compare each trial with the trial after next, etc.

% memory
trial1='trialTypes.led~=1'; 
test.trial1=trial1;
test.templateSequence2_cond=eval(trial1);
% memory
trial2=['trialTypes.led~=1' ' & trialTypes.led_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.optoGroup_1forward~=3']; % & trialTypes.led_2forward==0'];
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


tbt_filter.sortField='fractionThroughSess';
tbt_filter.range_values=part1_fracThroughSess;
tbt_filter.name=[tbt_filter.sortField num2str(tbt_filter.range_values(1)) 'to' num2str(tbt_filter.range_values(2))];
temp=tbt_filter.name;
temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))=''; 
temp=temp(~isspace(temp));
tbt_filter.name=temp;
tbt_filter.clock_progress=true;
[alltbt,trialTypes,metadata]=filtTbt(alltbt,trialTypes,tbt_filter.sortField,tbt_filter.range_values,metadata,tbt_filter.clock_progress);

trial1=[flankingTrials];
trial2=[flankingTrials];
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir);
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
[a,b]=getdprimes(reachrates);
if ~isempty(a)
    dprimes_lasttrial=a;
    dprimes_firsttrial=b;
end

[a,b]=getdprimes(reachrates);
if ~isempty(a)
    dprimes_lasttrial_=a;
    dprimes_firsttrial=b;
else
    dprimes_lasttrial=[];
    dprimes_firsttrial=[];
end





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
trial2=['trialTypes.led~=1' ' & trialTypes.led_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.optoGroup_1forward~=3']; % & trialTypes.led_2forward==0'];
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
reachratesettings.acrossSess_window1=[0.05 2]; % cued window [0.05 1]
% reachratesettings.acrossSess_window1=[4 7];
% note that after mouse gets a pellet, reaching is suppressed
reachratesettings.acrossSess_window2=[7 reachratesettings.maxTrialLength]; % beware reach suppression after a success
reachratesettings.acrossSess_window3=[reachratesettings.minTrialLength 0];
reachratesettings.scatterPointSize=50; % size for points in scatter plot
reachratesettings.addSatietyLines=false; % whether to add proportionality lines to figure
reachratesettings.stopPlottingTrialsAfterN=500; % will stop plotting
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
trialTypes.led_4forward=[trialTypes.led(5:end); 0; 0; 0; 0];
linker=' & trialTypes.led_1forward~=1'; % & trialTypes.led_2forward~=1';
% linker=' & trialTypes.led_1forward==0 & trialTypes.led_2forward==0 & trialTypes.led_3forward==0 & trialTypes.led_4forward==0';
trial2=['trialTypes.led==1 & trialTypes.optoGroup~=1 & trialTypes.optoGroup~=3' linker];
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
LED_all_uncued=uncued_inThisPartOfSess(1:end);
LED_all_cued=cued_inThisPartOfSess(1:end);
LED_all_fracs=fracs_inThisPartOfSess_LED(1:end);
LED_all_uncued=LED_all_uncued(~isnan(LED_all_fracs));
LED_all_cued=LED_all_cued(~isnan(LED_all_fracs));
LED_all_fracs=LED_all_fracs(~isnan(LED_all_fracs));

binsForFracs=0:0.02:1;
figure();
[n,x]=histcounts(fracs_inThisPartOfSess_noLED(~isnan(fracs_inThisPartOfSess_noLED)),binsForFracs);
[new_n,new_x]=cityscape_hist(n,x);
plot(new_x,new_n,'Color','k');
con.x=x;
con.n=n;
hold on;
[n,x]=histcounts(fracs_inThisPartOfSess_LED(~isnan(fracs_inThisPartOfSess_LED)),binsForFracs);
[new_n,new_x]=cityscape_hist(n,x);
plot(new_x,new_n,'Color','r');
led.x=x;
led.n=n;
title('Before matching fracs');

resampleToMatchFracs=false;
if resampleToMatchFracs==true
    nRuns=10;
    resampled_LEDuncued=[];
    resampled_LEDcued=[];
    resampled_LEDfracs=[];
    for i=1:length(con.x)-1
        startingpool_LEDuncued=LED_all_uncued(LED_all_fracs>=binsForFracs(i) & LED_all_fracs<binsForFracs(i+1));
        startingpool_LEDcued=LED_all_cued(LED_all_fracs>=binsForFracs(i) & LED_all_fracs<binsForFracs(i+1));
        startingpool_LEDfracs=LED_all_fracs(LED_all_fracs>=binsForFracs(i) & LED_all_fracs<binsForFracs(i+1));
        takeN=con.n(i);
        if takeN~= 0 && isempty(startingpool_LEDuncued)
            disp('Cannot resample, not enough trials');
            disp('Proceeding without resample');
            resampleToMatchFracs=false;
            break
        end
        if takeN>length(startingpool_LEDuncued)
            newpool_LEDuncued=[];
            newpool_LEDcued=[];
            newpool_LEDfracs=[];
            nTimesMore=ceil(takeN/length(startingpool_LEDuncued));
            newpool_LEDuncued=repmat(startingpool_LEDuncued,1,nTimesMore);
            newpool_LEDcued=repmat(startingpool_LEDcued,1,nTimesMore);
            newpool_LEDfracs=repmat(startingpool_LEDfracs,1,nTimesMore);
        else
            newpool_LEDuncued=startingpool_LEDuncued;
            newpool_LEDcued=startingpool_LEDcued;
            newpool_LEDfracs=startingpool_LEDfracs;
        end
        if isempty(newpool_LEDuncued) && takeN>0
            error('Need bigger binning for fracs');
        end
        if takeN==0
            continue
        end
        for j=1:nRuns
            disp(j);
            r=randsample(length(newpool_LEDuncued),takeN);
            resampled_LEDuncued=[resampled_LEDuncued newpool_LEDuncued(r)];
            resampled_LEDcued=[resampled_LEDcued newpool_LEDcued(r)];
            resampled_LEDfracs=[resampled_LEDfracs newpool_LEDfracs(r)];
        end
    end
    if resampleToMatchFracs==true
        figure();
        [n,x]=histcounts(fracs_inThisPartOfSess_noLED(~isnan(fracs_inThisPartOfSess_noLED)),binsForFracs);
        [new_n,new_x]=cityscape_hist(n,x);
        plot(new_x,new_n./nanmax(new_n),'Color','k');
        con.x=x;
        con.n=n;
        hold on;
        [n,x]=histcounts(resampled_LEDfracs(~isnan(resampled_LEDfracs)),binsForFracs);
        [new_n,new_x]=cityscape_hist(n,x);
        plot(new_x,new_n./nanmax(new_n),'Color','r');
        title('After matching fracs');
        mean_uncued_LED=nanmean(resampled_LEDuncued);
        mean_cued_LED=nanmean(resampled_LEDcued);
    end
end


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

function reachratesettings=getReachRateSettings()

reachratesettings.epsilon_cue=0; % in seconds
reachratesettings.epsilon_uncue=2; % in seconds
reachratesettings.epsilon_beforecue=1; % in seconds
reachratesettings.percentOfReachesFromSess_forInitCond=20; % use this fraction of trials from beginning of session to get initial windows
reachratesettings.percentOfReachesFromSess_forInitRate=20; % use this fraction of trials from beginning of session to get initial reach rates
reachratesettings.maxTrialLength=9.5; % in sec, wrt cue
reachratesettings.minTrialLength=-2; % wrt cue, in sec
reachratesettings.suppressPlots=true;
reachratesettings.acrossSess_window1=[0.05 2]; % cued window [0.05 1]
reachratesettings.acrossSess_window2=[7 reachratesettings.maxTrialLength]; % beware reach suppression after a success
reachratesettings.acrossSess_window3=[reachratesettings.minTrialLength -1]; 
reachratesettings.scatterPointSize=50; % size for points in scatter plot
reachratesettings.addSatietyLines=false; % whether to add proportionality lines to figure
reachratesettings.stopPlottingTrialsAfterN=500;
reachratesettings.showFitLine=false;
reachratesettings.useWindowsForUncued=[3]; % to use window2 or window3 or both for the uncued reach rate
reachratesettings.initWindows=[]; % empty if want to calculate from dataset
reachratesettings.addSessionLines=false; % for no averaging across sessions plot, whether to connect trial bins within same session with lines
reachratesettings.binTrialsForAvAcrossSess=true; % whether to bin multiple trials for first figure, will bin into binThisManyTrials
reachratesettings.binThisManyTrials=6; % how many trials to bin within each session
reachratesettings.nBinsForZones=40; % will be nBinsForZones squared total bins, this is # bins for each x and y axis
reachratesettings.useRateMethod=3; 

end

function [test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir)

tbt_filter.name='throwaway';
test.nInSequence=[nInSequence]; 
test.trial1=trial1;
test.templateSequence2_cond=eval(trial1);
test.trial2=trial2;
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=true; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name=['alltrials' tbt_filter.name 'inBetweenAnything' num2str(test.fillInBetweenWithAnything)];
saveDir2=[saveDir '\' test.event_name];
mkdir(saveDir2);
save([saveDir2 '\test_settings.mat'],'test');
fakeCueInd=50; 
skipCorrected=true;

end

function [dprimes_lasttrial,dprimes_firsttrial]=getdprimes(reachrates)

dprimes_lasttrial=[];
dprimes_firsttrial=[];
if isempty(reachrates)
    return
end
dprimes_lasttrial=calc_dprime_per_sess(reachrates.alltrials_uncued,reachrates.alltrials_cued);
dprimes_firsttrial=calc_dprime_per_sess(reachrates.trial1_alltrials_uncued,reachrates.trial1_alltrials_cued);

end

function dprimes=calc_dprime_per_sess(uncued_events,cued_events)

hit_rates=nansum(cued_events>0,2)./nansum(~isnan(cued_events),2);
fa_rates=nansum(uncued_events>0,2)./nansum(~isnan(uncued_events),2);
% closest we can get to 1 or zero is defined by number of trials
ns=nansum(~isnan(cued_events),2);
hit_rates(ns<3)=nan;
fa_rates(ns<3)=nan;
hit_rates(hit_rates==1)=1-(1./ns(hit_rates==1));
hit_rates(hit_rates==0)=0+(1./ns(hit_rates==0));
fa_rates(fa_rates==1)=1-(1./ns(fa_rates==1));
fa_rates(fa_rates==0)=0+(1./ns(fa_rates==0));
dprimes=dprime(hit_rates,fa_rates);

end

function out=dprime(hit_rates,FA_rates)

out=norminv(hit_rates)-norminv(FA_rates);

end