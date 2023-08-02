function [alltbt,unique_sessid,changeInDprimes]=findSessWhereMouseLearned(alltbt,metadata,trialTypes,part1_fracThroughSess,part2_fracThroughSess,learningThresh)

nInSequence=2;
% flankingTrials='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1'; % take every trial
flankingTrials='trialTypes.led==0'; % take no LED trials
shuffleTrialOrder=false; % if want to randomly permute trial order to test for ordering effects
alltbt.sessid=metadata.sessid;
trialTypes.sessid=metadata.sessid;
reachratesettings.epsilon_cue=0; % in seconds
reachratesettings.epsilon_uncue=2; % in seconds
reachratesettings.epsilon_beforecue=1; % in seconds
reachratesettings.percentOfReachesFromSess_forInitCond=20; % use this fraction of trials from beginning of session to get initial windows
reachratesettings.percentOfReachesFromSess_forInitRate=20; % use this fraction of trials from beginning of session to get initial reach rates
reachratesettings.maxTrialLength=9.5; % in sec, wrt cue
reachratesettings.minTrialLength=-2; % wrt cue, in sec
reachratesettings.suppressPlots=true;

dpset=settingsForDprimes(alltbt,'cueZone_onVoff',false);

reachratesettings.acrossSess_window1=[dpset.reachAfterCueWindow_start dpset.reachAfterCueWindow_end]; % cued window [0.05 1]
% reachratesettings.acrossSess_window1=[0 9.5]; % cued window [0.05 1]
% reachratesettings.acrossSess_window1=[4 7];
% note that after mouse gets a pellet, reaching is suppressed
reachratesettings.acrossSess_window2=[7 reachratesettings.maxTrialLength]; % beware reach suppression after a success
reachratesettings.acrossSess_window3=[dpset.preCueWindow_start1 dpset.preCueWindow_end1]; %[reachratesettings.minTrialLength 0]; 
reachratesettings.scatterPointSize=50; % size for points in scatter plot
reachratesettings.addSatietyLines=true; % whether to add proportionality lines to figure
reachratesettings.stopPlottingTrialsAfterN=500;
reachratesettings.showFitLine=false;
reachratesettings.useWindowsForUncued=[3]; % to use window2 or window3 or both for the uncued reach rate
reachratesettings.initWindows=[]; % empty if want to calculate from dataset
reachratesettings.addSessionLines=false; % for no averaging across sessions plot, whether to connect trial bins within same session with lines
reachratesettings.binThisManyTrials=25; % how many trials to bin within each session
reachratesettings.nBinsForZones=40; % will be nBinsForZones squared total bins, this is # bins for each x and y axis
reachratesettings.useRateMethod=3; % 1, 2 or 3 (see explanation below)

backup_alltbt=alltbt; backup_metadata=metadata; backup_trialTypes=trialTypes;
% FIRST PART OF SESSION
temp=datestr(datetime('now'));
temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))=''; 
temp=temp(~isspace(temp));
saveDir=['/Volumes/Neurobio/MICROSCOPE/Kim/RT pairs data sets/' temp]; % where to save details of alltbt filtering and RT pairs data set
tbt_filter.sortField='fractionThroughSess_adjusted';
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





alltbt=backup_alltbt; metadata=backup_metadata; trialTypes=backup_trialTypes;
% SECOND PART OF SESSION
temp=datestr(datetime('now'));
temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))=''; 
temp=temp(~isspace(temp));
saveDir=['/Volumes/Neurobio/MICROSCOPE/Kim/RT pairs data sets/' temp]; % where to save details of alltbt filtering and RT pairs data set
tbt_filter.sortField='fractionThroughSess_adjusted';
tbt_filter.range_values=part2_fracThroughSess;
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

matchesEventCond_trial_n=dataset.realDistributions.allTrialsSequence_isSeq{1}==1;
nth_sessions=metadata.sessid(matchesEventCond_trial_n);

[a,b]=getdprimes(reachrates);
if ~isempty(a)
    dprimes_lasttrial_secondPart=a;
    dprimes_firsttrial_secondPart=b;
end

alltbt=backup_alltbt; 
metadata=backup_metadata; 
trialTypes=backup_trialTypes;
unique_sessid=unique(metadata.sessid);
changeInDprimes=dprimes_lasttrial_secondPart-dprimes_lasttrial;

if length(changeInDprimes)~=length(unique_sessid)
    expand_changeInDprimes=nan(1,length(unique_sessid));
    % not an event in every session
    % need to fill missing sessions
    unique_nth_sessions=unique(nth_sessions);
    for i=1:length(unique_sessid)
        if ~ismember(unique_sessid(i),nth_sessions)
            % event did not happen in this session
            % already padded with nan
        else
            expand_changeInDprimes(i)=changeInDprimes(unique_nth_sessions==unique_sessid(i));
        end
    end
    changeInDprimes=expand_changeInDprimes;
end

alltbt.mouseLearned=nan(size(alltbt.sessid));
alltbt.withinSess_changeInDprime=nan(size(alltbt.sessid));
for i=1:length(unique_sessid)
    currsess=unique_sessid(i);
    alltbt.withinSess_changeInDprime(alltbt.sessid==currsess)=changeInDprimes(unique_sessid==currsess);
    alltbt.mouseLearned(alltbt.sessid==currsess)=changeInDprimes(unique_sessid==currsess)>learningThresh;
end

% plot histogram of dprimes changes within session
figure();
histogram(changeInDprimes,50);
hold on;
line([learningThresh learningThresh],[0 10],'Color','r');

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