function [f1,f2,returnout]=outcomeDependentShift_wrapper(alltbt,trialTypes,metadata,saveDir,f1,f2,reachratesettings,timeWindowOfEventReach,testEventReach,whichToPlot)

compareToFirstTrial=true;
linkSuccesses=false;
[~,cueindma]=nanmax(nanmean(alltbt.cueZone_onVoff,1));
timestep=mode(diff(nanmean(alltbt.times,1)));
if isnan(testEventReach.trial1)
    testEventReach.fillInBetweenWithAnything=true;
end

trialTypes.isLongITI_1forward=[trialTypes.isLongITI(2:end); 0];
trialTypes.optoGroup_1forward=[trialTypes.optoGroup(2:end); 0];
trialTypes.optoGroup_1back=[0; trialTypes.optoGroup(1:end-1)];
trialTypes.isLongITI_1back=[0; trialTypes.isLongITI(1:end-1)];
trialTypes.optoGroup_2back=[0; 0; trialTypes.optoGroup(1:end-2)];
trialTypes.noReach=~any(alltbt.all_reachBatch>0.05,2);
trialTypes.noReach_1forward=[trialTypes.noReach(2:end); 0];
trialTypes.noReach_1back=[0; trialTypes.noReach(1:end-1)];
trialTypes.reachedBeforeCue=any(alltbt.all_reachBatch(:,1:cueindma-1)>0.05,2);
trialTypes.reachedAfterCue=any(alltbt.all_reachBatch(:,cueindma:end)>0.05,2);
trialTypes.reachToPelletBeforeCue=any(alltbt.reachStarts_pelletPresent(:,1:cueindma-1)>0.05,2);
trialTypes.reachedBeforeCue_1forward=[trialTypes.reachedBeforeCue(2:end); 0];
trialTypes.reachToPelletBeforeCue_1forward=[trialTypes.reachToPelletBeforeCue(2:end); 0];
trialTypes.reachedAfterCue_1forward=[trialTypes.reachedAfterCue(2:end); 0];

% flankingTrials='trialTypes.optoGroup~=1 & trialTypes.chewing_at_trial_start==0';
% flankingTrials='trialTypes.optoGroup~=1';
flankingTrials='(trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1)';

% which to plot
% whichToPlot='success'; % can be 'success','delayed success','drop','cued touch','cued touch and switch color','failed cued reach','false alarm','no reach','basic','wildcard','backward success'
[plotset,trialTypes]=whichToPlotNow(whichToPlot,trialTypes,alltbt,timeWindowOfEventReach);

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
% reachratesettings.acrossSess_window1=[0.05 1]; % cued window [0.05 1]
% % reachratesettings.acrossSess_window1=[0 9.5]; % cued window [0.05 1]
% % reachratesettings.acrossSess_window1=[4 7];
% % note that after mouse gets a pellet, reaching is suppressed
% reachratesettings.acrossSess_window2=[7 reachratesettings.maxTrialLength]; % beware reach suppression after a success
% reachratesettings.acrossSess_window3=[reachratesettings.minTrialLength -1]; 
reachratesettings.scatterPointSize=50; % size for points in scatter plot
reachratesettings.addSatietyLines=true; % whether to add proportionality lines to figure
reachratesettings.stopPlottingTrialsAfterN=500;
reachratesettings.showFitLine=false;
% reachratesettings.useWindowsForUncued=[3]; % to use window2 or window3 or both for the uncued reach rate
reachratesettings.initWindows=[]; % empty if want to calculate from dataset
reachratesettings.addSessionLines=false; % for no averaging across sessions plot, whether to connect trial bins within same session with lines
reachratesettings.binThisManyTrials=25; % how many trials to bin within each session
reachratesettings.nBinsForZones=40; % will be nBinsForZones squared total bins, this is # bins for each x and y axis
reachratesettings.useRateMethod=3; % 1, 2 or 3 (see explanation below)
reachratesettings.binTrialsForAvAcrossSess=false; 
% see script_for_reaching_rate_analysis.m for explanation of rate methods
dprimes_noLED_lasttrial=[];
dprimes_LED_lasttrial=[];

% print settings
disp(['Using as time window to classify reach as cued (wrt cue onset in sec): ' num2str(reachratesettings.acrossSess_window1(1)) ' to ' num2str(reachratesettings.acrossSess_window1(2))]);
disp(['Using as time window to classify reach as UNCUED (wrt cue onset in sec): ' num2str(reachratesettings.acrossSess_window3(1)) ' to ' num2str(reachratesettings.acrossSess_window3(2))]);
pause;

if ~isempty(f1)
    set(0,'CurrentFigure',f1);
else
    f1=figure();
end
xlabel('Uncued reach rate (1/sec)');
ylabel('Depends on approach, if 3, then Cued reach rate (1/sec)');
title(['Approach ' num2str(reachratesettings.useWindowsForUncued)]);
% baseline all trials
nInSequence=3;
trial1=[flankingTrials ' & trialTypes.led_1forward==0 & trialTypes.optoGroup_1forward~=1']; % & trialTypes.isLongITI_1forward==1'];
trial2=[flankingTrials];
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir,true);
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
if isempty(reachrates)
    disp('No trials matching this criterion');
    return
end
if compareToFirstTrial==false
    [baseEffect_uncued_mean_out,baseEffect_cued_mean_out]=bootstrap(reachrates.alltrials_uncued,reachrates.alltrials_cued,'k','k',true);
end
[a,b]=doPlottingAndBootstrap(reachrates,[0.5 0.5 0.5],[0.5 0.5 0.5],'k',plotset.wildcard,compareToFirstTrial);
if ~isempty(a)
    dprimes_noLED_lasttrial=a;
    dprimes_noLED_firsttrial=b;
end
if plotset.wildcard==true
    returnout.reachrates_noLED=reachrates;
    returnout.dprimes_noLED_firsttrial=dprimes_noLED_firsttrial;
    returnout.dprimes_noLED_lasttrial=dprimes_noLED_lasttrial;
end
% success
nInSequence=3;
trial1=[flankingTrials ' & trialTypes.consumed_pellet_1back==1' ' & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.success_in_cued_window_1forward==1 & trialTypes.consumed_pellet_1forward==1 & trialTypes.led_1forward==0']; % & trialTypes.optoGroup_1forward~=1']; % & trialTypes.isLongITI_1forward==1'];
if linkSuccesses==false
    trial2=[flankingTrials];
else
    trial2=[flankingTrials ' & trialTypes.led_1forward==1 & trialTypes.optoGroup_1forward~=1'];
end
if plotset.success==true & ~isnan(testEventReach.trial1)
    disp('Using passed in test event conditions'); pause;
    trial1=testEventReach.trial1; trial2=testEventReach.trial2;
    if testEventReach.fillInBetweenWithAnything==false
        disp(['Doing nInSequence ' num2str(testEventReach.nInSequence)]); pause;
        nInSequence=testEventReach.nInSequence;
    end
end
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir,testEventReach.fillInBetweenWithAnything);
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
[a,b]=doPlottingAndBootstrap(reachrates,'g','g','k',plotset.success,compareToFirstTrial);
if ~isempty(a)
    dprimes_noLED_lasttrial=a;
    dprimes_noLED_firsttrial=b;
    noLED_rr=reachrates;
end
if plotset.success==true
    returnout.reachrates_noLED=reachrates;
    returnout.dprimes_noLED_firsttrial=dprimes_noLED_firsttrial;
    returnout.dprimes_noLED_lasttrial=dprimes_noLED_lasttrial;
end
% backwards success
nInSequence=3;
if linkSuccesses==false
    trial1=[flankingTrials];
else
    trial1=[flankingTrials ' & trialTypes.led_1back==1 & trialTypes.optoGroup_1back~=1'];
end
trial2=[flankingTrials ' & trialTypes.consumed_pellet_1forward==1' ' & trialTypes.reachedInTimeWindow_1back==1 & trialTypes.success_in_cued_window_1back==1 & trialTypes.consumed_pellet_1back==1 & trialTypes.led_1back==0 & trialTypes.optoGroup_1back~=1']; % & trialTypes.isLongITI_1forward==1'];
if plotset.backward_success==true & ~isnan(testEventReach.trial1)
    disp('Using passed in test event conditions'); pause;
    trial1=testEventReach.trial1; trial2=testEventReach.trial2;
    if testEventReach.fillInBetweenWithAnything==false
        disp(['Doing nInSequence ' num2str(testEventReach.nInSequence)]); pause;
        nInSequence=testEventReach.nInSequence;
    end
end
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir,testEventReach.fillInBetweenWithAnything);
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected);
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings);
[a,b]=doPlottingAndBootstrap(reachrates,'g','g','b',plotset.backward_success,compareToFirstTrial);
if ~isempty(a)
    dprimes_noLED_lasttrial=a;
    dprimes_noLED_firsttrial=b;
    noLED_rr=reachrates;
end
if plotset.backward_success==true
    returnout.reachrates_noLED=reachrates;
    returnout.dprimes_noLED_firsttrial=dprimes_noLED_firsttrial;
    returnout.dprimes_noLED_lasttrial=dprimes_noLED_lasttrial;
end
% delayed success
nInSequence=3;
trial1=[flankingTrials ' & trialTypes.consumed_pellet_1back==1' ' & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.after_cue_success_1forward==1 & trialTypes.consumed_pellet_1forward==1 & trialTypes.led_1forward==0 & trialTypes.optoGroup_1forward~=1']; % & trialTypes.isLongITI_1forward==1'];
if linkSuccesses==false
    trial2=[flankingTrials];
else
    trial2=[flankingTrials ' & trialTypes.led_1forward==1 & trialTypes.optoGroup_1forward~=1'];
end
if plotset.delayed==true & ~isnan(testEventReach.trial1)
    disp('Using passed in test event conditions'); pause;
    trial1=testEventReach.trial1; trial2=testEventReach.trial2;
    if testEventReach.fillInBetweenWithAnything==false
        disp(['Doing nInSequence ' num2str(testEventReach.nInSequence)]); pause;
        nInSequence=testEventReach.nInSequence;
    end
end
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir,testEventReach.fillInBetweenWithAnything);
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
[a,b]=doPlottingAndBootstrap(reachrates,[15 141 6]./255,[15 141 6]./255,'k',plotset.delayed,compareToFirstTrial);
if ~isempty(a)
    dprimes_noLED_lasttrial=a;
    dprimes_noLED_firsttrial=b;
    noLED_rr=reachrates;
end
if plotset.delayed==true
    returnout.reachrates_noLED=reachrates;
    returnout.dprimes_noLED_firsttrial=dprimes_noLED_firsttrial;
    returnout.dprimes_noLED_lasttrial=dprimes_noLED_lasttrial;
end
% drop
nInSequence=3;
trial1=[flankingTrials ' & trialTypes.consumed_pellet_1back==0' ' & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touch_in_cued_window_1forward==1 & trialTypes.consumed_pellet_1forward==0 & trialTypes.led_1forward==0 & trialTypes.optoGroup_1forward~=1']; % & trialTypes.isLongITI_1forward==1'];
trial2=[flankingTrials];
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir,testEventReach.fillInBetweenWithAnything);
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
[a,b]=doPlottingAndBootstrap(reachrates,'r','r','k',plotset.drop,compareToFirstTrial);
if ~isempty(a)
    dprimes_noLED_lasttrial=a;
    dprimes_noLED_firsttrial=b;
    noLED_rr=reachrates;
end
if plotset.drop==true
    returnout.reachrates_noLED=reachrates;
    returnout.dprimes_noLED_firsttrial=dprimes_noLED_firsttrial;
    returnout.dprimes_noLED_lasttrial=dprimes_noLED_lasttrial;
end
% touched after cue, i.e., success or drop
nInSequence=3;
trial1=[flankingTrials ' & trialTypes.touched_pellet_1back==1' ' & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touch_in_cued_window_1forward==1 & trialTypes.led_1forward==0 & trialTypes.optoGroup_1forward~=1']; % & trialTypes.isLongITI_1forward==1'];
trial2=[flankingTrials];
if plotset.cuedtouch==true & ~isnan(testEventReach.trial1)
    disp('Using passed in test event conditions'); pause;
    trial1=testEventReach.trial1; trial2=testEventReach.trial2;
    if testEventReach.fillInBetweenWithAnything==false
        disp(['Doing nInSequence ' num2str(testEventReach.nInSequence)]); pause;
        nInSequence=testEventReach.nInSequence;
    end
end
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir,testEventReach.fillInBetweenWithAnything);
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
[a,b]=doPlottingAndBootstrap(reachrates,plotset.cuedtouchcolor,plotset.cuedtouchcolor,'k',plotset.cuedtouch,compareToFirstTrial);
if ~isempty(a)
    dprimes_noLED_lasttrial=a;
    dprimes_noLED_firsttrial=b;
    noLED_rr=reachrates;
end
if plotset.cuedtouch==true
    returnout.reachrates_noLED=reachrates;
    returnout.dprimes_noLED_firsttrial=dprimes_noLED_firsttrial;
    returnout.dprimes_noLED_lasttrial=dprimes_noLED_lasttrial;
end
% did not touch after cue despite reaching
nInSequence=3;
trial1=[flankingTrials ' & trialTypes.touched_pellet_1back==0 & trialTypes.noReach_1back==0' ' & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touch_in_cued_window_1forward==0 & trialTypes.led_1forward==0 & trialTypes.optoGroup_1forward~=1']; % & trialTypes.isLongITI_1forward==1'];
if linkSuccesses==false
    trial2=[flankingTrials];
else
    trial2=[flankingTrials ' & trialTypes.led_1forward==1 & trialTypes.optoGroup_1forward~=1'];
end
if plotset.failedcued==true & ~isnan(testEventReach.trial1)
    disp('Using passed in test event conditions'); pause;
    trial1=testEventReach.trial1; trial2=testEventReach.trial2;
    if testEventReach.fillInBetweenWithAnything==false
        disp(['Doing nInSequence ' num2str(testEventReach.nInSequence)]); pause;
        nInSequence=testEventReach.nInSequence;
    end
end
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir,testEventReach.fillInBetweenWithAnything);
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
[a,b]=doPlottingAndBootstrap(reachrates,[8 41 175]./255,[8 41 175]./255,'k',plotset.failedcued,compareToFirstTrial);
if ~isempty(a)
    dprimes_noLED_lasttrial=a;
    dprimes_noLED_firsttrial=b;
    noLED_rr=reachrates;
end
if plotset.failedcued==true
    returnout.reachrates_noLED=reachrates;
    returnout.dprimes_noLED_firsttrial=dprimes_noLED_firsttrial;
    returnout.dprimes_noLED_lasttrial=dprimes_noLED_lasttrial;
end
% miss before cue
nInSequence=3;
trial1=['trialTypes.reachedBeforeCue_1forward==1 & trialTypes.reachToPelletBeforeCue_1forward==0 & trialTypes.led_1forward==0 & trialTypes.optoGroup_2back~=1']; % & trialTypes.isLongITI_1forward==1'];
trial2=[flankingTrials];
if plotset.falsealarm==true & ~isnan(testEventReach.trial1)
    disp('Using passed in test event conditions'); pause;
    trial1=testEventReach.trial1; trial2=testEventReach.trial2;
    if testEventReach.fillInBetweenWithAnything==false
        disp(['Doing nInSequence ' num2str(testEventReach.nInSequence)]); pause;
        nInSequence=testEventReach.nInSequence;
    end
end
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir,testEventReach.fillInBetweenWithAnything);
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
[a,b]=doPlottingAndBootstrap(reachrates,'c','c','k',plotset.falsealarm,compareToFirstTrial);
if ~isempty(a)
    dprimes_noLED_lasttrial=a;
    dprimes_noLED_firsttrial=b;
    noLED_rr=reachrates;
end
if plotset.falsealarm==true
    returnout.reachrates_noLED=reachrates;
    returnout.dprimes_noLED_firsttrial=dprimes_noLED_firsttrial;
    returnout.dprimes_noLED_lasttrial=dprimes_noLED_lasttrial;
end
% does not reach
nInSequence=3;
trial1=[flankingTrials ' & trialTypes.noReach_1forward==1 & trialTypes.touched_pellet_1forward==0 & trialTypes.led_1forward==0 & trialTypes.optoGroup_1forward~=1']; % & trialTypes.isLongITI_1forward==1'];
trial2=[flankingTrials];
if plotset.noreach==true & ~isnan(testEventReach.trial1)
    disp('Using passed in test event conditions'); pause;
    trial1=testEventReach.trial1; trial2=testEventReach.trial2;
    if testEventReach.fillInBetweenWithAnything==false
        disp(['Doing nInSequence ' num2str(testEventReach.nInSequence)]); pause;
        nInSequence=testEventReach.nInSequence;
    end
end
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir,testEventReach.fillInBetweenWithAnything);
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
[a,b]=doPlottingAndBootstrap(reachrates,[0.8 0.8 0.8],[0.8 0.8 0.8],'k',plotset.noreach,compareToFirstTrial);
if ~isempty(a)
    dprimes_noLED_lasttrial=a;
    dprimes_noLED_firsttrial=b;
    noLED_rr=reachrates;
end
if plotset.noreach==true
    returnout.reachrates_noLED=reachrates;
    returnout.dprimes_noLED_firsttrial=dprimes_noLED_firsttrial;
    returnout.dprimes_noLED_lasttrial=dprimes_noLED_lasttrial;
end

if ~strcmp(whichToPlot,'basic')
    if ~isempty(f2)
        set(0,'CurrentFigure',f2);
    else
        f2=figure();
    end
end

%%%% LEDs
% baseline all trials
nInSequence=3;
trial1=[flankingTrials ' & trialTypes.led_1forward==1 & trialTypes.optoGroup_1forward~=1']; % & trialTypes.isLongITI_1forward==1'];
trial2=[flankingTrials];
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir,true);
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
if compareToFirstTrial==false
    [baseEffect_uncued_mean_out,baseEffect_cued_mean_out]=bootstrap(reachrates.alltrials_uncued,reachrates.alltrials_cued,'k','m',true);
end
[a,b]=doPlottingAndBootstrap(reachrates,[0.5 0.5 0.5],'m','k',plotset.wildcard,compareToFirstTrial);
if ~isempty(a)
    dprimes_LED_lasttrial=a;
    dprimes_LED_firsttrial=b;
    LED_rr=reachrates;
end
if plotset.wildcard==true
    returnout.reachrates_LED=reachrates;
    returnout.dprimes_LED_firsttrial=dprimes_LED_firsttrial;
    returnout.dprimes_LED_lasttrial=dprimes_LED_lasttrial;
end
% success
nInSequence=3;
trial1=[flankingTrials ' & trialTypes.consumed_pellet_1back==1' ' & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.success_in_cued_window_1forward==1 & trialTypes.consumed_pellet_1forward==1 & trialTypes.led_1forward==1 & (trialTypes.optoGroup_1forward==3)']; % & trialTypes.optoGroup_1forward~=1 & trialTypes.optoGroup_1forward==2']; % & trialTypes.isLongITI_1forward==1'];
if linkSuccesses==false
    trial2=[flankingTrials];
else
    trial2=[flankingTrials ' & trialTypes.led_1forward==0'];
end
if plotset.success==true & ~isnan(testEventReach.trial1)
    if ~isfield(testEventReach,'trial1_LED')
        error('When passing in testEventReach conditions, need to receive both no LED and LED conditions');
    end
    trial1=testEventReach.trial1_LED; trial2=testEventReach.trial2_LED;
    if testEventReach.fillInBetweenWithAnything==false
        disp(['Doing nInSequence ' num2str(testEventReach.nInSequence)]); pause;
        nInSequence=testEventReach.nInSequence;
    end
end
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir,testEventReach.fillInBetweenWithAnything);
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
[a,b]=doPlottingAndBootstrap(reachrates,'g','m','k',plotset.success,compareToFirstTrial);
if ~isempty(a)
    dprimes_LED_lasttrial=a;
    dprimes_LED_firsttrial=b;
    LED_rr=reachrates;
end
if plotset.success==true
    returnout.reachrates_LED=reachrates;
    returnout.dprimes_LED_firsttrial=dprimes_LED_firsttrial;
    returnout.dprimes_LED_lasttrial=dprimes_LED_lasttrial;
end
% backwards success
nInSequence=3;
if linkSuccesses==false
    trial1=[flankingTrials];
else
    trial1=[flankingTrials ' & trialTypes.led_1back==0'];
end
trial2=[flankingTrials ' & trialTypes.consumed_pellet_1forward==1' ' & trialTypes.reachedInTimeWindow_1back==1 & trialTypes.success_in_cued_window_1back==1 & trialTypes.consumed_pellet_1back==1 & trialTypes.led_1back==1 & trialTypes.optoGroup_1back~=1']; % & trialTypes.isLongITI_1forward==1'];];
if plotset.backward_success==true & ~isnan(testEventReach.trial1)
    if ~isfield(testEventReach,'trial1_LED')
        error('When passing in testEventReach conditions, need to receive both no LED and LED conditions');
    end
    trial1=testEventReach.trial1_LED; trial2=testEventReach.trial2_LED;
    if testEventReach.fillInBetweenWithAnything==false
        disp(['Doing nInSequence ' num2str(testEventReach.nInSequence)]); pause;
        nInSequence=testEventReach.nInSequence;
    end
end
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir,testEventReach.fillInBetweenWithAnything);
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected);
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings);
[a,b]=doPlottingAndBootstrap(reachrates,'g','m','b',plotset.backward_success,compareToFirstTrial);
if ~isempty(a)
    dprimes_noLED_lasttrial=a;
    dprimes_noLED_firsttrial=b;
    LED_rr=reachrates;
end
if plotset.backward_success==true
    returnout.reachrates_LED=reachrates;
    returnout.dprimes_LED_firsttrial=dprimes_LED_firsttrial;
    returnout.dprimes_LED_lasttrial=dprimes_LED_lasttrial;
end
% delayed success
nInSequence=3;
trial1=[flankingTrials ' & trialTypes.consumed_pellet_1back==1' ' & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.after_cue_success_1forward==1 & trialTypes.consumed_pellet_1forward==1 & trialTypes.led_1forward==1 & trialTypes.optoGroup_1forward~=1']; % & trialTypes.isLongITI_1forward==1'];
if linkSuccesses==false
    trial2=[flankingTrials];
else
    trial2=[flankingTrials ' & trialTypes.led_1forward==0'];
end
if plotset.delayed==true & ~isnan(testEventReach.trial1)
    if ~isfield(testEventReach,'trial1_LED')
        error('When passing in testEventReach conditions, need to receive both no LED and LED conditions');
    end
    trial1=testEventReach.trial1_LED; trial2=testEventReach.trial2_LED;
    if testEventReach.fillInBetweenWithAnything==false
        disp(['Doing nInSequence ' num2str(testEventReach.nInSequence)]); pause;
        nInSequence=testEventReach.nInSequence;
    end
end
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir,testEventReach.fillInBetweenWithAnything);
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
[a,b]=doPlottingAndBootstrap(reachrates,[15 141 6]./255,'m','k',plotset.delayed,compareToFirstTrial);
if ~isempty(a)
    dprimes_LED_lasttrial=a;
    dprimes_LED_firsttrial=b;
    LED_rr=reachrates;
end
if plotset.delayed==true
    returnout.reachrates_LED=reachrates;
    returnout.dprimes_LED_firsttrial=dprimes_LED_firsttrial;
    returnout.dprimes_LED_lasttrial=dprimes_LED_lasttrial;
end
% drop
nInSequence=3;
trial1=[flankingTrials ' & trialTypes.consumed_pellet_1back==0' ' & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touch_in_cued_window_1forward==1 & trialTypes.consumed_pellet_1forward==0 & trialTypes.led_1forward==1 & trialTypes.optoGroup_1forward~=1']; % & trialTypes.isLongITI_1forward==1'];
trial2=[flankingTrials];
if plotset.drop==true & ~isnan(testEventReach.trial1)
    if ~isfield(testEventReach,'trial1_LED')
        error('When passing in testEventReach conditions, need to receive both no LED and LED conditions');
    end
    trial1=testEventReach.trial1_LED; trial2=testEventReach.trial2_LED;
    if testEventReach.fillInBetweenWithAnything==false
        disp(['Doing nInSequence ' num2str(testEventReach.nInSequence)]); pause;
        nInSequence=testEventReach.nInSequence;
    end
end
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir,testEventReach.fillInBetweenWithAnything);
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings);
[a,b]=doPlottingAndBootstrap(reachrates,'r','m','k',plotset.drop,compareToFirstTrial);
if ~isempty(a)
    dprimes_LED_lasttrial=a;
    dprimes_LED_firsttrial=b;
    LED_rr=reachrates;
end
if plotset.drop==true
    returnout.reachrates_LED=reachrates;
    returnout.dprimes_LED_firsttrial=dprimes_LED_firsttrial;
    returnout.dprimes_LED_lasttrial=dprimes_LED_lasttrial;
end
% touched after cue, i.e., success or drop
nInSequence=3;
trial1=[flankingTrials ' & trialTypes.touched_pellet_1back==1' ' & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touch_in_cued_window_1forward==1 & trialTypes.led_1forward==1 & trialTypes.optoGroup_1forward~=1']; % & trialTypes.isLongITI_1forward==1'];
trial2=[flankingTrials];
if plotset.cuedtouch==true & ~isnan(testEventReach.trial1)
    if ~isfield(testEventReach,'trial1_LED')
        error('When passing in testEventReach conditions, need to receive both no LED and LED conditions');
    end
    trial1=testEventReach.trial1_LED; trial2=testEventReach.trial2_LED;
    if testEventReach.fillInBetweenWithAnything==false
        disp(['Doing nInSequence ' num2str(testEventReach.nInSequence)]); pause;
        nInSequence=testEventReach.nInSequence;
    end
end
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir,testEventReach.fillInBetweenWithAnything);
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
[a,b]=doPlottingAndBootstrap(reachrates,plotset.cuedtouchcolor,'m','k',plotset.cuedtouch,compareToFirstTrial);
if ~isempty(a)
    dprimes_LED_lasttrial=a;
    dprimes_LED_firsttrial=b;
    LED_rr=reachrates;
end
if plotset.cuedtouch==true
    returnout.reachrates_LED=reachrates;
    returnout.dprimes_LED_firsttrial=dprimes_LED_firsttrial;
    returnout.dprimes_LED_lasttrial=dprimes_LED_lasttrial;
end
% did not touch after cue despite reaching
nInSequence=3;
trial1=[flankingTrials ' & trialTypes.touched_pellet_1back==0 & trialTypes.noReach_1back==0' ' & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touch_in_cued_window_1forward==0 & trialTypes.led_1forward==1 & trialTypes.optoGroup_1forward~=1']; % & trialTypes.isLongITI_1forward==1'];
if linkSuccesses==false
    trial2=[flankingTrials];
else
    trial2=[flankingTrials ' & trialTypes.led_1forward==0'];
end
if plotset.failedcued==true & ~isnan(testEventReach.trial1)
    if ~isfield(testEventReach,'trial1_LED')
        error('When passing in testEventReach conditions, need to receive both no LED and LED conditions');
    end
    trial1=testEventReach.trial1_LED; trial2=testEventReach.trial2_LED;
    if testEventReach.fillInBetweenWithAnything==false
        disp(['Doing nInSequence ' num2str(testEventReach.nInSequence)]); pause;
        nInSequence=testEventReach.nInSequence;
    end
end
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir,testEventReach.fillInBetweenWithAnything);
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
[a,b]=doPlottingAndBootstrap(reachrates,[8 41 175]./255,'m','k',plotset.failedcued,compareToFirstTrial);
if ~isempty(a)
    dprimes_LED_lasttrial=a;
    dprimes_LED_firsttrial=b;
    LED_rr=reachrates;
end
if plotset.failedcued==true
    returnout.reachrates_LED=reachrates;
    returnout.dprimes_LED_firsttrial=dprimes_LED_firsttrial;
    returnout.dprimes_LED_lasttrial=dprimes_LED_lasttrial;
end
% miss before cue
nInSequence=3;
trial1=['trialTypes.reachedBeforeCue_1forward==1 & trialTypes.reachToPelletBeforeCue_1forward==0 & trialTypes.led_1forward==1 & trialTypes.optoGroup_1forward~=1']; % & trialTypes.isLongITI_1forward==1'];
trial2=[flankingTrials];
if plotset.falsealarm==true & ~isnan(testEventReach.trial1)
    if ~isfield(testEventReach,'trial1_LED')
        error('When passing in testEventReach conditions, need to receive both no LED and LED conditions');
    end
    trial1=testEventReach.trial1_LED; trial2=testEventReach.trial2_LED;
    if testEventReach.fillInBetweenWithAnything==false
        disp(['Doing nInSequence ' num2str(testEventReach.nInSequence)]); pause;
        nInSequence=testEventReach.nInSequence;
    end
end
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir,testEventReach.fillInBetweenWithAnything);
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings);
[a,b]=doPlottingAndBootstrap(reachrates,'c','m','k',plotset.falsealarm,compareToFirstTrial);
if ~isempty(a)
    dprimes_LED_lasttrial=a;
    dprimes_LED_firsttrial=b;
    LED_rr=reachrates;
end
if plotset.falsealarm==true
    returnout.reachrates_LED=reachrates;
    returnout.dprimes_LED_firsttrial=dprimes_LED_firsttrial;
    returnout.dprimes_LED_lasttrial=dprimes_LED_lasttrial;
end
% does not reach
nInSequence=3;
trial1=[flankingTrials ' & trialTypes.noReach_1forward==1 & trialTypes.touched_pellet_1forward==0 & trialTypes.led_1forward==1 & trialTypes.optoGroup_1forward~=1']; % & trialTypes.isLongITI_1forward==1'];
trial2=[flankingTrials];
if plotset.noreach==true & ~isnan(testEventReach.trial1)
    if ~isfield(testEventReach,'trial1_LED')
        error('When passing in testEventReach conditions, need to receive both no LED and LED conditions');
    end
    trial1=testEventReach.trial1_LED; trial2=testEventReach.trial2_LED;
    if testEventReach.fillInBetweenWithAnything==false
        disp(['Doing nInSequence ' num2str(testEventReach.nInSequence)]); pause;
        nInSequence=testEventReach.nInSequence;
    end
end
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir,testEventReach.fillInBetweenWithAnything);
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
[a,b]=doPlottingAndBootstrap(reachrates,[0.8 0.8 0.8],'m','k',plotset.noreach,compareToFirstTrial);
if ~isempty(a)
    dprimes_LED_lasttrial=a;
    dprimes_LED_firsttrial=b;
    LED_rr=reachrates;
end
if plotset.noreach==true
    returnout.reachrates_LED=reachrates;
    returnout.dprimes_LED_firsttrial=dprimes_LED_firsttrial;
    returnout.dprimes_LED_lasttrial=dprimes_LED_lasttrial;
end

if ~strcmp(whichToPlot,'basic') && ~isempty(dprimes_noLED_lasttrial) && ~isempty(dprimes_LED_lasttrial)
    mi=nanmin([dprimes_noLED_lasttrial-dprimes_noLED_firsttrial; dprimes_LED_lasttrial-dprimes_LED_firsttrial]);
    ma=nanmax([dprimes_noLED_lasttrial-dprimes_noLED_firsttrial; dprimes_LED_lasttrial-dprimes_LED_firsttrial]);
    binstep=0.3;
    anyisnan=isnan(dprimes_noLED_lasttrial) | isnan(dprimes_LED_lasttrial);
    dprimes_noLED_lasttrial(anyisnan)=nan;
    dprimes_noLED_firsttrial(anyisnan)=nan;
    dprimes_LED_lasttrial(anyisnan)=nan;
    dprimes_LED_firsttrial(anyisnan)=nan;
    if isempty([fliplr(0-binstep/2:-binstep:mi) 0+binstep/2:binstep:ma]) || isinf(mi) || isinf(ma) || isnan(mi) || isnan(ma) || range([fliplr(0-binstep/2:-binstep:mi) 0+binstep/2:binstep:ma])==0
        disp('No bins');
        return
    end
    [n,x]=hist(dprimes_noLED_lasttrial-dprimes_noLED_firsttrial,[fliplr(0-binstep/2:-binstep:mi) 0+binstep/2:binstep:ma]);
    [noLED_n,noLED_x]=cityscape_hist(n,x);
    [n,x]=hist(dprimes_LED_lasttrial-dprimes_LED_firsttrial,[fliplr(0-binstep/2:-binstep:mi) 0+binstep/2:binstep:ma]);
    [LED_n,LED_x]=cityscape_hist(n,x);
    % figure();
    % plot(noLED_x,noLED_n./nansum(noLED_n),'Color','k');
    % hold on;
    % plot(LED_x,LED_n./nansum(LED_n),'Color','r');
    figure();
    plot(noLED_x,noLED_n,'Color','k');
    hold on;
    plot(LED_x,LED_n,'Color','r'); title('Session by session change in dprime');
    if ~any(~isnan(dprimes_noLED_lasttrial-dprimes_noLED_firsttrial+dprimes_LED_lasttrial-dprimes_LED_firsttrial))
        disp('all nans');
    else
        pval=signrank(dprimes_noLED_lasttrial-dprimes_noLED_firsttrial,dprimes_LED_lasttrial-dprimes_LED_firsttrial);
        disp('pval from signrank comparing change in dprime');
        disp(pval);
    end
    figure();
    [black_x,black_y]=plotMouseByMouseChangeTrialToTrial(noLED_rr,metadata,'k',true); hold on;
    [red_x,red_y]=plotMouseByMouseChangeTrialToTrial(LED_rr,metadata,'r',true);
    figure();
    [black_x,black_y,changeindprime_bymouse_black]=plotMouseByMouseChangeTrialToTrial(noLED_rr,metadata,'k',false); hold on;
    [red_x,red_y,changeindprime_bymouse_red]=plotMouseByMouseChangeTrialToTrial(LED_rr,metadata,'r',false);
    [n,x]=hist(changeindprime_bymouse_black,[fliplr(0-binstep/2:-binstep:mi) 0+binstep/2:binstep:ma]);
    [noLED_n,noLED_x]=cityscape_hist(n,x);
    [n,x]=hist(changeindprime_bymouse_red,[fliplr(0-binstep/2:-binstep:mi) 0+binstep/2:binstep:ma]);
    [LED_n,LED_x]=cityscape_hist(n,x);
    figure();
    plot(noLED_x,noLED_n,'Color','k');
    hold on;
    plot(LED_x,LED_n,'Color','r'); title('Mouse by mouse change in dprime');
    if ~any(~isnan(changeindprime_bymouse_black+changeindprime_bymouse_red))
        disp('all nans');
    else
        pval=signrank(changeindprime_bymouse_black,changeindprime_bymouse_red);
        disp('pval from signrank comparing change in dprime MOUSE BY MOUSE');
        disp(pval);
    end
    returnout.changeindprime_bymouse_black=changeindprime_bymouse_black;
    returnout.changeindprime_bymouse_red=changeindprime_bymouse_red;
end

end

function [x_mousebymouse,y_mousebymouse,changeindprime_bymouse]=plotMouseByMouseChangeTrialToTrial(reachrates,metadata,c,doSessBySess)

atleastthismanytrials=3; 
[~,ui]=unique(metadata.sessid);
mouseids=metadata.mouseid(ui);

u=unique(mouseids);
mousebymouse_trialn_uncued_all=cell(1,length(u));
mousebymouse_trialn_cued_all=cell(1,length(u));
mousebymouse_trialnplus_uncued_all=cell(1,length(u));
mousebymouse_trialnplus_cued_all=cell(1,length(u));
for i=1:length(ui)
    currmouse=mouseids(i);
    
    temp=mousebymouse_trialn_uncued_all{u==currmouse};
    temprr=reachrates.trial1_alltrials_uncued(i,:);
    if isempty(temp)
        temp=temprr(~isnan(temprr));
    else
        temp=[temp temprr(~isnan(temprr))];
    end
    mousebymouse_trialn_uncued_all{u==currmouse}=temp;
    
    temp=mousebymouse_trialn_cued_all{u==currmouse};
    temprr=reachrates.trial1_alltrials_cued(i,:);
    if isempty(temp)
        temp=temprr(~isnan(temprr));
    else
        temp=[temp temprr(~isnan(temprr))];
    end
    mousebymouse_trialn_cued_all{u==currmouse}=temp;
    
    temp=mousebymouse_trialnplus_uncued_all{u==currmouse};
    temprr=reachrates.alltrials_uncued(i,:);
    if isempty(temp)
        temp=temprr(~isnan(temprr));
    else
        temp=[temp temprr(~isnan(temprr))];
    end
    mousebymouse_trialnplus_uncued_all{u==currmouse}=temp;
    
    temp=mousebymouse_trialnplus_cued_all{u==currmouse};
    temprr=reachrates.alltrials_cued(i,:);
    if isempty(temp)
        temp=temprr(~isnan(temprr));
    else
        temp=[temp temprr(~isnan(temprr))];
    end
    mousebymouse_trialnplus_cued_all{u==currmouse}=temp;
end
maxlength=0;
for i=1:length(mousebymouse_trialn_uncued_all)
    if length(mousebymouse_trialn_uncued_all{i})>maxlength
        maxlength=length(mousebymouse_trialn_uncued_all{i});
    end
end
for i=1:length(mousebymouse_trialn_uncued_all)
    if length(mousebymouse_trialn_uncued_all{i})<maxlength
        mousebymouse_trialn_uncued_all_asmatrix(i,:)=[mousebymouse_trialn_uncued_all{i} nan(1,maxlength-length(mousebymouse_trialn_uncued_all{i}))];
        mousebymouse_trialn_cued_all_asmatrix(i,:)=[mousebymouse_trialn_cued_all{i} nan(1,maxlength-length(mousebymouse_trialn_cued_all{i}))];
        mousebymouse_trialnplus_uncued_all_asmatrix(i,:)=[mousebymouse_trialnplus_uncued_all{i} nan(1,maxlength-length(mousebymouse_trialnplus_uncued_all{i}))];
        mousebymouse_trialnplus_cued_all_asmatrix(i,:)=[mousebymouse_trialnplus_cued_all{i} nan(1,maxlength-length(mousebymouse_trialnplus_cued_all{i}))];
    else
        mousebymouse_trialn_uncued_all_asmatrix(i,:)=mousebymouse_trialn_uncued_all{i};
        mousebymouse_trialn_cued_all_asmatrix(i,:)=mousebymouse_trialn_cued_all{i};
        mousebymouse_trialnplus_uncued_all_asmatrix(i,:)=mousebymouse_trialnplus_uncued_all{i};
        mousebymouse_trialnplus_cued_all_asmatrix(i,:)=mousebymouse_trialnplus_cued_all{i};
    end
end
dprimes_bymouse_firsttrial=calc_dprime_per_sess(mousebymouse_trialn_uncued_all_asmatrix,mousebymouse_trialn_cued_all_asmatrix);
dprimes_bymouse_trialnplus=calc_dprime_per_sess(mousebymouse_trialnplus_uncued_all_asmatrix,mousebymouse_trialnplus_cued_all_asmatrix);
changeindprime_bymouse=dprimes_bymouse_trialnplus-dprimes_bymouse_firsttrial;
mousebymouse_trialn_uncued=nan(1,length(u));
mousebymouse_trialn_cued=nan(1,length(u));
mousebymouse_trialnplus_uncued=nan(1,length(u));
mousebymouse_trialnplus_cued=nan(1,length(u));
for i=1:length(u)
    mousebymouse_trialn_uncued(i)=nanmean(mousebymouse_trialn_uncued_all{i});
    mousebymouse_trialn_cued(i)=nanmean(mousebymouse_trialn_cued_all{i});
    mousebymouse_trialnplus_uncued(i)=nanmean(mousebymouse_trialnplus_uncued_all{i});
    mousebymouse_trialnplus_cued(i)=nanmean(mousebymouse_trialnplus_cued_all{i});
end

temp=reachrates.trial1_alltrials_uncued;
temp(nansum(~isnan(temp),2)<atleastthismanytrials,:)=nan;
sessbysess_trialn_uncued=nanmean(temp,2);
temp=reachrates.trial1_alltrials_cued;
temp(nansum(~isnan(temp),2)<atleastthismanytrials,:)=nan;
sessbysess_trialn_cued=nanmean(temp,2);
temp=reachrates.alltrials_uncued;
temp(nansum(~isnan(temp),2)<atleastthismanytrials,:)=nan;
sessbysess_trialnplus_uncued=nanmean(temp,2);
temp=reachrates.alltrials_cued;
temp(nansum(~isnan(temp),2)<atleastthismanytrials,:)=nan;
sessbysess_trialnplus_cued=nanmean(temp,2);
% SESSION BY SESSION
temp1=sessbysess_trialnplus_uncued-sessbysess_trialn_uncued;
temp2=sessbysess_trialnplus_cued-sessbysess_trialn_cued;
a=[temp1(~isnan(temp1) & ~isnan(temp2)) temp2(~isnan(temp1) & ~isnan(temp2))];
if doSessBySess==true
    for i=1:size(a,1)
        quiver(0,0,a(i,1),a(i,2),'Color',c);
        hold on;
    end
    me=nanmean(a,1); s=nanstd(a,[],1)./sqrt(size(a,1));
    line([me(1) me(1)],[me(2)-s(2) me(2)+s(2)]); line([me(1)-s(1) me(1)+s(1)],[me(2) me(2)]);
    quiver(0,0,nanmean(a(:,1)),nanmean(a(:,2)),'Color',c,'LineWidth',4); title('session by session');
end

% mousebymouse_trialn_uncued=nan(1,length(unique(mouseids)));
% mousebymouse_trialn_cued=nan(1,length(unique(mouseids)));
% mousebymouse_trialnplus_uncued=nan(1,length(unique(mouseids)));
% mousebymouse_trialnplus_cued=nan(1,length(unique(mouseids)));
% u=unique(mouseids);
% for i=1:length(u)
%     currmouse=u(i);
%     mousebymouse_trialn_uncued(i)=nanmean(sessbysess_trialn_uncued(ismember(mouseids,currmouse)));
%     mousebymouse_trialn_cued(i)=nanmean(sessbysess_trialn_cued(ismember(mouseids,currmouse)));
%     mousebymouse_trialnplus_uncued(i)=nanmean(sessbysess_trialnplus_uncued(ismember(mouseids,currmouse)));
%     mousebymouse_trialnplus_cued(i)=nanmean(sessbysess_trialnplus_cued(ismember(mouseids,currmouse)));
% end

% MOUSE BY MOUSE
if doSessBySess==false
    for i=1:length(u)
        quiver(0,0,mousebymouse_trialnplus_uncued(i)-mousebymouse_trialn_uncued(i),mousebymouse_trialnplus_cued(i)-mousebymouse_trialn_cued(i),'Color',c);
        hold on;
    end
    me=nanmean([mousebymouse_trialnplus_uncued-mousebymouse_trialn_uncued mousebymouse_trialnplus_cued-mousebymouse_trialn_cued],1); 
    s=nanstd([mousebymouse_trialnplus_uncued-mousebymouse_trialn_uncued mousebymouse_trialnplus_cued-mousebymouse_trialn_cued],[],1)./sqrt(size([mousebymouse_trialnplus_uncued-mousebymouse_trialn_uncued mousebymouse_trialnplus_cued-mousebymouse_trialn_cued],1));
    line([me(1) me(1)],[me(2)-s(2) me(2)+s(2)]); line([me(1)-s(1) me(1)+s(1)],[me(2) me(2)]);
    quiver(0,0,nanmean(mousebymouse_trialnplus_uncued-mousebymouse_trialn_uncued),nanmean(mousebymouse_trialnplus_cued-mousebymouse_trialn_cued),'Color',c,'LineWidth',4);
    title('mouse by mouse');
end
x_mousebymouse=mousebymouse_trialnplus_uncued-mousebymouse_trialn_uncued;
y_mousebymouse=mousebymouse_trialnplus_cued-mousebymouse_trialn_cued;

end

function [dprimes_lasttrial,dprimes_firsttrial]=doPlottingAndBootstrap(reachrates,color1,color2,colordefault,thisplotset,compareToFirstTrial)

dprimes_lasttrial=[];
dprimes_firsttrial=[];
if isempty(reachrates)
    return
end
[uncued_mean_out,cued_mean_out]=bootstrap(reachrates.alltrials_uncued,reachrates.alltrials_cued,color1,color2,~thisplotset);
if thisplotset
    plotMeAndSe(reachrates.alltrials_uncued,reachrates.alltrials_cued,color1,2,false);
	plotMeAndSe(reachrates.alltrials_uncued,reachrates.alltrials_cued,color2,0.5,false);
	dprimes_lasttrial=calc_dprime_per_sess(reachrates.alltrials_uncued,reachrates.alltrials_cued);
end
if compareToFirstTrial==true
    [baseEffect_uncued_mean_out,baseEffect_cued_mean_out]=bootstrap(reachrates.trial1_alltrials_uncued,reachrates.trial1_alltrials_cued,colordefault,color2,~thisplotset);
    if thisplotset
        plotMeAndSe(reachrates.trial1_alltrials_uncued,reachrates.trial1_alltrials_cued,colordefault,2,false);
        plotMeAndSe(reachrates.trial1_alltrials_uncued,reachrates.trial1_alltrials_cued,color2,0.5,false);
        dprimes_firsttrial=calc_dprime_per_sess(reachrates.trial1_alltrials_uncued,reachrates.trial1_alltrials_cued);
    end
end
if thisplotset
    quiver(baseEffect_uncued_mean_out,baseEffect_cued_mean_out,uncued_mean_out-baseEffect_uncued_mean_out,cued_mean_out-baseEffect_cued_mean_out,'Color',color1,'LineWidth',2);
    quiver(baseEffect_uncued_mean_out,baseEffect_cued_mean_out,uncued_mean_out-baseEffect_uncued_mean_out,cued_mean_out-baseEffect_cued_mean_out,'Color',color2,'LineWidth',0.5);
% line([0 baseEffect_uncued_mean_out],[0 baseEffect_cued_mean_out],'Color',[0.2 0.2 0.2]);
end

end

function [plotset,trialTypes]=whichToPlotNow(whichToPlot,trialTypes,alltbt,timeWindow)

plotset.cuedtouchcolor=[171 104 87]./255;
switch whichToPlot
    case 'success'
        if isnan(timeWindow)
            timeWindow=[0 0.5]; % WHISPER mice
        end
        % for opto grc, used time window 0.25 to 1.25
%         timeWindow=[0.1 0.5]; %[0.1 1]; % change this to select only trials with a reach in this window
        plotset.wildcard=false;
        plotset.success=true;
        plotset.delayed=false;
        plotset.drop=false;
        plotset.cuedtouch=false;
        plotset.failedcued=false;
        plotset.falsealarm=false;
        plotset.noreach=false;
        plotset.backward_success=false;
    case 'delayed success'
        if isnan(timeWindow)
%         timeWindow=[5 9]; % from cue, in seconds
%         timeWindow=[3 7.5]; % from cue, in seconds
            timeWindow=[5 7.5]; % from cue, in seconds
        end
        plotset.wildcard=false;
        plotset.success=false;
        plotset.delayed=true;
        plotset.drop=false;
        plotset.cuedtouch=false;
        plotset.failedcued=false;
        plotset.falsealarm=false;
        plotset.noreach=false;
        plotset.backward_success=false;
    case 'backward success'
        if isnan(timeWindow)
            timeWindow=[];
        end
        plotset.wildcard=false;
        plotset.success=false;
        plotset.delayed=false;
        plotset.drop=false;
        plotset.cuedtouch=false;
        plotset.failedcued=false;
        plotset.falsealarm=false;
        plotset.noreach=false;
        plotset.backward_success=true;
    case 'drop'
        if isnan(timeWindow)
            timeWindow=[];
        end
        plotset.wildcard=false;
        plotset.success=false;
        plotset.delayed=false;
        plotset.drop=true;
        plotset.cuedtouch=false;
        plotset.failedcued=false;
        plotset.falsealarm=false;
        plotset.noreach=false;
        plotset.backward_success=false;
    case 'cued touch'
        if isnan(timeWindow)
            timeWindow=[];
        end
        plotset.wildcard=false;
        plotset.success=false;
        plotset.delayed=false;
        plotset.drop=false;
        plotset.cuedtouch=true;
        plotset.failedcued=false;
        plotset.falsealarm=false;
        plotset.noreach=false;
        plotset.backward_success=false;
    case 'cued touch and switch color'
        if isnan(timeWindow)
            timeWindow=[];
        end
        plotset.wildcard=false;
        plotset.success=false;
        plotset.delayed=false;
        plotset.drop=false;
        plotset.cuedtouch=true;
        plotset.failedcued=false;
        plotset.falsealarm=false;
        plotset.noreach=false;
        plotset.cuedtouchcolor='g';
        plotset.backward_success=false;
    case 'failed cued reach'
        if isnan(timeWindow)
            timeWindow=[]; %1.5];
        end
        plotset.wildcard=false;
        plotset.success=false;
        plotset.delayed=false;
        plotset.drop=false;
        plotset.cuedtouch=false;
        plotset.failedcued=true;
        plotset.falsealarm=false;
        plotset.noreach=false;
        plotset.backward_success=false;
    case 'false alarm'
        if isnan(timeWindow)
            timeWindow=[];
        end
        plotset.wildcard=false;
        plotset.success=false;
        plotset.delayed=false;
        plotset.drop=false;
        plotset.cuedtouch=false;
        plotset.failedcued=false;
        plotset.falsealarm=true;
        plotset.noreach=false;
        plotset.backward_success=false;
    case 'no reach'
        if isnan(timeWindow)
            timeWindow=[];
        end
        plotset.wildcard=false;
        plotset.success=false;
        plotset.delayed=false;
        plotset.drop=false;
        plotset.cuedtouch=false;
        plotset.failedcued=false;
        plotset.falsealarm=false;
        plotset.noreach=true;
        plotset.backward_success=false;
    case 'wildcard'
        if isnan(timeWindow)
            timeWindow=[];
        end
        plotset.wildcard=true;
        plotset.success=false;
        plotset.delayed=false;
        plotset.drop=false;
        plotset.cuedtouch=false;
        plotset.failedcued=false;
        plotset.falsealarm=false;
        plotset.noreach=false;
        plotset.backward_success=false;
    case 'basic'
        if isnan(timeWindow)
            timeWindow=[];
        end
        plotset.wildcard=false;
        plotset.success=true;
        plotset.delayed=false;
        plotset.drop=false;
        plotset.cuedtouch=true;
        plotset.failedcued=false;
        plotset.falsealarm=true;
        plotset.noreach=true;  
        plotset.backward_success=false;
end
[~,cueindma]=nanmax(nanmean(alltbt.cueZone_onVoff,1));
timestep=mode(diff(nanmean(alltbt.times,1)));
if isempty(timeWindow)
    trialTypes.reachedInTimeWindow=ones(size(trialTypes.led));
    trialTypes.reachedInTimeWindow_1forward=[trialTypes.reachedInTimeWindow(2:end); 1];
    trialTypes.reachedInTimeWindow_1back=[1; trialTypes.reachedInTimeWindow(1:end-1)];
else
    timeWindowInds(1)=floor(timeWindow(1)/timestep);
    timeWindowInds(2)=floor(timeWindow(2)/timestep);
    trialTypes.reachedInTimeWindow=any(alltbt.all_reachBatch(:,cueindma+timeWindowInds(1)-1:cueindma+timeWindowInds(2))>0.05,2);
    trialTypes.reachedInTimeWindow_1forward=[trialTypes.reachedInTimeWindow(2:end); 0];
    trialTypes.reachedInTimeWindow_1back=[0; trialTypes.reachedInTimeWindow(1:end-1)];
end

end

function dprimes=calc_dprime_per_sess(uncued_events,cued_events)

useBayes=false; % Bayes estimator helps to ameliorate SOME of the shift in d-prime that results from simply having too few trials

if useBayes==true
    disp('USING BAYES!');
    hit_rates=(nansum(cued_events>0,2)+1)./(nansum(~isnan(cued_events),2)+2);
    fa_rates=(nansum(uncued_events>0,2)+1)./(nansum(~isnan(uncued_events),2)+2);
else
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
end
dprimes=dprime(hit_rates,fa_rates);

end

function out=dprime(hit_rates,FA_rates)

out=norminv(hit_rates)-norminv(FA_rates);

end

function plotMeAndSe(data1,data2,c,linewidth,suppressOutput)
% make inputs vectors if they are not
data1=data1(1:end);
data2=data2(1:end);
% average within each session
% data1=nanmean(data1,2); data1=data1';
% data2=nanmean(data2,2); data2=data2';
if suppressOutput==false
    line([nanmean(data1)-nanstd(data1,[],2)./sqrt(nansum(~isnan(data1))) nanmean(data1)+nanstd(data1,[],2)./sqrt(nansum(~isnan(data1)))],[nanmean(data2) nanmean(data2)],'Color',c,'LineWidth',linewidth);
    hold on;
    line([nanmean(data1) nanmean(data1)],[nanmean(data2)-nanstd(data2,[],2)./sqrt(nansum(~isnan(data2))) nanmean(data2)+nanstd(data2,[],2)./sqrt(nansum(~isnan(data2)))],'Color',c,'LineWidth',linewidth);
end

end

function [test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir,fillInBetweenWithAnything)

tbt_filter.name='throwaway';
test.nInSequence=[nInSequence]; 
test.trial1=trial1;
if ~isempty(regexp(trial1,'SPLIT'))
    test.templateSequence2_cond=trial1;
else
    test.templateSequence2_cond=eval(trial1);
end
test.trial2=trial2;
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=fillInBetweenWithAnything; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name=['alltrials' tbt_filter.name 'inBetweenAnything' num2str(test.fillInBetweenWithAnything)];
saveDir2=[saveDir '\' test.event_name];
mkdir(saveDir2);
save([saveDir2 '\test_settings.mat'],'test');
fakeCueInd=50; 
skipCorrected=true;

end

function [uncued_mean_out,cued_mean_out,bootMeans]=bootstrap(approach2_alltrials_uncued,approach2_alltrials_cued,colorForBootstrapPoints,scatterPointsEdgeColor,suppressOutput)

stillsuppressbootstrap=false;

% altogether_prob_cued=nanmean(approach2_alltrials_cued,2);
% altogether_prob_uncued=nanmean(approach2_alltrials_uncued,2);
altogether_prob_cued=approach2_alltrials_cued(1:end); % better to bootstrap across trials, not sessions, because in some sessions, mouse drops a lot
altogether_prob_uncued=approach2_alltrials_uncued(1:end);
takeTrials=~isnan(altogether_prob_cued) & ~isnan(altogether_prob_uncued);
% nan trials are from cases where a session had fewer trials, just nanned
% to fill in matrix
% disp(['dropping this many trials because of nan ' num2str(nansum(~takeTrials))]);
altogether_prob_cued=altogether_prob_cued(takeTrials==1);
altogether_prob_uncued=altogether_prob_uncued(takeTrials==1);
% Show bootstrapped 95% CI
takeFracForBootstrap=0.66;
takeIndsForBootstrap=ceil(takeFracForBootstrap*length(altogether_prob_cued));
nRuns=100;
bootMeans=nan(2,nRuns);
for i=1:nRuns
    takeTheseForBoot=randi(length(altogether_prob_cued),1,takeIndsForBootstrap); % with replacement
    sub_prob_cued=altogether_prob_cued(takeTheseForBoot);
    sub_prob_uncued=altogether_prob_uncued(takeTheseForBoot);
    bootMeans(1,i)=nanmean(sub_prob_uncued);
    bootMeans(2,i)=nanmean(sub_prob_cued);
end
if suppressOutput==false && stillsuppressbootstrap==false
    s=scatter(bootMeans(1,:),bootMeans(2,:),20,'filled','MarkerEdgeColor',scatterPointsEdgeColor,'MarkerFaceColor',colorForBootstrapPoints); hold on;
    s.AlphaData = 0.5*ones(1,size(bootMeans,2));
    s.MarkerFaceAlpha = 'flat';
    scatter(nanmean(altogether_prob_uncued),nanmean(altogether_prob_cued),50,'filled','MarkerEdgeColor',scatterPointsEdgeColor,'MarkerFaceColor',colorForBootstrapPoints);
end
uncued_mean_out=nanmean(altogether_prob_uncued);
cued_mean_out=nanmean(altogether_prob_cued);

end