function outcomeDependentShift_wrapper(alltbt,trialTypes,metadata,saveDir)

compareToFirstTrial=true;
[~,cueindma]=nanmax(nanmean(alltbt.cueZone_onVoff,1));
timestep=mode(diff(nanmean(alltbt.times,1)));
% flankingTrials='trialTypes.optoGroup~=1 & trialTypes.chewing_at_trial_start==0';
flankingTrials='trialTypes.optoGroup~=1';
trialTypes.isLongITI_1forward=[trialTypes.isLongITI(2:end); 0];
trialTypes.optoGroup_1forward=[trialTypes.optoGroup(2:end); 0];
trialTypes.noReach=~any(alltbt.all_reachBatch>0.05,2);
trialTypes.noReach_1forward=[trialTypes.noReach(2:end); 0];
trialTypes.reachedBeforeCue=any(alltbt.all_reachBatch(:,1:cueindma-1)>0.05,2);
trialTypes.reachedBeforeCue_1forward=[trialTypes.reachedBeforeCue(2:end); 0];
timeWindow=[];
% timeWindow=[0.5 1]; % from cue, in seconds
% timeWindowInds(1)=floor(timeWindow(1)/timestep);
% timeWindowInds(2)=floor(timeWindow(2)/timestep);
if isempty(timeWindow)
    trialTypes.reachedInTimeWindow=ones(size(trialTypes.led));
else
    trialTypes.reachedInTimeWindow=any(alltbt.all_reachBatch(:,cueindma+timeWindowInds(1)-1:cueindma+timeWindowInds(2))>0.05,2);
end
trialTypes.reachedInTimeWindow_1forward=[trialTypes.reachedInTimeWindow(2:end); 0];

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
reachratesettings.addSatietyLines=true; % whether to add proportionality lines to figure
reachratesettings.useWindowsForUncued=[3]; % to use window2 or window3 or both for the uncued reach rate
reachratesettings.initWindows=[]; % empty if want to calculate from dataset
reachratesettings.addSessionLines=false; % for no averaging across sessions plot, whether to connect trial bins within same session with lines
reachratesettings.binThisManyTrials=25; % how many trials to bin within each session
reachratesettings.nBinsForZones=40; % will be nBinsForZones squared total bins, this is # bins for each x and y axis
reachratesettings.useRateMethod=3; % 1, 2 or 3 (see explanation below)
% see script_for_reaching_rate_analysis.m for explanation of rate methods

figure();
xlabel('Uncued reach rate (1/sec)');
% ylabel('Probability that reached faster after cue on second trial');
% title('Approach 2, window 2 vs window 1');
% baseline all trials
nInSequence=3;
trial1=[flankingTrials ' & trialTypes.led_1forward==0 & trialTypes.optoGroup_1forward~=1']; % & trialTypes.isLongITI_1forward==1'];
trial2=[flankingTrials];
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir);
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
if compareToFirstTrial==false
    [baseEffect_uncued_mean_out,baseEffect_cued_mean_out]=bootstrap(reachrates.alltrials_uncued,reachrates.alltrials_cued,'k','k');
end
% success
nInSequence=3;
trial1=[flankingTrials ' & trialTypes.consumed_pellet_1back==1' ' & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.after_cue_success_1forward==1 & trialTypes.consumed_pellet_1forward==1 & trialTypes.led_1forward==0 & trialTypes.optoGroup_1forward~=1']; % & trialTypes.isLongITI_1forward==1'];
trial2=[flankingTrials];
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir);
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
[uncued_mean_out,cued_mean_out]=bootstrap(reachrates.alltrials_uncued,reachrates.alltrials_cued,'g','g');
if compareToFirstTrial==true
    [baseEffect_uncued_mean_out,baseEffect_cued_mean_out]=bootstrap(reachrates.trial1_alltrials_uncued,reachrates.trial1_alltrials_cued,'k','k');
end
quiver(baseEffect_uncued_mean_out,baseEffect_cued_mean_out,uncued_mean_out-baseEffect_uncued_mean_out,cued_mean_out-baseEffect_cued_mean_out,'Color','g','LineWidth',2);
line([0 baseEffect_uncued_mean_out],[0 baseEffect_cued_mean_out],'Color',[0.2 0.2 0.2]);
% drop
nInSequence=3;
trial1=[flankingTrials ' & trialTypes.consumed_pellet_1back==0' ' & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.after_cue_drop_1forward==1 & trialTypes.consumed_pellet_1forward==0 & trialTypes.led_1forward==0 & trialTypes.optoGroup_1forward~=1']; % & trialTypes.isLongITI_1forward==1'];
trial2=[flankingTrials];
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir);
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
[uncued_mean_out,cued_mean_out]=bootstrap(reachrates.alltrials_uncued,reachrates.alltrials_cued,'r','r');
if compareToFirstTrial==true
    [baseEffect_uncued_mean_out,baseEffect_cued_mean_out]=bootstrap(reachrates.trial1_alltrials_uncued,reachrates.trial1_alltrials_cued,'k','k');
end
quiver(baseEffect_uncued_mean_out,baseEffect_cued_mean_out,uncued_mean_out-baseEffect_uncued_mean_out,cued_mean_out-baseEffect_cued_mean_out,'Color','r','LineWidth',2);
% touched after cue, i.e., success or drop
nInSequence=3;
trial1=[flankingTrials ' & trialTypes.touched_pellet_1back==1' ' & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touch_in_cued_window_1forward==1 & trialTypes.led_1forward==0 & trialTypes.optoGroup_1forward~=1']; % & trialTypes.isLongITI_1forward==1'];
trial2=[flankingTrials];
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir);
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
[uncued_mean_out,cued_mean_out]=bootstrap(reachrates.alltrials_uncued,reachrates.alltrials_cued,[171 104 87]./255,[171 104 87]./255);
if compareToFirstTrial==true
    [baseEffect_uncued_mean_out,baseEffect_cued_mean_out]=bootstrap(reachrates.trial1_alltrials_uncued,reachrates.trial1_alltrials_cued,'k','k');
end
quiver(baseEffect_uncued_mean_out,baseEffect_cued_mean_out,uncued_mean_out-baseEffect_uncued_mean_out,cued_mean_out-baseEffect_cued_mean_out,'Color',[171 104 87]./255,'LineWidth',2);
% miss before cue
nInSequence=3;
trial1=[flankingTrials ' & trialTypes.reachedBeforeCue_1forward==1 & trialTypes.led_1forward==0 & trialTypes.optoGroup_1forward~=1']; % & trialTypes.isLongITI_1forward==1'];
trial2=[flankingTrials];
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir);
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
[uncued_mean_out,cued_mean_out]=bootstrap(reachrates.alltrials_uncued,reachrates.alltrials_cued,'c','c');
if compareToFirstTrial==true
    [baseEffect_uncued_mean_out,baseEffect_cued_mean_out]=bootstrap(reachrates.trial1_alltrials_uncued,reachrates.trial1_alltrials_cued,'k','k');
end
quiver(baseEffect_uncued_mean_out,baseEffect_cued_mean_out,uncued_mean_out-baseEffect_uncued_mean_out,cued_mean_out-baseEffect_cued_mean_out,'Color','c','LineWidth',2);
% does not reach
nInSequence=3;
trial1=[flankingTrials ' & trialTypes.noReach_1forward==1 & trialTypes.touched_pellet_1forward==0 & trialTypes.led_1forward==0 & trialTypes.optoGroup_1forward~=1']; % & trialTypes.isLongITI_1forward==1'];
trial2=[flankingTrials];
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir);
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
[uncued_mean_out,cued_mean_out]=bootstrap(reachrates.alltrials_uncued,reachrates.alltrials_cued,[0.8 0.8 0.8],[0.8 0.8 0.8]);
if compareToFirstTrial==true
    [baseEffect_uncued_mean_out,baseEffect_cued_mean_out]=bootstrap(reachrates.trial1_alltrials_uncued,reachrates.trial1_alltrials_cued,'k','k');
end
quiver(baseEffect_uncued_mean_out,baseEffect_cued_mean_out,uncued_mean_out-baseEffect_uncued_mean_out,cued_mean_out-baseEffect_cued_mean_out,'Color',[0.8 0.8 0.8],'LineWidth',2);

%%%% LEDs
% baseline all trials
nInSequence=3;
trial1=[flankingTrials ' & trialTypes.led_1forward==1 & trialTypes.optoGroup_1forward~=1']; % & trialTypes.isLongITI_1forward==1'];
trial2=[flankingTrials];
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir);
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
if compareToFirstTrial==false
    [baseEffect_uncued_mean_out,baseEffect_cued_mean_out]=bootstrap(reachrates.alltrials_uncued,reachrates.alltrials_cued,'k','m');
end
% success
nInSequence=3;
trial1=[flankingTrials ' & trialTypes.consumed_pellet_1back==1' ' & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.after_cue_success_1forward==1 & trialTypes.consumed_pellet_1forward==1 & trialTypes.led_1forward==1 & trialTypes.optoGroup_1forward~=1']; % & trialTypes.isLongITI_1forward==1'];
trial2=[flankingTrials];
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir);
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
[uncued_mean_out,cued_mean_out]=bootstrap(reachrates.alltrials_uncued,reachrates.alltrials_cued,'g','m');
if compareToFirstTrial==true
    [baseEffect_uncued_mean_out,baseEffect_cued_mean_out]=bootstrap(reachrates.trial1_alltrials_uncued,reachrates.trial1_alltrials_cued,'k','m');
end
quiver(baseEffect_uncued_mean_out,baseEffect_cued_mean_out,uncued_mean_out-baseEffect_uncued_mean_out,cued_mean_out-baseEffect_cued_mean_out,'Color','g','LineWidth',2);
quiver(baseEffect_uncued_mean_out,baseEffect_cued_mean_out,uncued_mean_out-baseEffect_uncued_mean_out,cued_mean_out-baseEffect_cued_mean_out,'Color','m','LineWidth',0.5);
line([0 baseEffect_uncued_mean_out],[0 baseEffect_cued_mean_out],'Color',[0.2 0.2 0.2]);
% drop
nInSequence=3;
trial1=[flankingTrials ' & trialTypes.consumed_pellet_1back==0' ' & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.after_cue_drop_1forward==1 & trialTypes.consumed_pellet_1forward==0 & trialTypes.led_1forward==1 & trialTypes.optoGroup_1forward~=1']; % & trialTypes.isLongITI_1forward==1'];
trial2=[flankingTrials];
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir);
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
[uncued_mean_out,cued_mean_out]=bootstrap(reachrates.alltrials_uncued,reachrates.alltrials_cued,'r','m');
if compareToFirstTrial==true
    [baseEffect_uncued_mean_out,baseEffect_cued_mean_out]=bootstrap(reachrates.trial1_alltrials_uncued,reachrates.trial1_alltrials_cued,'k','m');
end
quiver(baseEffect_uncued_mean_out,baseEffect_cued_mean_out,uncued_mean_out-baseEffect_uncued_mean_out,cued_mean_out-baseEffect_cued_mean_out,'Color','r','LineWidth',2);
quiver(baseEffect_uncued_mean_out,baseEffect_cued_mean_out,uncued_mean_out-baseEffect_uncued_mean_out,cued_mean_out-baseEffect_cued_mean_out,'Color','m','LineWidth',0.5);
% touched after cue, i.e., success or drop
nInSequence=3;
trial1=[flankingTrials ' & trialTypes.touched_pellet_1back==1' ' & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.touch_in_cued_window_1forward==1 & trialTypes.led_1forward==1 & trialTypes.optoGroup_1forward~=1']; % & trialTypes.isLongITI_1forward==1'];
trial2=[flankingTrials];
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir);
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
[uncued_mean_out,cued_mean_out]=bootstrap(reachrates.alltrials_uncued,reachrates.alltrials_cued,[171 104 87]./255,'m');
if compareToFirstTrial==true
    [baseEffect_uncued_mean_out,baseEffect_cued_mean_out]=bootstrap(reachrates.trial1_alltrials_uncued,reachrates.trial1_alltrials_cued,'k','m');
end
quiver(baseEffect_uncued_mean_out,baseEffect_cued_mean_out,uncued_mean_out-baseEffect_uncued_mean_out,cued_mean_out-baseEffect_cued_mean_out,'Color',[171 104 87]./255,'LineWidth',2);
quiver(baseEffect_uncued_mean_out,baseEffect_cued_mean_out,uncued_mean_out-baseEffect_uncued_mean_out,cued_mean_out-baseEffect_cued_mean_out,'Color','m','LineWidth',0.5);
% miss before cue
nInSequence=3;
trial1=[flankingTrials ' & trialTypes.reachedBeforeCue_1forward==1 & trialTypes.led_1forward==1 & trialTypes.optoGroup_1forward~=1']; % & trialTypes.isLongITI_1forward==1'];
trial2=[flankingTrials];
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir);
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
[uncued_mean_out,cued_mean_out]=bootstrap(reachrates.alltrials_uncued,reachrates.alltrials_cued,'c','m');
if compareToFirstTrial==true
    [baseEffect_uncued_mean_out,baseEffect_cued_mean_out]=bootstrap(reachrates.trial1_alltrials_uncued,reachrates.trial1_alltrials_cued,'k','m');
end
quiver(baseEffect_uncued_mean_out,baseEffect_cued_mean_out,uncued_mean_out-baseEffect_uncued_mean_out,cued_mean_out-baseEffect_cued_mean_out,'Color','c','LineWidth',2);
quiver(baseEffect_uncued_mean_out,baseEffect_cued_mean_out,uncued_mean_out-baseEffect_uncued_mean_out,cued_mean_out-baseEffect_cued_mean_out,'Color','m','LineWidth',0.5);
% does not reach
nInSequence=3;
trial1=[flankingTrials ' & trialTypes.noReach_1forward==1 & trialTypes.touched_pellet_1forward==0 & trialTypes.led_1forward==1 & trialTypes.optoGroup_1forward~=1']; % & trialTypes.isLongITI_1forward==1'];
trial2=[flankingTrials];
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir);
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
[uncued_mean_out,cued_mean_out]=bootstrap(reachrates.alltrials_uncued,reachrates.alltrials_cued,[0.8 0.8 0.8],'m');
if compareToFirstTrial==true
    [baseEffect_uncued_mean_out,baseEffect_cued_mean_out]=bootstrap(reachrates.trial1_alltrials_uncued,reachrates.trial1_alltrials_cued,'k','m');
end
quiver(baseEffect_uncued_mean_out,baseEffect_cued_mean_out,uncued_mean_out-baseEffect_uncued_mean_out,cued_mean_out-baseEffect_cued_mean_out,'Color',[0.8 0.8 0.8],'LineWidth',2);
quiver(baseEffect_uncued_mean_out,baseEffect_cued_mean_out,uncued_mean_out-baseEffect_uncued_mean_out,cued_mean_out-baseEffect_cued_mean_out,'Color','m','LineWidth',0.5);

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

function [uncued_mean_out,cued_mean_out,bootMeans]=bootstrap(approach2_alltrials_uncued,approach2_alltrials_cued,colorForBootstrapPoints,scatterPointsEdgeColor)

altogether_prob_cued=nanmean(approach2_alltrials_cued,2);
altogether_prob_uncued=nanmean(approach2_alltrials_uncued,2);
% altogether_prob_cued=approach2_alltrials_cued(1:end); % better to bootstrap across trials, not sessions, because in some sessions, mouse drops a lot
% altogether_prob_uncued=approach2_alltrials_uncued(1:end);
takeTrials=~isnan(altogether_prob_cued) & ~isnan(altogether_prob_uncued);
disp(['dropping this many trials because of nan ' num2str(nansum(~takeTrials))]);
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
s=scatter(bootMeans(1,:),bootMeans(2,:),20,'filled','MarkerEdgeColor',scatterPointsEdgeColor,'MarkerFaceColor',colorForBootstrapPoints); hold on;
s.AlphaData = 0.5*ones(1,size(bootMeans,2));
s.MarkerFaceAlpha = 'flat';
scatter(nanmean(altogether_prob_uncued),nanmean(altogether_prob_cued),50,'filled','MarkerEdgeColor',scatterPointsEdgeColor,'MarkerFaceColor',colorForBootstrapPoints);
uncued_mean_out=nanmean(altogether_prob_uncued);
cued_mean_out=nanmean(altogether_prob_cued);

end