function plotSummaryOfRTresults(alltbt,trialTypes,metadata,ledIsOn,useLongITIonly,excludeFromAllTrials,alsoExcludeFromSubset,useCorrected)

cueName='cueZone_onVoff';

% useCorrected==true if want to use corrected input distributions

% Set up figure layout
figure();
Nh=5; % number of rows
Nw=3; % number of columns
gap=[.01 .03]; % between plots
marg_h=[.1 .01]; % margin
marg_w=[.01 .01]; % margin
[ha,pos]=tight_subplot(Nh,Nw,gap,marg_h,marg_w);

if ~isempty(ledIsOn)
    trialTypes.ledIsOn=ledIsOn;
else
    trialTypes.ledIsOn=trialTypes.led;
end
if useLongITIonly==true
    trialTypes.ITItoUse=trialTypes.isLongITI==1;
else
    trialTypes.ITItoUse=ones(size(trialTypes.led));
end

if ~isempty(excludeFromAllTrials)
    % exclude these trials from "all"
    allTrialsToUse=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
    trialTypes.allTrialsToUse=allTrialsToUse & ~(excludeFromAllTrials==1);
else
    trialTypes.allTrialsToUse=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
end
if alsoExcludeFromSubset==true
    trialTypes.allTrialsForSubset=trialTypes.allTrialsToUse;
else
    trialTypes.allTrialsForSubset=ones(size(trialTypes.allTrialsToUse));
end

% 1. Hypothesis: touched pellet during cued reach on trial N makes reaction time (RT) on trial N+1 faster
% -------------- shift in all trials RT distribution
test.nInSequence=[2]; 
trial1='trialTypes.touch_in_cued_window==1 & trialTypes.led==0 & trialTypes.ITItoUse==1';
test.trial1=trial1;
test.trial1_name='touched pellet during cued reach no opto';
test.templateSequence2_cond=eval(trial1);
trial2='trialTypes.allTrialsToUse==1';
test.trial2=trial2;
test.trial2_name='any';
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=false; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name='temp';
% build paired RT data set
fakeCueInd=50; % in indices
if useCorrected==false
    dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,true);
    plotReactionTimeDotScatter(alltbt,'cueZone_onVoff','all_reachBatch','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts',trialTypes,metadata,true,dataset.realDistributions.event_isSeq{1},true,true);
    title(['First trial: ' test.trial1_name]);
    ylabel(['Second trial: ' test.trial2_name]);
    [event_cdf,alltrials_cdf]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rt_cdf',true);
else
    [~,dataset]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,false);
    [event_cdf,alltrials_cdf]=plotBehaviorEventFx(dataset,alltbt,[],'plot_rt_cdf',true);
end
axes(ha(1));
plot(alltrials_cdf.x,alltrials_cdf.y,'Color','k','LineWidth',1);
hold on;
plot(event_cdf.x,event_cdf.y,'Color',[1 0.6 0.6],'LineWidth',1); 
% opto
test.nInSequence=[2]; 
trial1='trialTypes.touch_in_cued_window==1 & trialTypes.ledIsOn==1 & trialTypes.ITItoUse==1';
test.trial1=trial1;
test.trial1_name='touched pellet during cued reach w opto';
test.templateSequence2_cond=eval(trial1);
trial2='trialTypes.allTrialsToUse==1';
test.trial2=trial2;
test.trial2_name='any';
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=false; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name='temp';
% build paired RT data set
fakeCueInd=50; % in indices
if useCorrected==false
    dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,true); 
    plotReactionTimeDotScatter(alltbt,'cueZone_onVoff','all_reachBatch','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts',trialTypes,metadata,true,dataset.realDistributions.event_isSeq{1},true,true);
    title(['First trial: ' test.trial1_name]);
    ylabel(['Second trial: ' test.trial2_name]);
    [event_cdf,~]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rt_cdf',true);
else
    [~,dataset]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,false);
    [event_cdf,~]=plotBehaviorEventFx(dataset,alltbt,[],'plot_rt_cdf',true);
end
axes(ha(1));
plot(event_cdf.x,event_cdf.y,'Color','r','LineWidth',1); 
h=legend({'All trials','All after cued reach touch','All aft cue reach touch opto'},'Location','best');
legend('boxoff');
set(h,'FontSize',6); 
xlabel('RT (sec)');
ylabel('CDF');
% -------------- shift in matched trial type distribution
test.nInSequence=[2]; 
trial1='trialTypes.touch_in_cued_window==1 & trialTypes.led==0 & trialTypes.ITItoUse==1';
test.trial1=trial1;
test.trial1_name='touched pellet during cued reach no opto';
test.templateSequence2_cond=eval(trial1);
trial2='trialTypes.touch_in_cued_window==1 & trialTypes.ITItoUse==1 & trialTypes.allTrialsForSubset==1';
test.trial2=trial2;
test.trial2_name='touched pellet during cued reach';
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=false; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name='temp';
% build paired RT data set
fakeCueInd=50; % in indices
if useCorrected==false
    dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,true);
    plotReactionTimeDotScatter(alltbt,'cueZone_onVoff','all_reachBatch','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts',trialTypes,metadata,true,dataset.realDistributions.event_isSeq{1},true,true);
    title(['First trial: ' test.trial1_name]);
    ylabel(['Second trial: ' test.trial2_name]);
    [event_cdf,alltrials_cdf]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rt_cdf',true);
    [delta_event_cdf,delta_alltrials_cdf]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_delta_rt_cdf',true);
else
    [~,dataset]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,false);
    [event_cdf,alltrials_cdf]=plotBehaviorEventFx(dataset,alltbt,[],'plot_rt_cdf',true);
    [delta_event_cdf,delta_alltrials_cdf]=plotBehaviorEventFx(dataset,alltbt,[],'plot_delta_rt_cdf',true);
end
axes(ha(2));
plot(alltrials_cdf.x,alltrials_cdf.y,'Color','k','LineWidth',1);
hold on;
plot(event_cdf.x,event_cdf.y,'Color',[1 0.6 0.6],'LineWidth',1); 
% opto
test.nInSequence=[2]; 
trial1='trialTypes.touch_in_cued_window==1 & trialTypes.ledIsOn==1 & trialTypes.ITItoUse==1';
test.trial1=trial1;
test.trial1_name='touched pellet during cued reach w opto';
test.templateSequence2_cond=eval(trial1);
trial2='trialTypes.touch_in_cued_window==1 & trialTypes.ITItoUse==1 & trialTypes.allTrialsForSubset==1';
test.trial2=trial2;
test.trial2_name='touched pellet during cued reach';
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=false; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name='temp';
% build paired RT data set
fakeCueInd=50; % in indices
if useCorrected==false
    dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,true);
    plotReactionTimeDotScatter(alltbt,'cueZone_onVoff','all_reachBatch','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts',trialTypes,metadata,true,dataset.realDistributions.event_isSeq{1},true,true);
    title(['First trial: ' test.trial1_name]);
    ylabel(['Second trial: ' test.trial2_name]);
    [event_cdf,~]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rt_cdf',true);
    [deltaopto_event_cdf,~]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_delta_rt_cdf',true);
else
    [~,dataset]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,false);
    [event_cdf,~]=plotBehaviorEventFx(dataset,alltbt,[],'plot_rt_cdf',true);
    [deltaopto_event_cdf,~]=plotBehaviorEventFx(dataset,alltbt,[],'plot_delta_rt_cdf',true);
end
axes(ha(2));
plot(event_cdf.x,event_cdf.y,'Color','r','LineWidth',1); 
h=legend({'All trials','Cued reach touch after cued reach touch','Cued reach touch aft cue reach touch opto'},'Location','best');
legend('boxoff');
set(h,'FontSize',6); 
xlabel('RT (sec)');
ylabel('CDF');
% change in rt
axes(ha(3));
plot(delta_alltrials_cdf.x,delta_alltrials_cdf.y,'Color','k','LineWidth',1);
hold on;
plot(delta_event_cdf.x,delta_event_cdf.y,'Color',[1 0.6 0.6],'LineWidth',1);
plot(deltaopto_event_cdf.x,deltaopto_event_cdf.y,'Color','r','LineWidth',1); 
h=legend({'All trials','Cued reach touch after cued reach touch','Cued reach touch aft cue reach touch opto'},'Location','best');
legend('boxoff');
set(h,'FontSize',6); 
xlabel('Change in RT (sec)');
ylabel('CDF');



% 2. Hypothesis: touched pellet after cue on trial N makes reaction time (RT) on trial N+1 faster
% -------------- shift in all trials RT distribution
test.nInSequence=[2]; 
trial1='(trialTypes.after_cue_drop==1 | trialTypes.after_cue_success==1) & trialTypes.led==0 & trialTypes.ITItoUse==1';
test.trial1=trial1;
test.trial1_name='touched pellet after cue no opto';
test.templateSequence2_cond=eval(trial1);
trial2='trialTypes.allTrialsToUse==1';
test.trial2=trial2;
test.trial2_name='any';
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=false; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name='temp';
% build paired RT data set
fakeCueInd=50; % in indices
if useCorrected==false
    dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,true);
    plotReactionTimeDotScatter(alltbt,'cueZone_onVoff','all_reachBatch','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts',trialTypes,metadata,true,dataset.realDistributions.event_isSeq{1},true,true);
    title(['First trial: ' test.trial1_name]);
    ylabel(['Second trial: ' test.trial2_name]);
    [event_cdf,alltrials_cdf]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rt_cdf',true);
else
    [~,dataset]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,false);
    [event_cdf,alltrials_cdf]=plotBehaviorEventFx(dataset,alltbt,[],'plot_rt_cdf',true);
end
axes(ha(4));
plot(alltrials_cdf.x,alltrials_cdf.y,'Color','k','LineWidth',1);
hold on;
plot(event_cdf.x,event_cdf.y,'Color',[1 0.6 0.6],'LineWidth',1); 
allCuedTouch=event_cdf;
% opto
test.nInSequence=[2]; 
trial1='(trialTypes.after_cue_drop==1 | trialTypes.after_cue_success==1) & trialTypes.ledIsOn==1 & trialTypes.ITItoUse==1';
test.trial1=trial1;
test.trial1_name='touched pellet after cue w opto';
test.templateSequence2_cond=eval(trial1);
trial2='trialTypes.allTrialsToUse==1';
test.trial2=trial2;
test.trial2_name='any';
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=false; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name='temp';
% build paired RT data set
fakeCueInd=50; % in indices
if useCorrected==false
    dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,true);
    plotReactionTimeDotScatter(alltbt,'cueZone_onVoff','all_reachBatch','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts',trialTypes,metadata,true,dataset.realDistributions.event_isSeq{1},true,true);
    title(['First trial: ' test.trial1_name]);
    ylabel(['Second trial: ' test.trial2_name]);
    [event_cdf,~]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rt_cdf',true);
else
    [~,dataset]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,false);
    [event_cdf,~]=plotBehaviorEventFx(dataset,alltbt,[],'plot_rt_cdf',true);
end
axes(ha(4));
plot(event_cdf.x,event_cdf.y,'Color','r','LineWidth',1); 
h=legend({'All trials','All after touch after cue','All aft touch after cue opto'},'Location','best');
legend('boxoff');
set(h,'FontSize',6); 
xlabel('RT (sec)');
ylabel('CDF');
% -------------- shift in matched trial type distribution
test.nInSequence=[2]; 
trial1='(trialTypes.after_cue_drop==1 | trialTypes.after_cue_success==1) & trialTypes.led==0 & trialTypes.ITItoUse==1';
test.trial1=trial1;
test.trial1_name='touched after cue no opto';
test.templateSequence2_cond=eval(trial1);
trial2='(trialTypes.after_cue_drop==1 | trialTypes.after_cue_success==1) & trialTypes.ITItoUse==1 & trialTypes.allTrialsForSubset==1';
test.trial2=trial2;
test.trial2_name='touched after cue';
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=false; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name='temp';
% build paired RT data set
fakeCueInd=50; % in indices
if useCorrected==false
    dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,true);
    plotReactionTimeDotScatter(alltbt,'cueZone_onVoff','all_reachBatch','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts',trialTypes,metadata,true,dataset.realDistributions.event_isSeq{1},true,true);
    title(['First trial: ' test.trial1_name]);
    ylabel(['Second trial: ' test.trial2_name]);
    [event_cdf,alltrials_cdf]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rt_cdf',true);
    [delta_event_cdf,delta_alltrials_cdf]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_delta_rt_cdf',true);
else
    [~,dataset]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,false);
    [event_cdf,alltrials_cdf]=plotBehaviorEventFx(dataset,alltbt,[],'plot_rt_cdf',true);
    [delta_event_cdf,delta_alltrials_cdf]=plotBehaviorEventFx(dataset,alltbt,[],'plot_delta_rt_cdf',true);
end
axes(ha(5));
plot(alltrials_cdf.x,alltrials_cdf.y,'Color','k','LineWidth',1);
hold on;
plot(event_cdf.x,event_cdf.y,'Color',[1 0.6 0.6],'LineWidth',1); 
matchedCuedTouch=event_cdf;
% opto
test.nInSequence=[2]; 
trial1='(trialTypes.after_cue_drop==1 | trialTypes.after_cue_success==1) & trialTypes.ledIsOn==1 & trialTypes.ITItoUse==1';
test.trial1=trial1;
test.trial1_name='touch after cue w opto';
test.templateSequence2_cond=eval(trial1);
trial2='(trialTypes.after_cue_drop==1 | trialTypes.after_cue_success==1) & trialTypes.ITItoUse==1 & trialTypes.allTrialsForSubset==1';
test.trial2=trial2;
test.trial2_name='touch after cue';
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=false; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name='temp';
% build paired RT data set
fakeCueInd=50; % in indices
if useCorrected==false
    dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,true);
    plotReactionTimeDotScatter(alltbt,'cueZone_onVoff','all_reachBatch','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts',trialTypes,metadata,true,dataset.realDistributions.event_isSeq{1},true,true);
    title(['First trial: ' test.trial1_name]);
    ylabel(['Second trial: ' test.trial2_name]);
    [event_cdf,~]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rt_cdf',true);
    [deltaopto_event_cdf,~]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_delta_rt_cdf',true);
else
    [~,dataset]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,false);
    [event_cdf,~]=plotBehaviorEventFx(dataset,alltbt,[],'plot_rt_cdf',true);
    [deltaopto_event_cdf,~]=plotBehaviorEventFx(dataset,alltbt,[],'plot_delta_rt_cdf',true);
end
axes(ha(5));
plot(event_cdf.x,event_cdf.y,'Color','r','LineWidth',1); 
h=legend({'All trials','Touch after cue after touch after cue','Touch after cue after touch after cue opto'},'Location','best');
legend('boxoff');
set(h,'FontSize',6); 
xlabel('RT (sec)');
ylabel('CDF');
% change in rt
axes(ha(6));
plot(delta_alltrials_cdf.x,delta_alltrials_cdf.y,'Color','k','LineWidth',1);
hold on; 
plot(delta_event_cdf.x,delta_event_cdf.y,'Color',[1 0.6 0.6],'LineWidth',1); 
deltaCuedTouch=delta_event_cdf;
plot(deltaopto_event_cdf.x,deltaopto_event_cdf.y,'Color','r','LineWidth',1); 
h=legend({'All trials','Touch after cue after touch after cue','Touch after cue aft touch after cue opto'},'Location','best');
legend('boxoff');
set(h,'FontSize',6); 
xlabel('Change in RT (sec)');
ylabel('CDF');



% 3. Hypothesis: successful reach after cue on trial N makes reaction time (RT) on trial N+1 faster
% -------------- shift in all trials RT distribution
test.nInSequence=[2]; 
trial1='(trialTypes.after_cue_success==1) & trialTypes.led==0 & trialTypes.ITItoUse==1';
test.trial1=trial1;
test.trial1_name='ate pellet after cue no opto';
test.templateSequence2_cond=eval(trial1);
trial2='trialTypes.allTrialsToUse==1';
test.trial2=trial2;
test.trial2_name='any';
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=false; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name='temp';
% build paired RT data set
fakeCueInd=50; % in indices
if useCorrected==false
    dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,true);
    plotReactionTimeDotScatter(alltbt,'cueZone_onVoff','all_reachBatch','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts',trialTypes,metadata,true,dataset.realDistributions.event_isSeq{1},true,true);
    title(['First trial: ' test.trial1_name]);
    ylabel(['Second trial: ' test.trial2_name]);
    [event_cdf,alltrials_cdf]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rt_cdf',true);
else
    [~,dataset]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,false);
    [event_cdf,alltrials_cdf]=plotBehaviorEventFx(dataset,alltbt,[],'plot_rt_cdf',true);
end
axes(ha(7));
plot(alltrials_cdf.x,alltrials_cdf.y,'Color','k','LineWidth',1);
hold on;
plot(event_cdf.x,event_cdf.y,'Color',[0.6 1 0.6],'LineWidth',1); 
% opto
test.nInSequence=[2]; 
trial1='(trialTypes.after_cue_success==1) & trialTypes.ledIsOn==1 & trialTypes.ITItoUse==1';
test.trial1=trial1;
test.trial1_name='ate pellet after cue w opto';
test.templateSequence2_cond=eval(trial1);
trial2='trialTypes.allTrialsToUse==1';
test.trial2=trial2;
test.trial2_name='any';
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=false; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name='temp';
% build paired RT data set
fakeCueInd=50; % in indices
if useCorrected==false
    dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,true);
    plotReactionTimeDotScatter(alltbt,'cueZone_onVoff','all_reachBatch','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts',trialTypes,metadata,true,dataset.realDistributions.event_isSeq{1},true,true);
    title(['First trial: ' test.trial1_name]);
    ylabel(['Second trial: ' test.trial2_name]);
    [event_cdf,~]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rt_cdf',true);
else
    [~,dataset]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,false);
    [event_cdf,~]=plotBehaviorEventFx(dataset,alltbt,[],'plot_rt_cdf',true);
end
axes(ha(7));
plot(event_cdf.x,event_cdf.y,'Color','r','LineWidth',1); 
h=legend({'All trials','All after ate after cue','All aft ate after cue opto'},'Location','best');
legend('boxoff');
set(h,'FontSize',6); 
xlabel('RT (sec)');
ylabel('CDF');
% -------------- shift in matched trial type distribution
test.nInSequence=[2]; 
trial1='(trialTypes.after_cue_success==1) & trialTypes.led==0 & trialTypes.ITItoUse==1';
test.trial1=trial1;
test.trial1_name='ate after cue no opto';
test.templateSequence2_cond=eval(trial1);
trial2='(trialTypes.after_cue_success==1) & trialTypes.ITItoUse==1 & trialTypes.allTrialsForSubset==1';
test.trial2=trial2;
test.trial2_name='ate after cue';
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=false; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name='temp';
% build paired RT data set
fakeCueInd=50; % in indices
if useCorrected==false
    dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,true);
    plotReactionTimeDotScatter(alltbt,'cueZone_onVoff','all_reachBatch','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts',trialTypes,metadata,true,dataset.realDistributions.event_isSeq{1},true,true);
    title(['First trial: ' test.trial1_name]);
    ylabel(['Second trial: ' test.trial2_name]);
    [event_cdf,alltrials_cdf]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rt_cdf',true);
    [delta_event_cdf,delta_alltrials_cdf]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_delta_rt_cdf',true);
else
    [~,dataset]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,false);
    [event_cdf,alltrials_cdf]=plotBehaviorEventFx(dataset,alltbt,[],'plot_rt_cdf',true);
    [delta_event_cdf,delta_alltrials_cdf]=plotBehaviorEventFx(dataset,alltbt,[],'plot_delta_rt_cdf',true);
end
axes(ha(8));
plot(alltrials_cdf.x,alltrials_cdf.y,'Color','k','LineWidth',1);
hold on;
plot(event_cdf.x,event_cdf.y,'Color',[0.6 1 0.6],'LineWidth',1); 
% opto
test.nInSequence=[2]; 
trial1='(trialTypes.after_cue_success==1) & trialTypes.ledIsOn==1 & trialTypes.ITItoUse==1';
test.trial1=trial1;
test.trial1_name='ate after cue w opto';
test.templateSequence2_cond=eval(trial1);
trial2='(trialTypes.after_cue_success==1) & trialTypes.ITItoUse==1 & trialTypes.allTrialsForSubset==1';
test.trial2=trial2;
test.trial2_name='ate after cue';
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=false; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name='temp';
% build paired RT data set
fakeCueInd=50; % in indices
if useCorrected==false
    dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,true);
    plotReactionTimeDotScatter(alltbt,'cueZone_onVoff','all_reachBatch','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts',trialTypes,metadata,true,dataset.realDistributions.event_isSeq{1},true,true);
    title(['First trial: ' test.trial1_name]);
    ylabel(['Second trial: ' test.trial2_name]);
    [event_cdf,~]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rt_cdf',true);
    [deltaopto_event_cdf,~]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_delta_rt_cdf',true);
else
    [~,dataset]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,false);
    [event_cdf,~]=plotBehaviorEventFx(dataset,alltbt,[],'plot_rt_cdf',true);
    [deltaopto_event_cdf,~]=plotBehaviorEventFx(dataset,alltbt,[],'plot_delta_rt_cdf',true);
end
axes(ha(8));
plot(event_cdf.x,event_cdf.y,'Color','r','LineWidth',1); 
h=legend({'All trials','Ate after cue after ate after cue','Ate after cue after ate after cue opto'},'Location','best');
legend('boxoff');
set(h,'FontSize',6); 
xlabel('RT (sec)');
ylabel('CDF');
% change in rt
axes(ha(9));
plot(delta_alltrials_cdf.x,delta_alltrials_cdf.y,'Color','k','LineWidth',1);
hold on;
plot(delta_event_cdf.x,delta_event_cdf.y,'Color',[0.6 1 0.6],'LineWidth',1); 
plot(deltaopto_event_cdf.x,deltaopto_event_cdf.y,'Color','r','LineWidth',1); 
h=legend({'All trials','Ate after cue after ate after cue','Ate after cue aft ate after cue opto'},'Location','best');
legend('boxoff');
set(h,'FontSize',6); 
xlabel('Change in RT (sec)');
ylabel('CDF');



% 4. Hypothesis: cued reach doesn't touch pellet on trial N makes reaction time
% (RT) on trial N+1 slower
% -------------- shift in all trials RT distribution
cueInd=find(nanmean(alltbt.(cueName),1)>0.05,1,'first');
cueInd=num2str(cueInd);
test.nInSequence=[2]; 
trial1=['(trialTypes.after_cue_miss==1 | trialTypes.after_cue_no_pellet==1) & trialTypes.touched_pellet==0 & trialTypes.led==0 & trialTypes.ITItoUse==1'];
% trial1=['(trialTypes.after_cue_miss==1 | trialTypes.after_cue_no_pellet==1) & trialTypes.led==0 & trialTypes.ITItoUse==1'];
% trial1=['~(trialTypes.after_cue_drop==1 | trialTypes.after_cue_success==1) & any(alltbt.all_reachBatch(:,' cueInd ':end)>0.5,2) & trialTypes.led==0 & trialTypes.ITItoUse==1'];
test.trial1=trial1;
test.trial1_name='reach after cue no touch no opto';
test.templateSequence2_cond=eval(trial1);
trial2='trialTypes.allTrialsToUse==1';
test.trial2=trial2;
test.trial2_name='any';
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=false; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name='temp';
% build paired RT data set
fakeCueInd=50; % in indices
if useCorrected==false
    dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,true);
    plotReactionTimeDotScatter(alltbt,'cueZone_onVoff','all_reachBatch','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts',trialTypes,metadata,true,dataset.realDistributions.event_isSeq{1},true,true);
    title(['First trial: ' test.trial1_name]);
    ylabel(['Second trial: ' test.trial2_name]);
    [event_cdf,~]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rt_cdf',true);
else
    [~,dataset]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,false);
    [event_cdf,~]=plotBehaviorEventFx(dataset,alltbt,[],'plot_rt_cdf',true);
end
axes(ha(10));
plot(allCuedTouch.x,allCuedTouch.y,'Color',[1 0.6 0.6],'LineWidth',1);
hold on;
plot(event_cdf.x,event_cdf.y,'Color',[0.6 0.6 1],'LineWidth',1); 
% opto
test.nInSequence=[2]; 
trial1=['(trialTypes.after_cue_miss==1 | trialTypes.after_cue_no_pellet==1) & trialTypes.touched_pellet==0 & trialTypes.ledIsOn==1 & trialTypes.ITItoUse==1'];
% trial1=['(trialTypes.after_cue_miss==1 | trialTypes.after_cue_no_pellet==1) & trialTypes.ledIsOn==1 & trialTypes.ITItoUse==1'];
% trial1=['~(trialTypes.after_cue_drop==1 | trialTypes.after_cue_success==1) & any(alltbt.all_reachBatch(:,' cueInd ':end)>0.5,2) & trialTypes.ledIsOn==1 & trialTypes.ITItoUse==1'];
test.trial1=trial1;
test.trial1_name='reach after cue no touch w opto';
test.templateSequence2_cond=eval(trial1);
trial2='trialTypes.allTrialsToUse==1';
test.trial2=trial2;
test.trial2_name='any';
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=false; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name='temp';
% build paired RT data set
fakeCueInd=50; % in indices
if useCorrected==false
    dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,true);
    plotReactionTimeDotScatter(alltbt,'cueZone_onVoff','all_reachBatch','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts',trialTypes,metadata,true,dataset.realDistributions.event_isSeq{1},true,true);
    title(['First trial: ' test.trial1_name]);
    ylabel(['Second trial: ' test.trial2_name]);
    [event_cdf,~]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rt_cdf',true);
else
    [~,dataset]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,false);
    [event_cdf,~]=plotBehaviorEventFx(dataset,alltbt,[],'plot_rt_cdf',true);
end
axes(ha(10));
plot(event_cdf.x,event_cdf.y,'Color','r','LineWidth',1); 
h=legend({'All after reach after cue touch','All after reach after cue no touch','All after reach after cue no touch opto'},'Location','best');
legend('boxoff');
set(h,'FontSize',6); 
xlabel('RT (sec)');
ylabel('CDF');
% -------------- shift in matched trial type distribution
test.nInSequence=[2]; 
trial1=['(trialTypes.after_cue_miss==1 | trialTypes.after_cue_no_pellet==1) & trialTypes.touched_pellet==0 & trialTypes.led==0 & trialTypes.ITItoUse==1'];
test.trial1=trial1;
test.trial1_name='reach after cue no touch no opto';
test.templateSequence2_cond=eval(trial1);
trial2=['(trialTypes.after_cue_miss==1 | trialTypes.after_cue_no_pellet==1) & trialTypes.touched_pellet==0 & trialTypes.ITItoUse==1 & trialTypes.allTrialsForSubset==1'];
test.trial2=trial2;
test.trial2_name='reach after cue no touch';
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=false; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name='temp';
% build paired RT data set
fakeCueInd=50; % in indices
if useCorrected==false
    dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,true);
    plotReactionTimeDotScatter(alltbt,'cueZone_onVoff','all_reachBatch','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts',trialTypes,metadata,true,dataset.realDistributions.event_isSeq{1},true,true);
    title(['First trial: ' test.trial1_name]);
    ylabel(['Second trial: ' test.trial2_name]);
    [event_cdf,alltrials_cdf]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rt_cdf',true);
    [delta_event_cdf,delta_alltrials_cdf]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_delta_rt_cdf',true);
else
    [~,dataset]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,false);
    [event_cdf,alltrials_cdf]=plotBehaviorEventFx(dataset,alltbt,[],'plot_rt_cdf',true);
    [delta_event_cdf,delta_alltrials_cdf]=plotBehaviorEventFx(dataset,alltbt,[],'plot_delta_rt_cdf',true);
end
axes(ha(11));
plot(matchedCuedTouch.x,matchedCuedTouch.y,'Color',[1 0.6 0.6],'LineWidth',1);
hold on;
plot(event_cdf.x,event_cdf.y,'Color',[0.6 0.6 1],'LineWidth',1); 
% opto
test.nInSequence=[2]; 
trial1=['(trialTypes.after_cue_miss==1 | trialTypes.after_cue_no_pellet==1) & trialTypes.touched_pellet==0 & trialTypes.ledIsOn==1 & trialTypes.ITItoUse==1'];
test.trial1=trial1;
test.trial1_name='reach after cue no touch w opto';
test.templateSequence2_cond=eval(trial1);
trial2=['(trialTypes.after_cue_miss==1 | trialTypes.after_cue_no_pellet==1) & trialTypes.touched_pellet==0 & trialTypes.ITItoUse==1 & trialTypes.allTrialsForSubset==1'];
test.trial2=trial2;
test.trial2_name='reach after cue no touch';
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=false; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name='temp';
% build paired RT data set
fakeCueInd=50; % in indices
if useCorrected==false
    dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,true);
    plotReactionTimeDotScatter(alltbt,'cueZone_onVoff','all_reachBatch','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts',trialTypes,metadata,true,dataset.realDistributions.event_isSeq{1},true,true);
    title(['First trial: ' test.trial1_name]);
    ylabel(['Second trial: ' test.trial2_name]);
    [event_cdf,~]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rt_cdf',true);
    [deltaopto_event_cdf,~]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_delta_rt_cdf',true);
else
    [~,dataset]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,false);
    [event_cdf,~]=plotBehaviorEventFx(dataset,alltbt,[],'plot_rt_cdf',true);
    [deltaopto_event_cdf,~]=plotBehaviorEventFx(dataset,alltbt,[],'plot_delta_rt_cdf',true);
end
axes(ha(11));
plot(event_cdf.x,event_cdf.y,'Color','r','LineWidth',1); 
h=legend({'Reach after cue touch after reach after cue touch','Reach after cue no touch after reach after cue no touch','Reach after cue no touch after reach after cue no touch opto'},'Location','best');
legend('boxoff');
set(h,'FontSize',6); 
xlabel('RT (sec)');
ylabel('CDF');
% change in rt
axes(ha(12));
plot(deltaCuedTouch.x,deltaCuedTouch.y,'Color',[1 0.6 0.6],'LineWidth',1);
hold on;
plot(delta_event_cdf.x,delta_event_cdf.y,'Color',[0.6 0.6 1],'LineWidth',1); 
plot(deltaopto_event_cdf.x,deltaopto_event_cdf.y,'Color','r','LineWidth',1); 
h=legend({'Reach after cue touch after reach after cue touch','Reach after cue no touch after reach after cue no touch','Reach after cue no touch after reach after cue no touch opto'},'Location','best');
legend('boxoff');
set(h,'FontSize',6); 
xlabel('Change in RT (sec)');
ylabel('CDF');




% 5. Hypothesis: failing to reach after cue on trial N does not affect reaction time
% (RT) on trial N+1 
% -------------- shift in all trials RT distribution
cueInd=find(nanmean(alltbt.(cueName),1)>0.05,1,'first');
cueInd=num2str(cueInd);
test.nInSequence=[2]; 
trial1=['~any(alltbt.all_reachBatch(:,' cueInd ':end)>0.5,2) & trialTypes.led==0 & trialTypes.ITItoUse==1'];
test.trial1=trial1;
test.trial1_name='no reach after cue no opto';
test.templateSequence2_cond=eval(trial1);
trial2='trialTypes.allTrialsToUse==1';
test.trial2=trial2;
test.trial2_name='any';
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=false; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name='temp';
% build paired RT data set
fakeCueInd=50; % in indices
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,true);
plotReactionTimeDotScatter(alltbt,'cueZone_onVoff','all_reachBatch','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts',trialTypes,metadata,true,dataset.realDistributions.event_isSeq{1},true,true);
title(['First trial: ' test.trial1_name]);
ylabel(['Second trial: ' test.trial2_name]);
[event_cdf,~]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rt_cdf',true);
axes(ha(13));
plot(allCuedTouch.x,allCuedTouch.y,'Color',[1 0.6 0.6],'LineWidth',1);
hold on;
plot(event_cdf.x,event_cdf.y,'Color',[0.6 0.6 1],'LineWidth',1); 
% opto
test.nInSequence=[2]; 
trial1=['~any(alltbt.all_reachBatch(:,' cueInd ':end)>0.5,2) & trialTypes.ledIsOn==1 & trialTypes.ITItoUse==1'];
test.trial1=trial1;
test.trial1_name='no reach after cue w opto';
test.templateSequence2_cond=eval(trial1);
trial2='trialTypes.allTrialsToUse==1';
test.trial2=trial2;
test.trial2_name='any';
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=false; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name='temp';
% build paired RT data set
fakeCueInd=50; % in indices
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,true);
plotReactionTimeDotScatter(alltbt,'cueZone_onVoff','all_reachBatch','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts',trialTypes,metadata,true,dataset.realDistributions.event_isSeq{1},true,true);
title(['First trial: ' test.trial1_name]);
ylabel(['Second trial: ' test.trial2_name]);
[event_cdf,~]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rt_cdf',true);
axes(ha(13));
plot(event_cdf.x,event_cdf.y,'Color','r','LineWidth',1); 
h=legend({'All after reach after cue touch','All after no reach after cue','All aft no reach after cue opto'},'Location','best');
legend('boxoff');
set(h,'FontSize',6); 
xlabel('RT (sec)');
ylabel('CDF');




% ---------- effects of opto on current trial
test.nInSequence=[2]; 
trial1='trialTypes.allTrialsToUse==1';
test.trial1=trial1;
test.trial1_name='any';
test.templateSequence2_cond=eval(trial1);
trial2='trialTypes.led==0';
test.trial2=trial2;
test.trial2_name='current trial opto off';
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=false; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name='temp';
% build paired RT data set
fakeCueInd=50; % in indices
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,true); 
plotReactionTimeDotScatter(alltbt,'cueZone_onVoff','all_reachBatch','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts',trialTypes,metadata,true,dataset.realDistributions.event_isSeq{1},true,true);
title(['First trial: ' test.trial1_name]);
ylabel(['Second trial: ' test.trial2_name]);
[event_cdf,alltrials_cdf]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rt_cdf',true);
[delta_event_cdf,delta_alltrials_cdf]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_delta_rt_cdf',true);
axes(ha(14));
plot(alltrials_cdf.x,alltrials_cdf.y,'Color','k','LineWidth',1);
hold on;
plot(event_cdf.x,event_cdf.y,'Color',[0.7 0.7 0.7],'LineWidth',1); 
% opto
test.nInSequence=[2]; 
trial1='trialTypes.allTrialsToUse==1';
test.trial1=trial1;
test.trial1_name='any';
test.templateSequence2_cond=eval(trial1);
trial2='trialTypes.ledIsOn==1';
test.trial2=trial2;
test.trial2_name='current trial opto on';
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=false; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name='temp';
% build paired RT data set
fakeCueInd=50; % in indices
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,true); 
plotReactionTimeDotScatter(alltbt,'cueZone_onVoff','all_reachBatch','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts',trialTypes,metadata,true,dataset.realDistributions.event_isSeq{1},true,true);
title(['First trial: ' test.trial1_name]);
ylabel(['Second trial: ' test.trial2_name]);
[event_cdf,~]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rt_cdf',true);
[deltaopto_event_cdf,~]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_delta_rt_cdf',true);
axes(ha(14));
plot(event_cdf.x,event_cdf.y,'Color','r','LineWidth',1); 
h=legend({'All trials','All during opto off','All during opto on'},'Location','best');
legend('boxoff');
set(h,'FontSize',6); 
xlabel('RT (sec)');
ylabel('CDF');
% change in rt
axes(ha(15));
plot(delta_alltrials_cdf.x,delta_alltrials_cdf.y,'Color','k','LineWidth',1);
hold on;
plot(delta_event_cdf.x,delta_event_cdf.y,'Color',[0.7 0.7 0.7],'LineWidth',1); 
plot(deltaopto_event_cdf.x,deltaopto_event_cdf.y,'Color','r','LineWidth',1); 
h=legend({'All trials','All during opto off','All during opto on'},'Location','best');
legend('boxoff');
set(h,'FontSize',6); 
xlabel('Change in RT (sec)');
ylabel('CDF');