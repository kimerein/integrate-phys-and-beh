function plotMatchedTrialTypeShiftOnly(alltbt,trialTypes,metadata,ledIsOn,useLongITIonly,excludeFromAllTrials,alsoExcludeFromSubset,useCorrected,suppressScatters)

cueName='cueZone_onVoff';

% useCorrected==true if want to use corrected input distributions

% Set up figure layout
figure();
Nh=3; % number of rows
Nw=4; % number of columns
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

% 2. Hypothesis: touched pellet after cue on trial N makes reaction time (RT) on trial N+1 faster
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
    if suppressScatters==false
        plotReactionTimeDotScatter(alltbt,'cueZone_onVoff','all_reachBatch','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts',trialTypes,metadata,true,dataset.realDistributions.event_isSeq{1},true,true);
        title(['First trial: ' test.trial1_name]);
        ylabel(['Second trial: ' test.trial2_name]);
    end
    [event_cdf,alltrials_cdf,eventTrial1_cdf]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rt_cdf',true);
    [delta_event_cdf,delta_alltrials_cdf]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_delta_rt_cdf',true);
else
    [~,dataset]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,false);
    if suppressScatters==false
        [event_cdf,alltrials_cdf]=plotBehaviorEventFx(dataset,alltbt,[],'plot_rt_cdf',true);
        % [delta_event_cdf,delta_alltrials_cdf]=plotBehaviorEventFx(dataset,alltbt,[],'plot_delta_rt_cdf',true);
    end
end
axes(ha(1));
plot(alltrials_cdf.x,alltrials_cdf.y,'Color','k','LineWidth',1);
hold on;
plot(event_cdf.x,event_cdf.y,'Color',[1 0.6 0.6],'LineWidth',1); 
plot(eventTrial1_cdf.x,eventTrial1_cdf.y,'Color',[1 0.6 0.2],'LineWidth',1); 
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
    if suppressScatters==false
        plotReactionTimeDotScatter(alltbt,'cueZone_onVoff','all_reachBatch','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts',trialTypes,metadata,true,dataset.realDistributions.event_isSeq{1},true,true);
        title(['First trial: ' test.trial1_name]);
        ylabel(['Second trial: ' test.trial2_name]);
    end
    [event_cdf,~]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rt_cdf',true);
    [deltaopto_event_cdf,~]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_delta_rt_cdf',true);
else
    [~,dataset]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,false);
    [event_cdf,~]=plotBehaviorEventFx(dataset,alltbt,[],'plot_rt_cdf',true);
    % [deltaopto_event_cdf,~]=plotBehaviorEventFx(dataset,alltbt,[],'plot_delta_rt_cdf',true);
end
axes(ha(1));
plot(event_cdf.x,event_cdf.y,'Color','r','LineWidth',1); 
h=legend({'All trials','Touch after cue after touch after cue','Touch after cue trial 1','Touch after cue after touch after cue opto'},'Location','best');
legend('boxoff');
set(h,'FontSize',6); 
xlabel('RT (sec)');
ylabel('CDF');
% change in rt
axes(ha(2));
plot(delta_alltrials_cdf.x,delta_alltrials_cdf.y,'Color','k','LineWidth',1);
hold on; 
plot(delta_event_cdf.x,delta_event_cdf.y,'Color',[1 0.6 0.6],'LineWidth',1); 
deltaCuedTouch=delta_event_cdf;
plot(deltaopto_event_cdf.x,deltaopto_event_cdf.y,'Color','r','LineWidth',1);
if useCorrected==false
    h=legend({'All trials','Touch after cue after touch after cue','Touch after cue aft touch after cue opto'},'Location','best');
else
    h=legend({'All trials','All after touch after cue','All after touch after cue opto'},'Location','best');
end
legend('boxoff');
set(h,'FontSize',6); 
xlabel('Change in RT (sec)');
ylabel('CDF');
% difference between CDFs of RT changes
axes(ha(3));
plot(delta_alltrials_cdf.x,deltaopto_event_cdf.y-delta_event_cdf.y,'Color',[1 0.2 0.2],'LineWidth',1);
hold on;
line([0 0],[nanmin(deltaopto_event_cdf.y-delta_event_cdf.y) nanmax(deltaopto_event_cdf.y-delta_event_cdf.y)],'Color',[0.8 0.8 0.8]);
line([5 5],[nanmin(deltaopto_event_cdf.y-delta_event_cdf.y) nanmax(deltaopto_event_cdf.y-delta_event_cdf.y)],'Color',[0.8 0.8 0.8]);
touchchange_y=deltaopto_event_cdf.y-delta_event_cdf.y;
touchchange_x=delta_alltrials_cdf.x;





% 4. Hypothesis: cued reach doesn't touch pellet on trial N makes reaction time
% (RT) on trial N+1 slower
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
    if suppressScatters==false
        plotReactionTimeDotScatter(alltbt,'cueZone_onVoff','all_reachBatch','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts',trialTypes,metadata,true,dataset.realDistributions.event_isSeq{1},true,true);
        title(['First trial: ' test.trial1_name]);
        ylabel(['Second trial: ' test.trial2_name]);
    end
    [event_cdf,alltrials_cdf,eventTrial1_cdf]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rt_cdf',true);
    [delta_event_cdf,delta_alltrials_cdf]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_delta_rt_cdf',true);
else
    [~,dataset]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,false);
    [event_cdf,alltrials_cdf]=plotBehaviorEventFx(dataset,alltbt,[],'plot_rt_cdf',true);
    % [delta_event_cdf,delta_alltrials_cdf]=plotBehaviorEventFx(dataset,alltbt,[],'plot_delta_rt_cdf',true);
end
axes(ha(5));
plot(matchedCuedTouch.x,matchedCuedTouch.y,'Color',[1 0.6 0.6],'LineWidth',1);
hold on;
plot(event_cdf.x,event_cdf.y,'Color',[0.6 0.6 1],'LineWidth',1); 
plot(eventTrial1_cdf.x,eventTrial1_cdf.y,'Color','c','LineWidth',1); 
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
    if suppressScatters==false
        plotReactionTimeDotScatter(alltbt,'cueZone_onVoff','all_reachBatch','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts',trialTypes,metadata,true,dataset.realDistributions.event_isSeq{1},true,true);
        title(['First trial: ' test.trial1_name]);
        ylabel(['Second trial: ' test.trial2_name]);
    end
    [event_cdf,~]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rt_cdf',true);
    [deltaopto_event_cdf,~]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_delta_rt_cdf',true);
else
    [~,dataset]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,false);
    [event_cdf,~]=plotBehaviorEventFx(dataset,alltbt,[],'plot_rt_cdf',true);
    % [deltaopto_event_cdf,~]=plotBehaviorEventFx(dataset,alltbt,[],'plot_delta_rt_cdf',true);
end
axes(ha(5));
plot(event_cdf.x,event_cdf.y,'Color','r','LineWidth',1); 
h=legend({'Reach after cue touch after reach after cue touch','Reach after cue no touch after reach after cue no touch','Reach after cue no touch first trial','Reach after cue no touch after reach after cue no touch opto'},'Location','best');
legend('boxoff');
set(h,'FontSize',6); 
xlabel('RT (sec)');
ylabel('CDF');
% change in rt
axes(ha(6));
plot(deltaCuedTouch.x,deltaCuedTouch.y,'Color',[1 0.6 0.6],'LineWidth',1);
hold on;
plot(delta_event_cdf.x,delta_event_cdf.y,'Color',[0.6 0.6 1],'LineWidth',1); 
plot(deltaopto_event_cdf.x,deltaopto_event_cdf.y,'Color','r','LineWidth',1); 
if useCorrected==false
    h=legend({'Reach after cue touch after reach after cue touch','Reach after cue no touch after reach after cue no touch','Reach after cue no touch after reach after cue no touch opto'},'Location','best');
else
    h=legend({'All after reach after cue touch','All after reach after cue no touch','All after reach after cue no touch opto'},'Location','best');
end
legend('boxoff');
set(h,'FontSize',6); 
xlabel('Change in RT (sec)');
ylabel('CDF');
axes(ha(7));
plot(delta_alltrials_cdf.x,deltaopto_event_cdf.y-delta_event_cdf.y,'Color',[0.2 0.2 1],'LineWidth',1);
misschange_y=deltaopto_event_cdf.y-delta_event_cdf.y;
misschange_x=delta_alltrials_cdf.x;



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
if suppressScatters==false
    plotReactionTimeDotScatter(alltbt,'cueZone_onVoff','all_reachBatch','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts',trialTypes,metadata,true,dataset.realDistributions.event_isSeq{1},true,true);
    title(['First trial: ' test.trial1_name]);
    ylabel(['Second trial: ' test.trial2_name]);
end
[event_cdf,alltrials_cdf]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rt_cdf',true);
[delta_event_cdf,delta_alltrials_cdf]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_delta_rt_cdf',true);
axes(ha(9));
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
if suppressScatters==false
    plotReactionTimeDotScatter(alltbt,'cueZone_onVoff','all_reachBatch','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts',trialTypes,metadata,true,dataset.realDistributions.event_isSeq{1},true,true);
    title(['First trial: ' test.trial1_name]);
    ylabel(['Second trial: ' test.trial2_name]);
end
[event_cdf,~]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rt_cdf',true);
[deltaopto_event_cdf,~]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_delta_rt_cdf',true);
axes(ha(9));
plot(event_cdf.x,event_cdf.y,'Color','r','LineWidth',1); 
h=legend({'All trials','All during opto off','All during opto on'},'Location','best');
legend('boxoff');
set(h,'FontSize',6); 
xlabel('RT (sec)');
ylabel('CDF');
% change in rt
axes(ha(10));
plot(delta_alltrials_cdf.x,delta_alltrials_cdf.y,'Color','k','LineWidth',1);
hold on;
plot(delta_event_cdf.x,delta_event_cdf.y,'Color',[0.7 0.7 0.7],'LineWidth',1); 
plot(deltaopto_event_cdf.x,deltaopto_event_cdf.y,'Color','r','LineWidth',1); 
h=legend({'All trials','All during opto off','All during opto on'},'Location','best');
legend('boxoff');
set(h,'FontSize',6); 
xlabel('Change in RT (sec)');
ylabel('CDF');







% ---------- effects of opto on next trial
test.nInSequence=[2]; 
trial1='trialTypes.led==0';
test.trial1=trial1;
test.trial1_name='any';
test.templateSequence2_cond=eval(trial1);
trial2='trialTypes.allTrialsToUse==1';
test.trial2=trial2;
test.trial2_name='current trial opto off';
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=false; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name='temp';
% build paired RT data set
fakeCueInd=50; % in indices
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,true);
if suppressScatters==false
    plotReactionTimeDotScatter(alltbt,'cueZone_onVoff','all_reachBatch','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts',trialTypes,metadata,true,dataset.realDistributions.event_isSeq{1},true,true);
    title(['First trial: ' test.trial1_name]);
    ylabel(['Second trial: ' test.trial2_name]);
end
% [event_cdf,alltrials_cdf]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rt_cdf',true);
[delta_event_cdf,delta_alltrials_cdf]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_delta_rt_cdf',true);
% axes(ha(9));
% plot(alltrials_cdf.x,alltrials_cdf.y,'Color','k','LineWidth',1);
% hold on;
% plot(event_cdf.x,event_cdf.y,'Color',[0.7 0.7 0.7],'LineWidth',1); 
% opto
test.nInSequence=[2]; 
trial1='trialTypes.ledIsOn==1';
test.trial1=trial1;
test.trial1_name='any';
test.templateSequence2_cond=eval(trial1);
trial2='trialTypes.allTrialsToUse==1';
test.trial2=trial2;
test.trial2_name='current trial opto on';
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=false; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name='temp';
% build paired RT data set
fakeCueInd=50; % in indices
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,[],test,true); 
if suppressScatters==false
    plotReactionTimeDotScatter(alltbt,'cueZone_onVoff','all_reachBatch','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts',trialTypes,metadata,true,dataset.realDistributions.event_isSeq{1},true,true);
    title(['First trial: ' test.trial1_name]);
    ylabel(['Second trial: ' test.trial2_name]);
end
% [event_cdf,~]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rt_cdf',true);
[deltaopto_event_cdf,~]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_delta_rt_cdf',true);
% axes(ha(9));
% plot(event_cdf.x,event_cdf.y,'Color','r','LineWidth',1); 
% h=legend({'All trials','All during opto off','All during opto on'},'Location','best');
% legend('boxoff');
% set(h,'FontSize',6); 
% xlabel('RT (sec)');
% ylabel('CDF');
% change in rt
axes(ha(11));
plot(delta_alltrials_cdf.x,delta_alltrials_cdf.y,'Color','k','LineWidth',1);
hold on;
plot(delta_event_cdf.x,delta_event_cdf.y,'Color',[0.7 0.7 0.7],'LineWidth',1); 
plot(deltaopto_event_cdf.x,deltaopto_event_cdf.y,'Color','r','LineWidth',1); 
h=legend({'All trials','All after opto off','All after opto on'},'Location','best');
legend('boxoff');
set(h,'FontSize',6); 
xlabel('Change in RT (sec)');
ylabel('CDF');
allchange_y=deltaopto_event_cdf.y-delta_event_cdf.y;
allchange_x=delta_alltrials_cdf.x;




% touch change minus all change
axes(ha(4));
[new_n,new_x]=cityscape_hist(touchchange_y-allchange_y,touchchange_x);
plot(new_x,new_n,'Color',[1 0.4 0.4]);

% miss change minus all change
axes(ha(8));
[new_n,new_x]=cityscape_hist(misschange_y-allchange_y,touchchange_x);
plot(new_x,new_n,'Color',[0.4 0.4 1]);
hold on;
line([0 0],[nanmin(misschange_y-allchange_y) nanmax(misschange_y-allchange_y)],'Color',[0.8 0.8 0.8]);
line([5 5],[nanmin(misschange_y-allchange_y) nanmax(misschange_y-allchange_y)],'Color',[0.8 0.8 0.8]);
hold on;
line([-15 -15],[nanmin(misschange_y-allchange_y) nanmax(misschange_y-allchange_y)],'Color',[0.8 0.8 0.8]);
line([0 0],[nanmin(misschange_y-allchange_y) nanmax(misschange_y-allchange_y)],'Color',[0.8 0.8 0.8]);
disp(['Value in 0 to 5 sec window for RT change CDF difference: ' num2str(nanmean(touchchange_y(touchchange_x>0 & touchchange_x<=5)))]);
disp(['Value in 0 to 5 sec window for RT change CDF difference (miss): ' num2str(nanmean(misschange_y(touchchange_x>0 & touchchange_x<=5)-allchange_y(touchchange_x>0 & touchchange_x<=5)))]);
disp(['Value in -15 to 0 sec window for RT change CDF difference (miss): ' num2str(nanmean(misschange_y(touchchange_x>=-15 & touchchange_x<0)-allchange_y(touchchange_x>=-15 & touchchange_x<0)))]);


% all change
axes(ha(12));
[new_n,new_x]=cityscape_hist(allchange_y,allchange_x);
plot(new_x,new_n,'Color','r');