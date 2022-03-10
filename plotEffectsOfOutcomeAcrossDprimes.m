function plotEffectsOfOutcomeAcrossDprimes(alltbt,trialTypes,metadata)

dprimes_steps={[-100 0.3+0.01],[0.3-0.01 0.6+0.01],[0.6-0.01 0.9+0.01],[0.9-0.01 1.2+0.01],[1.2-0.01 1.6+0.01],[1.6-0.01 100+0.01]};
trialTypes.isLongITI_1back=[0; trialTypes.isLongITI(1:end-1)];
% dprimes_steps={[-100 0.25],[0.25 0.5],[0.5 0.75],[0.75 1],[1 1.25],[1.25 1.5],[1.5 1.75],[1.75 2],[2 2.25],[2.25 2.5],[2.5 100]};
%trial1backup='trialTypes.touch_in_cued_window_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.optoGroup~=1 & trialTypes.isLongITI_1forward==1';
%trial1backup='trialTypes.cued_reach_1forward==1 & trialTypes.touched_pellet_1forward==0 & trialTypes.optoGroup~=1 & trialTypes.optoGroup_1forward~=1';
%trial1backup='trialTypes.cued_reach_1forward==0  & trialTypes.touched_pellet_1forward==1 & trialTypes.optoGroup~=1 & trialTypes.isLongITI_1forward==1';
trial1backup='trialTypes.touch_in_cued_window_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.optoGroup~=1 & trialTypes.isLongITI_1forward==1';

backup.alltbt=alltbt; backup.trialTypes=trialTypes; backup.metadata=metadata;

boots_trial1_realmean_x=cell(1,length(dprimes_steps));
boots_trial1_realmean_y=cell(1,length(dprimes_steps));
boots_trialn_realmean_x=cell(1,length(dprimes_steps));
boots_trialn_realmean_y=cell(1,length(dprimes_steps));
boots_trial1=cell(1,length(dprimes_steps));
boots_trialn=cell(1,length(dprimes_steps));
boots_trial1_nooutliers=cell(1,length(dprimes_steps));
boots_trial1_realmean_x_nooutliers=cell(1,length(dprimes_steps));
boots_trial1_realmean_y_nooutliers=cell(1,length(dprimes_steps));

LEDboots_trial1_realmean_x=cell(1,length(dprimes_steps));
LEDboots_trial1_realmean_y=cell(1,length(dprimes_steps));
LEDboots_trialn_realmean_x=cell(1,length(dprimes_steps));
LEDboots_trialn_realmean_y=cell(1,length(dprimes_steps));
LEDboots_trial1=cell(1,length(dprimes_steps));
LEDboots_trialn=cell(1,length(dprimes_steps));

for i=1:length(dprimes_steps)
    curr_dprimes_steps=dprimes_steps{i};
    
    alltbt=backup.alltbt; trialTypes=backup.trialTypes; metadata=backup.metadata;
    
    %%%%%%%%%%%%%%%%%%%%% perform any filtering on alltbt
    % for example, filter by d-prime
    
    temp=datestr(datetime('now'));
    temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))='';
    temp=temp(~isspace(temp));
    saveDir=['/Volumes/Neurobio/MICROSCOPE/Kim/RT pairs data sets/' temp]; % where to save details of alltbt filtering and RT pairs data set
    
    % filter settings
    tbt_filter.sortField='dprimes';
    tbt_filter.range_values=curr_dprimes_steps;
    tbt_filter.name=[tbt_filter.sortField num2str(tbt_filter.range_values(1)) 'to' num2str(tbt_filter.range_values(2))];
    temp=tbt_filter.name;
    temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))='';
    temp=temp(~isspace(temp));
    tbt_filter.name=temp;
    tbt_filter.clock_progress=true;
    
    % filter alltbt
    [alltbt,trialTypes,metadata]=filtTbt(alltbt,trialTypes,tbt_filter.sortField,tbt_filter.range_values,metadata,tbt_filter.clock_progress);
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%% build relevant data sets
    
    % settings for paired RT data set
    test.nInSequence=[3]; % defines trial pairs, e.g., 2 means will compare each trial with its subsequent trial, 3 means will compare each trial with the trial after next, etc.
    % requirement for first trial in pair
    % trial1='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1';
    trialTypes.isLongITI_1forward=[trialTypes.isLongITI(2:end); 0];
    trialTypes.optoGroup_1forward=[trialTypes.optoGroup(2:end); 0];
    trial1=[trial1backup ' & (trialTypes.led_1forward==0)'];
    test.trial1=trial1;
    test.templateSequence2_cond=eval(trial1);
    % trial2='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1';
    trial2='trialTypes.optoGroup~=1';
    test.trial2=trial2;
    test.templateSequence2_end=eval(trial2);
    test.fillInBetweenWithAnything=true; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
    test.event_name=['alltrials' tbt_filter.name 'inBetweenAnything' num2str(test.fillInBetweenWithAnything)];
    saveDir2=[];
    % build paired RT data set
    fakeCueInd=50; % in indices, this is not relevant if not using PCA-based RT model
    skipCorrected=true;
    % this function builds the dataset using the trial type sequences specified above
    [dataset,correctedDistributions]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir2,test,skipCorrected);
    
    %%%%%%%%%%%%%%%%%%%%% measure how reach rate changes over course of session
    
    shuffleTrialOrder=false; % if want to randomly permute trial order to test for ordering effects
    reachratesettings.epsilon_cue=0; % in seconds
    reachratesettings.epsilon_uncue=2; % in seconds
    reachratesettings.epsilon_beforecue=1; % in seconds
    reachratesettings.percentOfReachesFromSess_forInitCond=20; % use this fraction of trials from beginning of session to get initial windows
    reachratesettings.percentOfReachesFromSess_forInitRate=20; % use this fraction of trials from beginning of session to get initial reach rates
    reachratesettings.maxTrialLength=9.5; % in sec, wrt cue
    reachratesettings.minTrialLength=-2; % wrt cue, in sec
    reachratesettings.suppressPlots=false;
    % sec wrt cue onset
    reachratesettings.acrossSess_window1=[0.05 1]; % cued window [0.05 1]
    % reachratesettings.acrossSess_window1=[4 7];
    % note that after mouse gets a pellet, reaching is suppressed
    reachratesettings.acrossSess_window2=[7 reachratesettings.maxTrialLength]; % beware reach suppression after a success
    reachratesettings.acrossSess_window3=[reachratesettings.minTrialLength -1];
    reachratesettings.scatterPointSize=50; % size for points in scatter plot
    reachratesettings.addSatietyLines=false; % whether to add proportionality lines to figure
    reachratesettings.useWindowsForUncued=[3]; % to use window2 or window3 or both for the uncued reach rate
    reachratesettings.initWindows=[]; % empty if want to calculate from dataset
    reachratesettings.useRateMethod=1; % 1, 2 or 3 (see explanation below)
    
    reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings);
    close all;
    boots_trial1_realmean_x{i}=reachrates.boot_trial1_realmean_x;
    boots_trial1_realmean_y{i}=reachrates.boot_trial1_realmean_y;
    boots_trialn_realmean_x{i}=reachrates.boot_trialn_realmean_x;
    boots_trialn_realmean_y{i}=reachrates.boot_trialn_realmean_y;
    boots_trial1{i}=reachrates.boot_trial1;
    boots_trialn{i}=reachrates.boot_trialn;
    
    boots_trial1_nooutliers{i}=rmoutliers(reachrates.boot_trial1);
    temp=boots_trial1_nooutliers{i};
    boots_trial1_realmean_x_nooutliers{i}=nanmean(temp(1,:));
    boots_trial1_realmean_y_nooutliers{i}=nanmean(temp(2,:));
    
    
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%% build relevant data sets
    trial1=[trial1backup ' & (trialTypes.led_1forward==1)'];
    test.trial1=trial1;
    test.templateSequence2_cond=eval(trial1);
    [dataset,correctedDistributions]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir2,test,skipCorrected);
    
    %%%%%%%%%%%%%%%%%%%%% measure how reach rate changes over course of session
    reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings);
    close all;
    LEDboots_trial1_realmean_x{i}=reachrates.boot_trial1_realmean_x;
    LEDboots_trial1_realmean_y{i}=reachrates.boot_trial1_realmean_y;
    LEDboots_trialn_realmean_x{i}=reachrates.boot_trialn_realmean_x;
    LEDboots_trialn_realmean_y{i}=reachrates.boot_trialn_realmean_y;
    LEDboots_trial1{i}=reachrates.boot_trial1;
    LEDboots_trialn{i}=reachrates.boot_trialn;
    
end

% Plot effects of trial outcome
figure();
cmap=colormap('cool');
k=1;
kstep=ceil(size(cmap,1)/length(dprimes_steps));
for i=1:length(dprimes_steps)
    bootMeans_trial1=boots_trial1_nooutliers{i};
    bootMeans_trialn=boots_trialn{i};
    %s=scatter(bootMeans_trial1(1,:),bootMeans_trial1(2,:),20,'MarkerEdgeColor','k','MarkerFaceColor',cmap(k,:),'LineWidth',0.1);
    hold on;
    %s.AlphaData = 0.5*ones(1,size(bootMeans_trial1,2));
    %s.MarkerFaceAlpha = 'flat';
    s=scatter(boots_trial1_realmean_x_nooutliers{i},boots_trial1_realmean_y_nooutliers{i},50,'k','filled');
    
    s=scatter(bootMeans_trialn(1,:),bootMeans_trialn(2,:),20,cmap(k,:),'filled');
    hold on;
    s.AlphaData = 0.5*ones(1,size(bootMeans_trialn,2));
    s.MarkerFaceAlpha = 'flat';
    s=scatter(boots_trialn_realmean_x{i},boots_trialn_realmean_y{i},50,cmap(k,:),'filled');
    line([boots_trial1_realmean_x_nooutliers{i} boots_trialn_realmean_x{i}],[boots_trial1_realmean_y_nooutliers{i} boots_trialn_realmean_y{i}],'Color','k');
    
    k=k+kstep;
    if k>size(cmap,1)
        k=size(cmap,1);
    end
end
xlabel('Uncued reach rate (1/sec)');
ylabel('Cued reach rate (1/sec)');

% Plot effects of LED
figure();
cmap=colormap('cool');
k=1;
kstep=ceil(size(cmap,1)/length(dprimes_steps));
for i=1:length(dprimes_steps)
    LEDbootMeans_trialn=LEDboots_trialn{i};
    bootMeans_trialn=boots_trialn{i};
    %s=scatter(bootMeans_trialn(1,:),bootMeans_trialn(2,:),20,cmap(k,:),'filled');
    hold on;
    %s.AlphaData = 0.5*ones(1,size(bootMeans_trialn,2));
    %s.MarkerFaceAlpha = 'flat';
    %s=scatter(boots_trialn_realmean_x{i},boots_trialn_realmean_y{i},50,cmap(k,:),'filled');
  
    s=scatter(LEDbootMeans_trialn(1,:),LEDbootMeans_trialn(2,:),20,'MarkerEdgeColor','r','MarkerFaceColor',cmap(k,:),'LineWidth',0.1);
    hold on;
    s.AlphaData = 0.5*ones(1,size(LEDbootMeans_trialn,2));
    s.MarkerFaceAlpha = 'flat';
    s=scatter(LEDboots_trialn_realmean_x{i},LEDboots_trialn_realmean_y{i},50,'r','filled');
    line([boots_trialn_realmean_x{i} LEDboots_trialn_realmean_x{i}],[boots_trialn_realmean_y{i} LEDboots_trialn_realmean_y{i}],'Color','r');
    
    k=k+kstep;
    if k>size(cmap,1)
        k=size(cmap,1);
    end
end
xlabel('Uncued reach rate (1/sec)');
ylabel('Cued reach rate (1/sec)');

end




