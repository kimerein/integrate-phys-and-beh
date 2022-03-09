function plotWithinSession_and_dayByDay(alltbt,metadata,trialTypes,startSessId,endSessId,saveFigDir,imgType,outlierTest)

days=startSessId:endSessId;
daybyday_uncued=cell(1,length(days));
daybyday_cued=cell(1,length(days));
backup_alltbt=alltbt;
backup_metadata=metadata;
backup_trialTypes=trialTypes;
maxuncued=0;
maxcued=0;
for i=1:length(days)
    alltbt=backup_alltbt;
    metadata=backup_metadata;
    trialTypes=backup_trialTypes;
    day=days(i);
    disp(['Day ' i]);

    % perform any filtering on alltbt
    % for example, filter by d-prime
    
    temp=datestr(datetime('now'));
    temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))='';
    temp=temp(~isspace(temp));
    saveDir=['/Volumes/Neurobio/MICROSCOPE/Kim/RT pairs data sets/' temp]; % where to save details of alltbt filtering and RT pairs data set
    
    % filter settings
    % tbt_filter.sortField='mouseid';
    tbt_filter.sortField='sessid';
    % tbt_filter.range_values=[1 2 6 9 10 11 12 18];
    % tbt_filter.range_values=[1     2     3     6     7     8     9    10    11    12    14    15    17    18];
    tbt_filter.range_values=[day-0.5 day+0.5];
    % tbt_filter.range_values=[10.5 11.5];
    tbt_filter.name=[tbt_filter.sortField num2str(tbt_filter.range_values(1)) 'to' num2str(tbt_filter.range_values(2))];
    temp=tbt_filter.name;
    temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))='';
    temp=temp(~isspace(temp));
    tbt_filter.name=temp;
    tbt_filter.clock_progress=true;
    
    % filter alltbt
    [alltbt,trialTypes,metadata]=filtTbt(alltbt,trialTypes,tbt_filter.sortField,tbt_filter.range_values,metadata,tbt_filter.clock_progress);
    
    % build relevant data sets
    
    % settings for paired RT data set
    test.nInSequence=[2]; % defines trial pairs, e.g., 2 means will compare each trial with its subsequent trial, 3 means will compare each trial with the trial after next, etc.
    % requirement for first trial in pair
    % trial1='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1';
    trialTypes.isLongITI_1forward=[trialTypes.isLongITI(2:end); 0];
    trialTypes.optoGroup_1forward=[trialTypes.optoGroup(2:end); 0];
    % trial1='trialTypes.optoGroup~=1';
    trial1='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1';
    % trial1='trialTypes.after_cue_success_1forward==1 & trialTypes.consumed_pellet_1forward==1 & trialTypes.led_1forward==0 & trialTypes.optoGroup_1forward~=1'; % & trialTypes.isLongITI_1forward==1'];
    % trial1='trialTypes.touch_in_cued_window_1forward==1 & trialTypes.led_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.optoGroup~=1  & trialTypes.isLongITI_1forward==1';
    % trial1='trialTypes.cued_reach_1forward==1 & trialTypes.touched_pellet_1forward==1 & (trialTypes.led_1forward==0) & trialTypes.optoGroup~=1  & trialTypes.optoGroup_1forward~=1';
    % trial1='trialTypes.cued_reach_1forward==0  & trialTypes.touched_pellet_1forward==1 & (trialTypes.led_1forward==0) & trialTypes.optoGroup~=1 & trialTypes.isLongITI_1forward==1';
    % trial1='trialTypes.cued_reach_1forward==1 & trialTypes.consumed_pellet_1forward==0 & trialTypes.led_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.optoGroup~=1 & trialTypes.isLongITI_1forward==1';
    % trial1='trialTypes.optoGroup~=1 & trialTypes.consumed_pellet_1back==1 & trialTypes.after_cue_success_1forward==1 & trialTypes.consumed_pellet_1forward==1 & trialTypes.led_1forward==1 & trialTypes.optoGroup_1forward~=1';
    test.trial1=trial1;
    test.templateSequence2_cond=eval(trial1);
    trial2='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1';
    % trial2='trialTypes.optoGroup~=1';
    test.trial2=trial2;
    test.templateSequence2_end=eval(trial2);
    test.fillInBetweenWithAnything=true; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
    test.event_name=['alltrials' tbt_filter.name 'inBetweenAnything' num2str(test.fillInBetweenWithAnything)];
    
    saveDir2=[saveDir '\' test.event_name];
    mkdir(saveDir2);
    
    % save settings for paired RT data set
    save([saveDir2 '\test_settings.mat'],'test');
    
    % build paired RT data set
    fakeCueInd=50; % in indices, this is not relevant if not using PCA-based RT model
    skipCorrected=true;
    % this function builds the dataset using the trial type sequences specified above
    [dataset,correctedDistributions]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected);
    
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
%     reachratesettings.acrossSess_window3=[reachratesettings.minTrialLength -0.35];
    reachratesettings.scatterPointSize=50; % size for points in scatter plot
    reachratesettings.addSatietyLines=false; % whether to add proportionality lines to figure
    % reachratesettings.stopPlottingTrialsAfterN=500; % will stop plotting after this nth trial in session, also only use this many trials for regression fit -- see next line, also controls colormap
%     reachratesettings.stopPlottingTrialsAfterN=175; % will stop plotting
    reachratesettings.stopPlottingTrialsAfterN=150; % will stop plotting
    % after this nth trial in session, also only use this many trials for
    % regression fit -- see next line, also controls colormap
    reachratesettings.showFitLine=true; % whether to show linear fit to change across trials
    reachratesettings.useWindowsForUncued=[3]; % to use window2 or window3 or both for the uncued reach rate
    reachratesettings.initWindows=[]; % empty if want to calculate from dataset
    reachratesettings.addSessionLines=false; % for no averaging across sessions plot, whether to connect trial bins within same session with lines
    reachratesettings.binTrialsForAvAcrossSess=true; % whether to bin multiple trials for first figure, will bin into binThisManyTrials
%     reachratesettings.binThisManyTrials=30; % how many trials to bin within each session
    reachratesettings.binThisManyTrials=70; % how many trials to bin within each session
    reachratesettings.nBinsForZones=40; % will be nBinsForZones squared total bins, this is # bins for each x and y axis
    reachratesettings.useRateMethod=3; % 1, 2 or 3 (see explanation below)
    
    reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings);
    close all;
    if isempty(reachrates)
        continue
    end
    if ~isempty(reachratesettings.stopPlottingTrialsAfterN)
        tempuncued=reachrates.uncued;
        tempcued=reachrates.cued;
        tempuncued(reachratesettings.stopPlottingTrialsAfterN+1:end)=nan;
        tempcued(reachratesettings.stopPlottingTrialsAfterN+1:end)=nan;
    end
    daybyday_uncued{i}=downSampAv(tempuncued,reachratesettings.binThisManyTrials);
    daybyday_cued{i}=downSampAv(tempcued,reachratesettings.binThisManyTrials);
    tempuncued=daybyday_uncued{i};
    tempcued=daybyday_cued{i};
    if any(tempuncued(~isnan(tempuncued))>maxuncued)
        maxuncued=nanmax(daybyday_uncued{i});
    end
    if any(tempcued(~isnan(tempcued))>maxcued)
        maxcued=nanmax(daybyday_cued{i});
    end
end

cmap=colormap('cool');
k=1;
kstep=ceil(size(cmap,1)/length(days));
% f=figure('Position',[50 50 1000 1000]);
f=figure();
xlim([0 maxuncued]);
ylim([0 maxcued]);
line([0 maxuncued],[0 maxuncued],'Color',[0.95 0.95 0.95],'LineWidth',3); hold on;
axis square;
% daspect([1 1 1]);
imgcounter=1;
for i=1:length(days)
    withinday_maptocmap=1:ceil(size(cmap,1)/nansum(~isnan(daybyday_uncued{i}))):size(cmap,1);
    tempuncued=daybyday_uncued{i};
    tempcued=daybyday_cued{i};
    if isempty(tempuncued)
        continue
    end
    tempuncued=tempuncued(~isnan(tempuncued));
    tempcued=tempcued(~isnan(tempcued));
    if ~isempty(saveFigDir)
        templine_uncued=daybyday_uncued{i};
        templine_cued=daybyday_cued{i};
        for j=1:length(tempuncued)
            if ~isempty(outlierTest)
                if eval(outlierTest) % e.g., tempcued(j)>1 && i<10
                    continue 
                end
            end
            scatter(tempuncued(j),tempcued(j),50,cmap(withinday_maptocmap(j),:),'filled');
            hold on;
            if j~=1
                plot(templine_uncued(1:j),templine_cued(1:j),'Color',cmap(k,:),'LineWidth',1.5);
            end
            saveas(f,[saveFigDir 'image' sprintf('%03d',imgcounter) imgType]);
            imgcounter=imgcounter+1;
        end
    else
        if ~isempty(outlierTest)
            if eval(outlierTest) % e.g., any(tempcued>1) && i<10
                continue
            end
        end
        scatter(tempuncued,tempcued,50,cmap(withinday_maptocmap(1:nansum(~isnan(daybyday_uncued{i}))),:),'filled');
        hold on;
        plot(daybyday_uncued{i},daybyday_cued{i},'Color',cmap(k,:),'LineWidth',1.5);
    end
    k=k+kstep;
    if k>size(cmap,1)
        k=size(cmap,1);
    end
end
xlabel('Uncued reach rate (1/sec)');
ylabel('Cued reach rate (1/sec)');
end