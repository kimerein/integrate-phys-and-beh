function [dprime_firsthalf,dprime_secondhalf]=getMouseByMouseChangeInDprime(alltbt,metadata,trialTypes,startSessId,endSessId,saveFigDir,imgType,outlierTest,plotDprimes)

days=startSessId:endSessId;
daybyday_uncued=cell(1,length(days));
daybyday_cued=cell(1,length(days));
backup_alltbt=alltbt;
backup_metadata=metadata;
backup_trialTypes=trialTypes;
maxuncued=0;
maxcued=0;
dprime_firsthalf=nan(1,length(days));
dprime_secondhalf=nan(1,length(days));
for i=1:length(days)
    alltbt=backup_alltbt;
    metadata=backup_metadata;
    trialTypes=backup_trialTypes;
    day=days(i);
    disp(['Day ' num2str(i)]);

    % perform any filtering on alltbt
    % for example, filter by d-prime
    
    temp=datestr(datetime('now'));
    temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))='';
    temp=temp(~isspace(temp));
    saveDir=['/Volumes/Neurobio/MICROSCOPE/Kim/RT pairs data sets/' temp]; % where to save details of alltbt filtering and RT pairs data set
    
    % filter settings
    tbt_filter.sortField='mouseid';
%     tbt_filter.sortField='sessid';
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
%     reachratesettings.acrossSess_window3=[reachratesettings.minTrialLength -1];
    reachratesettings.acrossSess_window3=[reachratesettings.minTrialLength -1];
%     reachratesettings.acrossSess_window3=[reachratesettings.minTrialLength -0.35];
    reachratesettings.scatterPointSize=50; % size for points in scatter plot
    reachratesettings.addSatietyLines=false; % whether to add proportionality lines to figure
    % reachratesettings.stopPlottingTrialsAfterN=500; % will stop plotting after this nth trial in session, also only use this many trials for regression fit -- see next line, also controls colormap
%     reachratesettings.stopPlottingTrialsAfterN=190; % will stop plotting
    reachratesettings.stopPlottingTrialsAfterN=200; % will stop plotting
    % after this nth trial in session, also only use this many trials for
    % regression fit -- see next line, also controls colormap
    reachratesettings.showFitLine=true; % whether to show linear fit to change across trials
    reachratesettings.useWindowsForUncued=[3]; % to use window2 or window3 or both for the uncued reach rate
    reachratesettings.initWindows=[]; % empty if want to calculate from dataset
    reachratesettings.addSessionLines=false; % for no averaging across sessions plot, whether to connect trial bins within same session with lines
    reachratesettings.binTrialsForAvAcrossSess=true; % whether to bin multiple trials for first figure, will bin into binThisManyTrials
%     reachratesettings.binThisManyTrials=30; % how many trials to bin within each session
%     reachratesettings.binThisManyTrials=70; % how many trials to bin within each session
    reachratesettings.binThisManyTrials=30; % how many trials to bin within each session
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
    maxReachRate=2; % in reaches per sec
    tempcued(tempcued>maxReachRate)=maxReachRate;
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
    
    plotVersusFrac=false;
    [dprimes,fracs_over_sess]=plotDprimesFromRR(reachrates,true,plotVersusFrac);
    avdprimes=nanmean(dprimes,1);
    halfwaypoint=floor(length(avdprimes)/2);
    dprime_firsthalf(i)=nanmean(avdprimes(1:halfwaypoint));
    dprime_secondhalf(i)=nanmean(avdprimes(halfwaypoint+1:end));
%     thirdwaypoint=floor(length(avdprimes)/3);
%     dprime_firsthalf(i)=nanmean(avdprimes(1:thirdwaypoint));
%     dprime_secondhalf(i)=nanmean(avdprimes(thirdwaypoint+1:2*thirdwaypoint));
end

cmap=colormap('cool');
k=1;
kstep=ceil(size(cmap,1)/(length(days)-1));
% f=figure('Position',[50 50 1000 1000]);
f=figure();
if plotDprimes==false
    xlim([0 maxuncued]);
    ylim([0 maxcued]);
    line([0 maxuncued],[0 maxuncued],'Color',[0.95 0.95 0.95],'LineWidth',3); hold on;
    axis square;
    ylim([0 0.9]);
    xlim([0 0.45]);
    % daspect([1 1 1]);
end
imgcounter=1;
for i=1:length(days)
    if plotDprimes
        withinday_maptocmap=[1 ceil(size(cmap,1)*0.75)];
    else
        withinday_maptocmap=1:ceil(size(cmap,1)/nansum(~isnan(daybyday_uncued{i}))):size(cmap,1);
    end
    tempuncued=daybyday_uncued{i};
    tempcued=daybyday_cued{i};
    if isempty(tempuncued)
        k=k+kstep;
        if k>size(cmap,1)
            k=size(cmap,1);
        end
        continue
    end
%     if ismember(i,[2 4 6 8 10 12 14 16 18 20 22 24 26])
%         k=k+kstep;
%         if k>size(cmap,1)
%             k=size(cmap,1);
%         end
%         continue
%     end
    tempuncued=tempuncued(~isnan(tempuncued));
    tempcued=tempcued(~isnan(tempcued));
    if plotDprimes
        tempuncued=[i-0.3 i+0.3];
        tempcued=[dprime_firsthalf(i) dprime_secondhalf(i)];
    end
    if ~isempty(saveFigDir)
        templine_uncued=daybyday_uncued{i};
        templine_cued=daybyday_cued{i};
        if plotDprimes
            templine_uncued=[i-0.3 i+0.3];
            templine_cued=[dprime_firsthalf(i) dprime_secondhalf(i)];
        end
        for j=1:length(tempuncued)
            if ~isempty(outlierTest) && plotDprimes==false
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
        if ~isempty(outlierTest) && plotDprimes==false
            if eval(outlierTest) % e.g., any(tempcued>1) && i<10
                continue
            end
        end
        if plotDprimes
            scatter([i-0.3 i+0.3]',[dprime_firsthalf(i) dprime_secondhalf(i)]',50,cmap(withinday_maptocmap(1:2),:),'filled');
            hold on;
            plot([i-0.3 i+0.3],[dprime_firsthalf(i) dprime_secondhalf(i)],'Color',cmap(k,:),'LineWidth',1.5);
        else
            scatter(tempuncued,tempcued,50,cmap(withinday_maptocmap(1:nansum(~isnan(daybyday_uncued{i}))),:),'filled');
            hold on;
            plot(daybyday_uncued{i},daybyday_cued{i},'Color',cmap(k,:),'LineWidth',1.5);
        end
    end
    k=k+kstep;
    if k>size(cmap,1)
        k=size(cmap,1);
    end
end
xlabel('Uncued reach rate (1/sec)');
ylabel('Cued reach rate (1/sec)');
end

function [dprimes,fracs_over_sess]=plotDprimesFromRR(reachrates,suppressPlots,plotVersusFrac)

% settings.stopPlottingTrialsAfterN=286;
settings.stopPlottingTrialsAfterN=5000;
settings.binTrialsForAvAcrossSess=true;
settings.binThisManyTrials=1000; %4; %10; % somehow this makes dprime bigger, SO sensitive to this
settings.stopPlottingBinsAfterN=100000;
settings.furtherBinBins=false;
settings.binThisManyBins=2;
settings.plotVersusFrac=plotVersusFrac; % if is true, will plot dprime versus fraction through session instead of trial count

[dprimes,fracs_over_sess]=getdprimes(reachrates,settings.binThisManyTrials,settings);
if suppressPlots==false
    makePlot(settings,dprimes,fracs_over_sess);
    if plotVersusFrac==true
        xlabel('Fraction through session');
    else
        xlabel('Trial #');
    end
    ylabel('d-prime');
end

end

function makePlot(settings,approach_alltrials_dprime,fracs_over_sess)

plotVersusFrac=settings.plotVersusFrac;
figure();
cmap=colormap('cool');
hold on;
k=1;
if ~isempty(settings.stopPlottingBinsAfterN)
    kstep=ceil(size(cmap,1)/settings.stopPlottingBinsAfterN);
else
    kstep=ceil(size(cmap,1)/nansum(~isnan(nanmean(approach_alltrials_dprime,1))));
end
currbincued=[];
currbinfracs=[];
currbintrialnum=0;
currbincounter=0;
trial_nums=1:settings.binThisManyTrials:settings.binThisManyTrials*(size(approach_alltrials_dprime,2)+1);
for i=1:size(approach_alltrials_dprime,2) % across trials
    % rows are different sessions, columns are different trials in each
    % session
    % so ACROSS ALL SESSIONS, take each trial in session
    temp_cued=approach_alltrials_dprime(:,i); % reach rate in trial n+i (last trial) of sequence
    temp_fracs=fracs_over_sess(:,i);
    if i==1 
        if plotVersusFrac==true
            xToPlot=nanmean(temp_fracs);
        else
            xToPlot=trial_nums(i);
        end
        scatter(xToPlot,nanmean(temp_cued),[],'k'); % first trial in SESSION, last trial in sequence
    end
    if ~isempty(settings.stopPlottingBinsAfterN)
        if i<=settings.stopPlottingBinsAfterN
            goAhead=true;
        else
            goAhead=false;
        end
    else
        goAhead=true;
    end
    if goAhead
        if settings.furtherBinBins==true
            currbincued=[currbincued temp_cued];
            currbinfracs=[currbinfracs temp_fracs];
            currbintrialnum=currbintrialnum+trial_nums(i);
            currbincounter=currbincounter+1;
            if currbincounter==settings.binThisManyBins
                % plot and reset
                currbincued=nanmean(currbincued,2);
                currbinfracs=nanmean(currbinfracs,2);
                currbintrialnum=currbintrialnum/settings.binThisManyBins;
                currbincounter=0;
                backup_temp_cued=temp_cued;
                backup_temp_fracs=temp_fracs;
                temp_cued=currbincued;
                temp_fracs=currbinfracs;
                temp_currbintrialnum=currbintrialnum;
                currbincued=[];
                currbinfracs=[];
                currbintrialnum=0;
                if plotVersusFrac==true
                    xToPlot=nanmean(temp_fracs);
                else
                    xToPlot=temp_currbintrialnum;
                end
                scatter(xToPlot,nanmean(temp_cued),[],cmap(k,:),'filled'); % later trials in SESSION, last trial in sequence
                hold on;
                % plot mean and s.e. across session
                line([xToPlot xToPlot],...
                     [nanmean(temp_cued)-nanstd(temp_cued)./sqrt(nansum(~isnan(temp_cued))) nanmean(temp_cued)+nanstd(temp_cued)./sqrt(nansum(~isnan(temp_cued)))],'Color',cmap(k,:),'LineWidth',1);
                temp_cued=backup_temp_cued;
                temp_fracs=backup_temp_fracs;
            end
        else
            if plotVersusFrac==true
                xToPlot=nanmean(temp_fracs);
            else
                xToPlot=trial_nums(i);
            end
            scatter(xToPlot,nanmean(temp_cued),[],cmap(k,:),'filled'); % later trials in SESSION, last trial in sequence
            hold on;
            line([xToPlot xToPlot],...
                [nanmean(temp_cued)-nanstd(temp_cued)./sqrt(nansum(~isnan(temp_cued))) nanmean(temp_cued)+nanstd(temp_cued)./sqrt(nansum(~isnan(temp_cued)))],'Color',cmap(k,:),'LineWidth',1);
        end
    end
    k=k+kstep;
    if k>size(cmap,1)
        k=size(cmap,1);
    end
end

end

function [dprimes_over_sess,fracs_over_sess]=getdprimes(reachrates,trialBinSize,settings)

% get dprimes per average trial in session
if isempty(reachrates)
    return
end
dprimes_over_sess=nan(size(reachrates.alltrials_uncued,1),ceil(size(reachrates.alltrials_uncued,2)/trialBinSize));
fracs_over_sess=nan(size(reachrates.alltrials_uncued,1),ceil(size(reachrates.alltrials_uncued,2)/trialBinSize));
currbincounter=0;
currtrial=0;
goAhead=true;
for i=1:size(reachrates.alltrials_uncued,2)
    if goAhead==false
        break
    end
    currbincounter=currbincounter+1;
    takeInds=i:i+trialBinSize-1;
    if any(takeInds>size(reachrates.alltrials_uncued,2))
        takeInds=i:size(reachrates.alltrials_uncued,2);
    end
    currtrial=currtrial+length(takeInds);
    if ~isempty(settings.stopPlottingTrialsAfterN)
        if currtrial>settings.stopPlottingTrialsAfterN
            % stop after this bin
            goAhead=false;
        end
    end            
    currbincued=reachrates.alltrials_cued(:,takeInds); % reach rate in trial n+i (last trial) of sequence
    currbinuncued=reachrates.alltrials_uncued(:,takeInds); % reach rate in trial n+i (last trial) of sequence
    curr_bin_fracs=reachrates.fracsThroughSess(:,takeInds);
    if all(isnan(currbincued(1:end))) && all(isnan(currbinuncued(1:end)))
        currbincounter=currbincounter-1;
        break
    end
    fracs_av=nanmean(curr_bin_fracs,2);
    dprimes_lasttrial=calc_dprimes(currbinuncued,currbincued);
    dprimes_over_sess(:,currbincounter)=dprimes_lasttrial;
    fracs_over_sess(:,currbincounter)=fracs_av;
end
dprimes_over_sess=dprimes_over_sess(:,1:currbincounter);
fracs_over_sess=fracs_over_sess(:,1:currbincounter);
fi=find(all(isnan(dprimes_over_sess),1),2,'first');
if ~isempty(fi)
    dprimes_over_sess=dprimes_over_sess(:,1:fi(1)-1);
    fracs_over_sess=fracs_over_sess(:,1:fi(1)-1);
end

end

function dprimes=calc_dprimes(uncued_events,cued_events)

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
