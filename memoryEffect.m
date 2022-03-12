function memoryEffect(alltbt,metadata,trialTypes,nInSequence,useFractionThroughSession,useReachType,plotCDFUpTo)

flankingTrials='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1'; % take every trial
plotVersusFrac=true;

trialTypes.isLongITI_1forward=[trialTypes.isLongITI(2:end); 0];
trialTypes.optoGroup_1forward=[trialTypes.optoGroup(2:end); 0];
trialTypes.led_4forward=[trialTypes.led(5:end); 0; 0; 0; 0];
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
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir);
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected);
reachratesettings.suppressPlots=false;
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
title('All trial types over session');
[dprimes,fracs_over_sess]=plotDprimesFromReachRates(reachrates,false,plotVersusFrac);
fracsThroughSessBins=[0:0.01:nanmax(fracs_over_sess(1:end)) nanmax(fracs_over_sess(1:end))+0.001];
k=[]; kstep=[]; cmap=[];
for i=1:length(fracsThroughSessBins)-1
    [me_dprime,s_dprime,~,fracs_inThisPart]=getDprimeMeanSE(dprimes,fracs_over_sess,[fracsThroughSessBins(i) fracsThroughSessBins(i+1)]);
    [cmap,k,kstep]=plotMeAndSE_dprime(i,fracsThroughSessBins,me_dprime,s_dprime,fracs_inThisPart,cmap,k,kstep);
end

% THEN
% NO LED FIRST
test.nInSequence=[nInSequence]; % defines trial pairs, e.g., 2 means will compare each trial with its subsequent trial, 3 means will compare each trial with the trial after next, etc.
trial1='trialTypes.led~=1'; 
test.trial1=trial1;
linker=' & trialTypes.led_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.optoGroup_1forward~=3'; 
trial2=['trialTypes.led~=1' linker];
test.trial2=trial2;
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir);
alltbt=useDifferentReachType(alltbt,useReachType,'all_reachBatch','switch');
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
alltbt=useDifferentReachType(alltbt,useReachType,'all_reachBatch','switch back');
plotChangeInReachCDF(dataset.realDistributions,alltbt,plotCDFUpTo); title('No LED');
reachratesettings.suppressPlots=false;
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
[dprimes_noLED,fracs_over_sess_noLED]=plotDprimesFromReachRates(reachrates,false,plotVersusFrac);
title('No led sequence');
[mean_uncued,se_uncued,mean_cued,se_cued,~,~,fracs_inThisPartOfSess_noLED]=getCuedUncuedMeanSE(reachrates,useFractionThroughSession);
[mean_dprime,se_dprime,dprimes_inThisPartOfSess,fracs_inThisPartOfSess_dprimeNoLED]=getDprimeMeanSE(dprimes_noLED,fracs_over_sess_noLED,useFractionThroughSession);

% LED SECOND
test.nInSequence=[nInSequence]; % defines trial pairs, e.g., 2 means will compare each trial with its subsequent trial, 3 means will compare each trial with the trial after next, etc.
linker=' & trialTypes.led_1back~=1'; 
trial1=['trialTypes.led~=1' linker]; 
test.trial1=trial1;
trial2=['trialTypes.led==1 & trialTypes.optoGroup~=1 & trialTypes.optoGroup~=3'];
test.trial2=trial2;
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir);
alltbt=useDifferentReachType(alltbt,useReachType,'all_reachBatch','switch');
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
alltbt=useDifferentReachType(alltbt,useReachType,'all_reachBatch','switch back');
plotChangeInReachCDF(dataset.realDistributions,alltbt,plotCDFUpTo); title('LED');
reachratesettings.suppressPlots=false;
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
[dprimes_LED,fracs_over_sess_LED]=plotDprimesFromReachRates(reachrates,false,plotVersusFrac);
title('Led sequence');
[mean_uncued_LED,se_uncued_LED,mean_cued_LED,se_cued_LED,uncued_inThisPartOfSess,cued_inThisPartOfSess,fracs_inThisPartOfSess_LED]=getCuedUncuedMeanSE(reachrates,useFractionThroughSession);
[mean_dprime_LED,se_dprime_LED,dprimes_inThisPartOfSess,fracs_inThisPartOfSess_dprimeLED]=getDprimeMeanSE(dprimes_LED,fracs_over_sess_LED,useFractionThroughSession);

LED_all_uncued=uncued_inThisPartOfSess(1:end);
LED_all_cued=cued_inThisPartOfSess(1:end);
LED_all_fracs=fracs_inThisPartOfSess_LED(1:end);
LED_all_uncued=LED_all_uncued(~isnan(LED_all_fracs));
LED_all_cued=LED_all_cued(~isnan(LED_all_fracs));
LED_all_fracs=LED_all_fracs(~isnan(LED_all_fracs));
binsForFracs=0:0.05:1;
con=plotHistoOfFracs(binsForFracs,fracs_inThisPartOfSess_noLED,fracs_inThisPartOfSess_LED);
resampleToMatchFracs=false;
if resampleToMatchFracs==true
    [mean_uncued_LED,mean_cued_LED]=resampleMatchingFracs(con,binsForFracs,LED_all_uncued,LED_all_cued,LED_all_fracs);
end

% plot output in 2D cued vs. uncued space
figure();
scatter(mean_uncued,mean_cued,[],'k','filled');
hold on; 
scatter(mean_uncued_LED,mean_cued_LED,[],'r','filled');
line([mean_uncued-se_uncued mean_uncued+se_uncued],[mean_cued mean_cued],'Color','k');
line([mean_uncued mean_uncued],[mean_cued-se_cued mean_cued+se_cued],'Color','k');
line([mean_uncued_LED-se_uncued_LED mean_uncued_LED+se_uncued_LED],[mean_cued_LED mean_cued_LED],'Color','r');
line([mean_uncued_LED mean_uncued_LED],[mean_cued_LED-se_cued_LED mean_cued_LED+se_cued_LED],'Color','r');

% plot output in 1D dprime space
figure();
scatter(nanmean(fracs_inThisPartOfSess_dprimeNoLED(1:end)),mean_dprime,[],'k','filled');
hold on; 
scatter(nanmean(fracs_inThisPartOfSess_dprimeLED(1:end)),mean_dprime_LED,[],'r','filled');
line([nanmean(fracs_inThisPartOfSess_dprimeNoLED(1:end)) nanmean(fracs_inThisPartOfSess_dprimeNoLED(1:end))],[mean_dprime-se_dprime mean_dprime+se_dprime],'Color','k');
line([nanmean(fracs_inThisPartOfSess_dprimeLED(1:end)) nanmean(fracs_inThisPartOfSess_dprimeLED(1:end))],[mean_dprime_LED-se_dprime_LED mean_dprime_LED+se_dprime_LED],'Color','r');

end

function [cmap,k,kstep]=plotMeAndSE_dprime(i,fracsThroughSessBins,me_dprime,s_dprime,fracs_inThisPart,cmap,k,kstep)

if i==1
    figure();
    cmap=colormap('cool');
    hold on;
    k=1;
    kstep=ceil(size(cmap,1)/(length(fracsThroughSessBins)-1));
end
scatter(nanmean(fracs_inThisPart(1:end)),me_dprime,[],cmap(k,:),'filled');    
line([nanmean(fracs_inThisPart(1:end)) nanmean(fracs_inThisPart(1:end))],[me_dprime-s_dprime me_dprime+s_dprime],'Color',cmap(k,:));   
k=k+kstep;
if k>size(cmap,1)
    k=size(cmap,1);
end

end

function con=plotHistoOfFracs(binsForFracs,fracs_inThisPartOfSess_noLED,fracs_inThisPartOfSess_LED)

figure();
[n,x]=histcounts(fracs_inThisPartOfSess_noLED(~isnan(fracs_inThisPartOfSess_noLED)),binsForFracs);
[new_n,new_x]=cityscape_hist(n,x);
plot(new_x,new_n,'Color','k'); hold on;
con.x=x; con.n=n; 
[n,x]=histcounts(fracs_inThisPartOfSess_LED(~isnan(fracs_inThisPartOfSess_LED)),binsForFracs);
[new_n,new_x]=cityscape_hist(n,x);
plot(new_x,new_n,'Color','r');
led.x=x; led.n=n;
title('Before matching fracs');

end

function [mean_uncued_LED,mean_cued_LED]=resampleMatchingFracs(con,binsForFracs,LED_all_uncued,LED_all_cued,LED_all_fracs)

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
    plot(new_x,new_n./nanmax(new_n),'Color','k'); hold on;
    con.x=x; con.n=n;
    [n,x]=histcounts(resampled_LEDfracs(~isnan(resampled_LEDfracs)),binsForFracs);
    [new_n,new_x]=cityscape_hist(n,x);
    plot(new_x,new_n./nanmax(new_n),'Color','r');
    title('After matching fracs');
    mean_uncued_LED=nanmean(resampled_LEDuncued);
    mean_cued_LED=nanmean(resampled_LEDcued);
end

end

function alltbt=useDifferentReachType(alltbt,useReachType,defaultAnalyzedReachType,switchOrSwitchBack)

if isempty(useReachType)
    return
end

switch switchOrSwitchBack
    case 'switch'
        alltbt.(['backup_' defaultAnalyzedReachType])=alltbt.(defaultAnalyzedReachType);
        alltbt.(defaultAnalyzedReachType)=alltbt.(useReachType);
    case 'switch back'
        alltbt.(defaultAnalyzedReachType)=alltbt.(['backup_' defaultAnalyzedReachType]);
end

end

function [mean_dprime,se_dprime,dprimes_inThisPartOfSess,fracs_inThisPartOfSess]=getDprimeMeanSE(dprimes,fracsThroughSess,useFractionThroughSession)

% get average and s.e. reach rate for useFractionThroughSession
dprimes_inThisPartOfSess=dprimes;
fracs_inThisPartOfSess=fracsThroughSess;
for i=1:size(dprimes,1)
    for j=1:size(dprimes,2)
        if fracsThroughSess(i,j)>=useFractionThroughSession(1) && fracsThroughSess(i,j)<useFractionThroughSession(2) 
        else
            dprimes_inThisPartOfSess(i,j)=nan;
            fracs_inThisPartOfSess(i,j)=nan;
        end
    end
end
mean_dprime=nanmean(dprimes_inThisPartOfSess(1:end));
se_dprime=nanstd(dprimes_inThisPartOfSess(1:end),[],2)./sqrt(nansum(~isnan(dprimes_inThisPartOfSess(1:end))));

end

function [mean_uncued,se_uncued,mean_cued,se_cued,uncued_inThisPartOfSess,cued_inThisPartOfSess,fracs_inThisPartOfSess]=getCuedUncuedMeanSE(reachrates,useFractionThroughSession)

% get average and s.e. reach rate for useFractionThroughSession
uncued_inThisPartOfSess=reachrates.alltrials_uncued;
cued_inThisPartOfSess=reachrates.alltrials_cued;
fracs_inThisPartOfSess=reachrates.fracsThroughSess;
for i=1:size(reachrates.alltrials_uncued,1)
    for j=1:size(reachrates.alltrials_uncued,2)
        if reachrates.fracsThroughSess(i,j)>=useFractionThroughSession(1) && reachrates.fracsThroughSess(i,j)<useFractionThroughSession(2) 
        else
            uncued_inThisPartOfSess(i,j)=nan;
            cued_inThisPartOfSess(i,j)=nan;
            fracs_inThisPartOfSess(i,j)=nan;
        end
    end
end
mean_uncued=nanmean(uncued_inThisPartOfSess(1:end));
se_uncued=nanstd(uncued_inThisPartOfSess(1:end),[],2)./sqrt(nansum(~isnan(uncued_inThisPartOfSess(1:end))));
mean_cued=nanmean(cued_inThisPartOfSess(1:end));
se_cued=nanstd(cued_inThisPartOfSess(1:end),[],2)./sqrt(nansum(~isnan(cued_inThisPartOfSess(1:end))));

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
% reachratesettings.acrossSess_window1=[0.05 2]; % cued window [0.05 1]
reachratesettings.acrossSess_window1=[0.05 1]; % cued window [0.05 1]
reachratesettings.acrossSess_window2=[7 reachratesettings.maxTrialLength]; % beware reach suppression after a success
reachratesettings.acrossSess_window3=[reachratesettings.minTrialLength -1]; 
reachratesettings.scatterPointSize=50; % size for points in scatter plot
reachratesettings.addSatietyLines=false; % whether to add proportionality lines to figure
reachratesettings.stopPlottingTrialsAfterN=[];
reachratesettings.showFitLine=false;
reachratesettings.useWindowsForUncued=[3]; % to use window2 or window3 or both for the uncued reach rate
reachratesettings.initWindows=[]; % empty if want to calculate from dataset
reachratesettings.addSessionLines=false; % for no averaging across sessions plot, whether to connect trial bins within same session with lines
reachratesettings.binTrialsForAvAcrossSess=false; % whether to bin multiple trials for first figure, will bin into binThisManyTrials
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
test.fillInBetweenWithAnything=false; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name=['alltrials' tbt_filter.name 'inBetweenAnything' num2str(test.fillInBetweenWithAnything)];
saveDir2=[saveDir '\' test.event_name];
mkdir(saveDir2);
save([saveDir2 '\test_settings.mat'],'test');
fakeCueInd=50; 
skipCorrected=true;

end

function dprimes_lasttrial=getdprimes(reachrates)

% get dprimes per average trial in session
dprimes_lasttrial=[];
if isempty(reachrates)
    return
end
dprimes_lasttrial=calc_dprimes(reachrates.alltrials_uncued,reachrates.alltrials_cued);

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