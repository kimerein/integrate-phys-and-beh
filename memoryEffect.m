function memoryEffect(alltbt,metadata,trialTypes,nInSequence,useFractionThroughSession,useReachType,plotCDFUpTo,withinCueTimeWindow)

flankingTrials='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1'; % take every trial
%flankingTrials='trialTypes.led~=1'; % take no LED trials only
plotVersusFrac=false;

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
reachratesettings=getReachRateSettings(withinCueTimeWindow);

% FIRST PLOT CHANGE OVER ALL TRIALS WITHIN SESSION
trial1=flankingTrials; 
test.trial1=trial1;
trial2=flankingTrials;
test.trial2=trial2;
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir);
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected);
% playing with an alternate metric: given reach, probability that cue preceded reach
% fracthrubins=0:0.1:1.0001;
fracthrubins={[0:0.05:0.95],[0.1:0.05:1.001]};
[dprime_given_reach,dprime_given_reachPLUSsd,dprime_given_reachMINUSsd]=dprimes_given_reach(alltbt,dataset,false,fracthrubins,'rawReaching_event_trialiInSeq',withinCueTimeWindow); % two formats for fracthrubins, see function
% continuing to plot change over all trials within session
reachratesettings.suppressPlots=false;
reachratesettings.binTrialsForAvAcrossSess=true;
reachratesettings.binThisManyTrials=20; 
reachratesettings.stopPlottingTrialsAfterN=200; %120;
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
title('All trial types over session');
[dprimes,fracs_over_sess,initialDprimes]=plotDprimesFromReachRates(reachrates,false,plotVersusFrac);
% fracsThroughSessBins=[0:0.01:nanmax(fracs_over_sess(1:end)) nanmax(fracs_over_sess(1:end))+0.001];
fracsThroughSessBins=[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
k=[]; kstep=[]; cmap=[];
for i=1:length(fracsThroughSessBins)-1
    [me_dprime,s_dprime,~,fracs_inThisPart]=getDprimeMeanSE(dprimes,fracs_over_sess,[fracsThroughSessBins(i) fracsThroughSessBins(i+1)]);
    [cmap,k,kstep]=plotMeAndSE_dprime(i,fracsThroughSessBins,me_dprime,s_dprime,fracs_inThisPart,cmap,k,kstep);
    xlabel('Fraction through session');
end
% Plot 2D sorted by fracThroughSession
fracsThroughSessBins=[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7];
k=[]; kstep=[]; cmap=[];
for i=1:length(fracsThroughSessBins)-1
    [me_uncued,s_uncued,me_cued,s_cued]=getCuedUncuedMeanSE(reachrates,[fracsThroughSessBins(i) fracsThroughSessBins(i+1)]);
    [cmap,k,kstep]=plotMeAndSE_2D(i,fracsThroughSessBins,me_uncued,s_uncued,me_cued,s_cued,cmap,k,kstep);
end

reachratesettings=getReachRateSettings(withinCueTimeWindow);
% THEN
% NO LED FIRST
test.nInSequence=[nInSequence]; % defines trial pairs, e.g., 2 means will compare each trial with its subsequent trial, 3 means will compare each trial with the trial after next, etc.
% trial1='trialTypes.led~=1';
trial1=flankingTrials;
test.trial1=trial1;
linker=' & trialTypes.led_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.optoGroup_1forward~=3'; 
% trial2=['trialTypes.led~=1' linker];
trial2=['trialTypes.led~=1'];
test.trial2=trial2;
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir);
alltbt=useDifferentReachType(alltbt,useReachType,'all_reachBatch','switch');
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
alltbt=useDifferentReachType(alltbt,useReachType,'all_reachBatch','switch back');
fracthrubins={[0:0.05:0.95],[0.1:0.05:1.001]};
dprimes_given_reach(alltbt,dataset,false,fracthrubins,'rawReaching_event_trialiInSeq',withinCueTimeWindow);
fracthrubins=0:useFractionThroughSession(1):useFractionThroughSession(2)+0.001;
[dprime_given_reach_noLED,dprime_given_reachPLUSsd_noLED,dprime_given_reachMINUSsd_noLED]=dprimes_given_reach(alltbt,dataset,false,fracthrubins,'rawReaching_event_trialiInSeq',withinCueTimeWindow);
title('dprime NO LED given reach');
returnThis=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rawReaching'); title('No LED');
plotChangeInReachCDF(dataset.realDistributions,alltbt,plotCDFUpTo); title('No LED');
reachratesettings.suppressPlots=false;
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
[dprimes_noLED,fracs_over_sess_noLED]=plotDprimesFromReachRates(reachrates,false,plotVersusFrac,initialDprimes);
title('No led sequence');
[mean_uncued,se_uncued,mean_cued,se_cued,uncued_inThisPartOfSess_noLED,cued_inThisPartOfSess_noLED,fracs_inThisPartOfSess_noLED]=getCuedUncuedMeanSE(reachrates,useFractionThroughSession);
[mean_dprime,se_dprime,dprimes_inThisPartOfSess,fracs_inThisPartOfSess_dprimeNoLED]=getDprimeMeanSE(dprimes_noLED,fracs_over_sess_noLED,useFractionThroughSession);

% LED SECOND
test.nInSequence=[nInSequence]; % defines trial pairs, e.g., 2 means will compare each trial with its subsequent trial, 3 means will compare each trial with the trial after next, etc.
linker=' & trialTypes.led_1back~=1'; 
% trial1=['trialTypes.led~=1' linker]; 
% trial1=['trialTypes.led~=1'];
trial1=flankingTrials;
test.trial1=trial1;
% trial2=['trialTypes.led==1 & trialTypes.optoGroup~=1 & trialTypes.optoGroup~=3'];
trial2=['trialTypes.optoGroup==2'];
test.trial2=trial2;
[test,fakeCueInd,skipCorrected]=fillInRestOfTest(nInSequence,trial1,trial2,trialTypes,saveDir);
alltbt=useDifferentReachType(alltbt,useReachType,'all_reachBatch','switch');
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
alltbt=useDifferentReachType(alltbt,useReachType,'all_reachBatch','switch back');
fracthrubins={[0:0.05:0.95],[0.1:0.05:1.001]};
dprimes_given_reach(alltbt,dataset,false,fracthrubins,'rawReaching_event_trialiInSeq',withinCueTimeWindow);
fracthrubins=0:useFractionThroughSession(1):useFractionThroughSession(2)+0.001;
[dprime_given_reach_LED,dprime_given_reachPLUSsd_LED,dprime_given_reachMINUSsd_LED]=dprimes_given_reach(alltbt,dataset,false,fracthrubins,'rawReaching_event_trialiInSeq',withinCueTimeWindow);
title('dprime LED given reach');
returnThis=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rawReaching'); title('LED');
plotChangeInReachCDF(dataset.realDistributions,alltbt,plotCDFUpTo); title('LED');
reachratesettings.suppressPlots=false;
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
[dprimes_LED,fracs_over_sess_LED]=plotDprimesFromReachRates(reachrates,false,plotVersusFrac,initialDprimes);
title('Led sequence');
[mean_uncued_LED,se_uncued_LED,mean_cued_LED,se_cued_LED,uncued_inThisPartOfSess_LED,cued_inThisPartOfSess_LED,fracs_inThisPartOfSess_LED]=getCuedUncuedMeanSE(reachrates,useFractionThroughSession);
[mean_dprime_LED,se_dprime_LED,dprimes_inThisPartOfSess,fracs_inThisPartOfSess_dprimeLED]=getDprimeMeanSE(dprimes_LED,fracs_over_sess_LED,useFractionThroughSession);

noLED_all_uncued=uncued_inThisPartOfSess_noLED(1:end);
noLED_all_cued=cued_inThisPartOfSess_noLED(1:end);
noLED_all_fracs=fracs_inThisPartOfSess_noLED(1:end);
noLED_all_uncued=noLED_all_uncued(~isnan(noLED_all_fracs));
noLED_all_cued=noLED_all_cued(~isnan(noLED_all_fracs));
noLED_all_fracs=noLED_all_fracs(~isnan(noLED_all_fracs));

LED_all_uncued=uncued_inThisPartOfSess_LED(1:end);
LED_all_cued=cued_inThisPartOfSess_LED(1:end);
LED_all_fracs=fracs_inThisPartOfSess_LED(1:end);
LED_all_uncued=LED_all_uncued(~isnan(LED_all_fracs));
LED_all_cued=LED_all_cued(~isnan(LED_all_fracs));
LED_all_fracs=LED_all_fracs(~isnan(LED_all_fracs));
binsForFracs=0:0.01:1;
con=plotHistoOfFracs(binsForFracs,fracs_inThisPartOfSess_noLED,fracs_inThisPartOfSess_LED);
% what's more important is to the sample the actual number of trials,
% because it's binomial
resampleToMatchNumberTrials=true; %false; 
if resampleToMatchNumberTrials
    % whichever has fewer trials
    nTrials_con=nansum(~isnan(noLED_all_uncued));
    nTrials_led=nansum(~isnan(LED_all_uncued));
    if nTrials_con<nTrials_led
        takeN=nTrials_con;
    else
        takeN=nTrials_led;
    end
    [mean_uncued,mean_cued]=resampleMatchNumTrials(noLED_all_uncued,noLED_all_cued,noLED_all_fracs,takeN);
    [mean_uncued_LED,mean_cued_LED]=resampleMatchNumTrials(LED_all_uncued,LED_all_cued,LED_all_fracs,takeN);    
end
resampleToMatchFracs=false;
if resampleToMatchFracs
    [mean_uncued_LED,mean_cued_LED]=resampleMatchingFracs(con,binsForFracs,LED_all_uncued,LED_all_cued,LED_all_fracs,fracs_inThisPartOfSess_noLED);
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
title('dprime useFractionThroughSession');

% plot output in dprime GIVEN REACH space
figure();
scatter(nanmean(useFractionThroughSession),dprime_given_reach_noLED(end),[],'k','filled');
hold on; 
scatter(nanmean(useFractionThroughSession),dprime_given_reach_LED(end),[],'r','filled');
line([nanmean(useFractionThroughSession) nanmean(useFractionThroughSession)],[dprime_given_reachMINUSsd_noLED(end) dprime_given_reachPLUSsd_noLED(end)],'Color','k');
line([nanmean(useFractionThroughSession) nanmean(useFractionThroughSession)],[dprime_given_reachMINUSsd_LED(end) dprime_given_reachPLUSsd_LED(end)],'Color','r');
title('dprime GIVEN REACH');


end

function [dprime_given_reach,dprime_given_reachPLUSsd,dprime_given_reachMINUSsd,p_reach_preceded_by_cue_tbt,sample_std_for_binomial_tbt,p_reach_followed_by_cue_tbt,sample_std_for_binomial_follow_tbt]=dprimes_given_reach(alltbt,dataset,suppressPlots,fracthrubins,whichEvent,withinCueTimeWindow)

if ~iscell(fracthrubins)
    fracthrubins_ends=[];
else
    temp=fracthrubins;
    fracthrubins=[temp{1} 1];
    fracthrubins_ends=[temp{2} 1];
end

% playing with an alternate metric: given reach, probability that cue preceded reach
temp=nanmean(alltbt.times,1); timestep=mode(diff(temp(~isnan(temp)))); [~,cueAtInd]=nanmax(nanmean(alltbt.cueZone_onVoff,1));
[was_reach_preceded_by_cue,trial_number,reach_time]=givenReach_probThatWasPrecededByCue(dataset.realDistributions, whichEvent, withinCueTimeWindow, timestep, cueAtInd); % if third arg is negative, will get prob cue AFTER reach
fracthru_allevs=alltbt.fractionThroughSess_adjusted(dataset.realDistributions.event_isSeq{1}==1);
fracthru=fracthru_allevs(trial_number);
% fracthrubins=0:0.1:1.0001;
p_reach_preceded_by_cue=nan(length(fracthrubins)-1,1);
sample_std_for_binomial=nan(length(fracthrubins)-1,1);
ns_reach_preceded_by_cue=nan(length(fracthrubins)-1,1);
% trial by trial
tbtu=unique(trial_number);
p_reach_preceded_by_cue_tbt=nan(length(tbtu),1);
sample_std_for_binomial_tbt=nan(length(tbtu),1);
for i=1:length(tbtu)
    p_reach_preceded_by_cue_tbt(i)=nanmean(was_reach_preceded_by_cue(trial_number==tbtu(i)));
    sample_std_for_binomial_tbt(i)=sqrt((p_reach_preceded_by_cue_tbt(i)*(1-p_reach_preceded_by_cue_tbt(i)))/nansum(trial_number==tbtu(i)));
end
% frac through sess
for i=1:length(fracthrubins)-1
    if isempty(fracthrubins_ends)
        p_reach_preceded_by_cue(i)=nanmean(was_reach_preceded_by_cue(fracthru>=fracthrubins(i) & fracthru<fracthrubins(i+1)));
        sample_std_for_binomial(i)=sqrt((p_reach_preceded_by_cue(i)*(1-p_reach_preceded_by_cue(i)))/nansum(fracthru>=fracthrubins(i) & fracthru<fracthrubins(i+1)));
        ns_reach_preceded_by_cue(i)=nansum(fracthru>=fracthrubins(i) & fracthru<fracthrubins(i+1));
    else
        p_reach_preceded_by_cue(i)=nanmean(was_reach_preceded_by_cue(fracthru>=fracthrubins(i) & fracthru<fracthrubins_ends(i)));
        sample_std_for_binomial(i)=sqrt((p_reach_preceded_by_cue(i)*(1-p_reach_preceded_by_cue(i)))/nansum(fracthru>=fracthrubins(i) & fracthru<fracthrubins_ends(i)));
        ns_reach_preceded_by_cue(i)=nansum(fracthru>=fracthrubins(i) & fracthru<fracthrubins_ends(i));
    end
end
if suppressPlots~=true
    figure(); plot(fracthrubins(1:end-1),p_reach_preceded_by_cue); hold on;
    for i=1:length(fracthrubins)-1
        line([fracthrubins(i) fracthrubins(i)],[p_reach_preceded_by_cue(i)-sample_std_for_binomial(i) p_reach_preceded_by_cue(i)+sample_std_for_binomial(i)]);
    end
    xlabel('Fraction through session'); ylabel('P (m sd given binom) reach preceded by cue');

    ds=10;
    figure(); plot(downSampAv(p_reach_preceded_by_cue_tbt,ds)); hold on;
    xlabel('Bin'); ylabel('P reach preceded by cue TRIALS BINNED');
end
[was_reach_followed_by_cue,trial_number,reach_time]=givenReach_probThatWasPrecededByCue(dataset.realDistributions, whichEvent, -1.5, timestep, cueAtInd); % if third arg is negative, will get prob cue AFTER reach
fracthru_allevs=alltbt.fractionThroughSess_adjusted(dataset.realDistributions.event_isSeq{1}==1);
fracthru=fracthru_allevs(trial_number);
p_reach_followed_by_cue=nan(length(fracthrubins)-1,1);
sample_std_for_binomial_follow=nan(length(fracthrubins)-1,1);
ns_reach_followed_by_cue=nan(length(fracthrubins)-1,1);
% trial by trial
tbtu=unique(trial_number);
p_reach_followed_by_cue_tbt=nan(length(tbtu),1);
sample_std_for_binomial_follow_tbt=nan(length(tbtu),1);
for i=1:length(tbtu)
    p_reach_followed_by_cue_tbt(i)=nanmean(was_reach_followed_by_cue(trial_number==tbtu(i)));
    sample_std_for_binomial_follow_tbt(i)=sqrt((p_reach_followed_by_cue_tbt(i)*(1-p_reach_followed_by_cue_tbt(i)))/nansum(trial_number==tbtu(i)));
end
% frac thru
for i=1:length(fracthrubins)-1
    if isempty(fracthrubins_ends)
        p_reach_followed_by_cue(i)=nanmean(was_reach_followed_by_cue(fracthru>=fracthrubins(i) & fracthru<fracthrubins(i+1)));
        sample_std_for_binomial_follow(i)=sqrt((p_reach_followed_by_cue(i)*(1-p_reach_followed_by_cue(i)))/nansum(fracthru>=fracthrubins(i) & fracthru<fracthrubins(i+1)));
        ns_reach_followed_by_cue(i)=nansum(fracthru>=fracthrubins(i) & fracthru<fracthrubins(i+1));
    else
        p_reach_followed_by_cue(i)=nanmean(was_reach_followed_by_cue(fracthru>=fracthrubins(i) & fracthru<fracthrubins_ends(i)));
        sample_std_for_binomial_follow(i)=sqrt((p_reach_followed_by_cue(i)*(1-p_reach_followed_by_cue(i)))/nansum(fracthru>=fracthrubins(i) & fracthru<fracthrubins_ends(i)));
        ns_reach_followed_by_cue(i)=nansum(fracthru>=fracthrubins(i) & fracthru<fracthrubins_ends(i));
    end
end
if suppressPlots~=true
    figure(); plot(fracthrubins(1:end-1),p_reach_followed_by_cue); hold on;
    for i=1:length(fracthrubins)-1
        line([fracthrubins(i) fracthrubins(i)],[p_reach_followed_by_cue(i)-sample_std_for_binomial_follow(i) p_reach_followed_by_cue(i)+sample_std_for_binomial_follow(i)]);
    end
    xlabel('Fraction through session'); ylabel('P (m sd given binom) reach followed by cue');

    ds=10;
    figure(); plot(downSampAv(p_reach_followed_by_cue_tbt,ds));
    xlabel('Bin'); ylabel('P reach followed by cue TRIALS BINNED');
end
p_reach_preceded_by_cue(p_reach_preceded_by_cue==0)=1./(2.*ns_reach_preceded_by_cue(p_reach_preceded_by_cue==0));
p_reach_preceded_by_cue(p_reach_preceded_by_cue==1)=1-(1./(2.*ns_reach_preceded_by_cue(p_reach_preceded_by_cue==1)));
p_reach_followed_by_cue(p_reach_followed_by_cue==0)=1./(2.*ns_reach_followed_by_cue(p_reach_followed_by_cue==0));
p_reach_followed_by_cue(p_reach_followed_by_cue==1)=1-(1./(2.*ns_reach_followed_by_cue(p_reach_followed_by_cue==1)));
dprime_given_reach=norminv(p_reach_preceded_by_cue)-norminv(p_reach_followed_by_cue);
dprime_given_reachPLUSsd=norminv(p_reach_preceded_by_cue+sample_std_for_binomial)-norminv(p_reach_followed_by_cue-sample_std_for_binomial_follow);
dprime_given_reachMINUSsd=norminv(p_reach_preceded_by_cue-sample_std_for_binomial)-norminv(p_reach_followed_by_cue+sample_std_for_binomial_follow);
if suppressPlots~=true
    figure(); plot(fracthrubins(1:end-1),dprime_given_reach);
    for i=1:length(fracthrubins)-1
        line([fracthrubins(i) fracthrubins(i)],[dprime_given_reachMINUSsd(i) dprime_given_reachPLUSsd(i)]);
    end
    xlabel('Fraction through session'); ylabel('Dprime given reach');
end

end

function [cmap,k,kstep]=plotMeAndSE_2D(i,fracsThroughSessBins,me_uncued,s_uncued,me_cued,s_cued,cmap,k,kstep)

if i==1
    figure();
    cmap=colormap('cool');
    hold on;
    k=1;
    kstep=ceil(size(cmap,1)/(length(fracsThroughSessBins)-1));
end
scatter(me_uncued,me_cued,[],cmap(k,:),'filled');  
line([me_uncued-s_uncued me_uncued+s_uncued],[me_cued me_cued],'Color',cmap(k,:));   
line([me_uncued me_uncued],[me_cued-s_cued me_cued+s_cued],'Color',cmap(k,:));   
k=k+kstep;
if k>size(cmap,1)
    k=size(cmap,1);
end

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

function [mean_uncued_LED,mean_cued_LED,mean_fracs_LED]=resampleMatchNumTrials(LED_all_uncued,LED_all_cued,LED_all_fracs,takeN)

nRuns=10;
resampled_LEDuncued=[];
resampled_LEDcued=[];
resampled_LEDfracs=[];
newpool_LEDuncued=LED_all_uncued;
newpool_LEDcued=LED_all_cued;
newpool_LEDfracs=LED_all_fracs;
for j=1:nRuns
    r=randsample(length(newpool_LEDuncued),takeN);
    resampled_LEDuncued=[resampled_LEDuncued newpool_LEDuncued(r)];
    resampled_LEDcued=[resampled_LEDcued newpool_LEDcued(r)];
    resampled_LEDfracs=[resampled_LEDfracs newpool_LEDfracs(r)];
end

mean_uncued_LED=nanmean(resampled_LEDuncued);
mean_cued_LED=nanmean(resampled_LEDcued);
mean_fracs_LED=nanmean(resampled_LEDfracs);

end

function [mean_uncued_LED,mean_cued_LED]=resampleMatchingFracs(con,binsForFracs,LED_all_uncued,LED_all_cued,LED_all_fracs,fracs_inThisPartOfSess_noLED)

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
% mean_uncued=nanmean(nanmean(uncued_inThisPartOfSess,1),2);
mean_uncued=nanmean(uncued_inThisPartOfSess(1:end));
se_uncued=nanstd(uncued_inThisPartOfSess(1:end),[],2)./sqrt(nansum(~isnan(uncued_inThisPartOfSess(1:end))));
% mean_cued=nanmean(nanmean(cued_inThisPartOfSess,1),2);
mean_cued=nanmean(cued_inThisPartOfSess(1:end));
se_cued=nanstd(cued_inThisPartOfSess(1:end),[],2)./sqrt(nansum(~isnan(cued_inThisPartOfSess(1:end))));

end

function reachratesettings=getReachRateSettings(withinCueTimeWindow)

reachratesettings.epsilon_cue=0; % in seconds
reachratesettings.epsilon_uncue=2; % in seconds
reachratesettings.epsilon_beforecue=1; % in seconds
reachratesettings.percentOfReachesFromSess_forInitCond=20; % use this fraction of trials from beginning of session to get initial windows
reachratesettings.percentOfReachesFromSess_forInitRate=20; % use this fraction of trials from beginning of session to get initial reach rates
reachratesettings.maxTrialLength=9.5; % in sec, wrt cue
reachratesettings.minTrialLength=-2; % wrt cue, in sec
reachratesettings.suppressPlots=true;
% reachratesettings.acrossSess_window1=[0.05 2]; % cued window [0.05 1]
if isempty(withinCueTimeWindow)
    reachratesettings.acrossSess_window1=[0.05 1.05]; % cued window [0.05 1]
else
    reachratesettings.acrossSess_window1=[0 0+withinCueTimeWindow]; 
end
reachratesettings.acrossSess_window2=[7 reachratesettings.maxTrialLength]; % beware reach suppression after a success
if isempty(withinCueTimeWindow)
    reachratesettings.acrossSess_window3=[reachratesettings.minTrialLength -1]; 
else
    reachratesettings.acrossSess_window3=[0-withinCueTimeWindow 0]; 
end
reachratesettings.scatterPointSize=50; % size for points in scatter plot
reachratesettings.addSatietyLines=false; % whether to add proportionality lines to figure
reachratesettings.stopPlottingTrialsAfterN=286;
reachratesettings.showFitLine=false;
%reachratesettings.discardNoisySessions=false; % will discard sessions with fewer than 3 trials
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