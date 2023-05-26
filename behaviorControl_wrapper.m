function [alltog,postOutcome_reaches,chewDurations]=behaviorControl_wrapper(data_loc_array)

cueOffset=-0.16; % match what used for physiology
behReadoutTimeWindow=[0 5]; % in sec from alignCompanion
maxTrialsPerSess=250; % cannot be more than this many trials of any type

% Cued success
event='success_fromPerchOrWheel';
timeWindow=[0+cueOffset 3]; % in seconds from cue onset
[cueSucc_chewendings,cueSucc_postoutcome_reaches,cueSucc_fromwhichsess_reaches,cueSucc_fromwhichsess_chews]=getChewsAndReachFromAllSess(data_loc_array,maxTrialsPerSess,event,cueOffset,timeWindow,behReadoutTimeWindow);
alltog.cueSucc_chewendings=cueSucc_chewendings;
alltog.cueSucc_postoutcome_reaches=cueSucc_postoutcome_reaches;
alltog.cueSucc_fromwhichsess_reaches=cueSucc_fromwhichsess_reaches;
alltog.cueSucc_fromwhichsess_chews=cueSucc_fromwhichsess_chews;

% Cued failure and no reaching afterward
event='failure_noSuccessBeforeAndNoReachingAfter';
timeWindow=[0+cueOffset 3]; % in seconds from cue onset
[cueFail_chewendings,cueFail_postoutcome_reaches,cueFail_fromwhichsess_reaches,cueFail_fromwhichsess_chews]=getChewsAndReachFromAllSess(data_loc_array,maxTrialsPerSess,event,cueOffset,timeWindow,behReadoutTimeWindow);
alltog.cueFail_chewendings=cueFail_chewendings;
alltog.cueFail_postoutcome_reaches=cueFail_postoutcome_reaches;
alltog.cueFail_fromwhichsess_reaches=cueFail_fromwhichsess_reaches;
alltog.cueFail_fromwhichsess_chews=cueFail_fromwhichsess_chews;

% Uncued success
event='success_fromPerchOrWheel';
timeWindow=[3 16]; % in seconds from cue onset
[uncueSucc_chewendings,uncueSucc_postoutcome_reaches,uncueSucc_fromwhichsess_reaches,uncueSucc_fromwhichsess_chews]=getChewsAndReachFromAllSess(data_loc_array,maxTrialsPerSess,event,cueOffset,timeWindow,behReadoutTimeWindow);
alltog.uncueSucc_chewendings=uncueSucc_chewendings;
alltog.uncueSucc_postoutcome_reaches=uncueSucc_postoutcome_reaches;
alltog.uncueSucc_fromwhichsess_reaches=uncueSucc_fromwhichsess_reaches;
alltog.uncueSucc_fromwhichsess_chews=uncueSucc_fromwhichsess_chews;

% Uncued failure and no reaching afterward
event='failure_noSuccessBeforeAndNoReachingAfter';
timeWindow=[3 16]; % in seconds from cue onset
[uncueFail_chewendings,uncueFail_postoutcome_reaches,uncueFail_fromwhichsess_reaches,uncueFail_fromwhichsess_chews]=getChewsAndReachFromAllSess(data_loc_array,maxTrialsPerSess,event,cueOffset,timeWindow,behReadoutTimeWindow);
alltog.uncueFail_chewendings=uncueFail_chewendings;
alltog.uncueFail_postoutcome_reaches=uncueFail_postoutcome_reaches;
alltog.uncueFail_fromwhichsess_reaches=uncueFail_fromwhichsess_reaches;
alltog.uncueFail_fromwhichsess_chews=uncueFail_fromwhichsess_chews;

usess=unique([cueSucc_fromwhichsess_reaches; cueFail_fromwhichsess_reaches; uncueSucc_fromwhichsess_reaches; uncueFail_fromwhichsess_reaches]);
dp_per_sess_cuedSucc_v_cuedFail=nan(length(usess),1);
dp_per_sess_uncuedSucc_v_uncuedFail=nan(length(usess),1);
dp_per_sess_cuedSucc_v_uncuedSucc=nan(length(usess),1);
dp_per_sess_cuedFail_v_uncuedFail=nan(length(usess),1);
for i=1:length(usess)
    currsess=usess(i);
    dp_per_sess_cuedSucc_v_cuedFail(i)=rms_dprime(cueSucc_postoutcome_reaches(cueSucc_fromwhichsess_reaches==currsess),cueFail_postoutcome_reaches(cueFail_fromwhichsess_reaches==currsess));
    dp_per_sess_uncuedSucc_v_uncuedFail(i)=rms_dprime(uncueSucc_postoutcome_reaches(uncueSucc_fromwhichsess_reaches==currsess),uncueFail_postoutcome_reaches(uncueFail_fromwhichsess_reaches==currsess));
    dp_per_sess_cuedSucc_v_uncuedSucc(i)=rms_dprime(cueSucc_postoutcome_reaches(cueSucc_fromwhichsess_reaches==currsess),uncueSucc_postoutcome_reaches(uncueSucc_fromwhichsess_reaches==currsess));
    dp_per_sess_cuedFail_v_uncuedFail(i)=rms_dprime(cueFail_postoutcome_reaches(cueFail_fromwhichsess_reaches==currsess),uncueFail_postoutcome_reaches(uncueFail_fromwhichsess_reaches==currsess));
end
postOutcome_reaches.dp_per_sess_cuedSucc_v_cuedFail=dp_per_sess_cuedSucc_v_cuedFail;
postOutcome_reaches.dp_per_sess_uncuedSucc_v_uncuedFail=dp_per_sess_uncuedSucc_v_uncuedFail;
postOutcome_reaches.dp_per_sess_cuedSucc_v_uncuedSucc=dp_per_sess_cuedSucc_v_uncuedSucc;
postOutcome_reaches.dp_per_sess_cuedFail_v_uncuedFail=dp_per_sess_cuedFail_v_uncuedFail;

usess=unique([cueSucc_fromwhichsess_chews; cueFail_fromwhichsess_chews; uncueSucc_fromwhichsess_chews; uncueFail_fromwhichsess_chews]);
dp_per_sess_cuedSucc_v_cuedFail=nan(length(usess),1);
dp_per_sess_uncuedSucc_v_uncuedFail=nan(length(usess),1);
dp_per_sess_cuedSucc_v_uncuedSucc=nan(length(usess),1);
dp_per_sess_cuedFail_v_uncuedFail=nan(length(usess),1);
for i=1:length(usess)
    currsess=usess(i);
    dp_per_sess_cuedSucc_v_cuedFail(i)=rms_dprime(cueSucc_chewendings(cueSucc_fromwhichsess_chews==currsess),cueFail_chewendings(cueFail_fromwhichsess_chews==currsess));
    dp_per_sess_uncuedSucc_v_uncuedFail(i)=rms_dprime(uncueSucc_chewendings(uncueSucc_fromwhichsess_chews==currsess),uncueFail_chewendings(uncueFail_fromwhichsess_chews==currsess));
    dp_per_sess_cuedSucc_v_uncuedSucc(i)=rms_dprime(cueSucc_chewendings(cueSucc_fromwhichsess_chews==currsess),uncueSucc_chewendings(uncueSucc_fromwhichsess_chews==currsess));
    dp_per_sess_cuedFail_v_uncuedFail(i)=rms_dprime(cueFail_chewendings(cueFail_fromwhichsess_chews==currsess),uncueFail_chewendings(uncueFail_fromwhichsess_chews==currsess));
end
chewDurations.dp_per_sess_cuedSucc_v_cuedFail=dp_per_sess_cuedSucc_v_cuedFail;
chewDurations.dp_per_sess_uncuedSucc_v_uncuedFail=dp_per_sess_uncuedSucc_v_uncuedFail;
chewDurations.dp_per_sess_cuedSucc_v_uncuedSucc=dp_per_sess_cuedSucc_v_uncuedSucc;
chewDurations.dp_per_sess_cuedFail_v_uncuedFail=dp_per_sess_cuedFail_v_uncuedFail;

% Some plots
figure();
plotSingleHist(alltog.cueSucc_chewendings,[0:0.1:5],'g'); hold on;
plotSingleHist(alltog.cueFail_chewendings,[0:0.1:5],'r');
plotSingleHist(alltog.uncueSucc_chewendings,[0:0.1:5],'k');
plotSingleHist(alltog.uncueFail_chewendings,[0:0.1:5],'b');
title('Chew durations (sec)');

figure();
plotSingleHist(alltog.cueSucc_postoutcome_reaches,[0:1:15],'g'); hold on;
plotSingleHist(alltog.cueFail_postoutcome_reaches,[0:1:15],'r');
plotSingleHist(alltog.uncueSucc_postoutcome_reaches,[0:1:15],'k');
plotSingleHist(alltog.uncueFail_postoutcome_reaches,[0:1:15],'b');
title('Number of confirmatory reaches');

figure(); plotSingleHist(postOutcome_reaches.dp_per_sess_cuedSucc_v_cuedFail,[-1-0.005:0.01:3+0.005],'k'); title('Confirmatory reaches Dprimes cuedSucc v cuedFail');
figure(); plotSingleHist(postOutcome_reaches.dp_per_sess_uncuedSucc_v_uncuedFail,[-3-0.005:0.01:3+0.005],'k'); title('Confirmatory reaches Dprimes uncuedSucc v uncuedFail');
figure(); plotSingleHist(postOutcome_reaches.dp_per_sess_cuedSucc_v_uncuedSucc,[-3-0.005:0.01:3+0.005],'k'); title('Confirmatory reaches Dprimes cuedSucc v uncuedSucc');
figure(); plotSingleHist(postOutcome_reaches.dp_per_sess_cuedFail_v_uncuedFail,[-3-0.005:0.01:3+0.005],'k'); title('Confirmatory reaches Dprimes cuedFail v uncuedFail');

figure(); plotSingleHist(chewDurations.dp_per_sess_cuedSucc_v_cuedFail,[-3-0.05:0.1:3+0.05],'k'); title('Chew durations Dprimes cuedSucc v cuedFail');
figure(); plotSingleHist(chewDurations.dp_per_sess_uncuedSucc_v_uncuedFail,[-3-0.05:0.1:3+0.05],'k'); title('Chew durations Dprimes uncuedSucc v uncuedFail');
figure(); plotSingleHist(chewDurations.dp_per_sess_cuedSucc_v_uncuedSucc,[-3-0.05:0.1:3+0.05],'k'); title('Chew durations Dprimes cuedSucc v uncuedSucc');
figure(); plotSingleHist(chewDurations.dp_per_sess_cuedFail_v_uncuedFail,[-3-0.05:0.1:3+0.05],'k'); title('Chew durations Dprimes cuedFail v uncuedFail');

end

function plotSingleHist(data,bins,col)

[n,x]=histcounts(data,bins);
[n,x]=cityscape_hist(n,x); 
plot(x,n,'Color',col);

end

function [chewendings,postoutcome_reaches,fromwhichsess_reaches,fromwhichsess_chews]=getChewsAndReachFromAllSess(data_loc_array,maxTrialsPerSess,whichReach,cueOffset,behTriggerTimeWindow,behReadoutTimeWindow)

chewendings=nan(size(data_loc_array,1)*maxTrialsPerSess,1);
postoutcome_reaches=nan(size(data_loc_array,1)*maxTrialsPerSess,1);
fromwhichsess_reaches=nan(size(data_loc_array,1)*maxTrialsPerSess,1);
fromwhichsess_chews=nan(size(data_loc_array,1)*maxTrialsPerSess,1);
chewendings_counter=1;
reach_counter=1;
reach_sess_counter=1;
chew_sess_counter=1;
for i=1:size(data_loc_array,1)
    if strcmp(data_loc_array{i,3},'no_spikes')
        continue
    end
    if strcmp(data_loc_array{i,6},'no_tbt')
        continue
    end
    disp(['loading ' data_loc_array{i,6} sep 'beh2_tbt.mat']);
    % load behavior tbt
    load([data_loc_array{i,6} sep 'beh2_tbt.mat']);
    beh2_tbt=getChewingEnds(beh2_tbt);
    beh2_tbt.cue_times=beh2_tbt.times_wrt_trial_start;
    beh2_tbt.red_time=beh2_tbt.times_wrt_trial_start;
    beh2_tbt.green_time=beh2_tbt.times_wrt_trial_start;
    % for this session
    % get number of confirmatory reaches
    reaches=getConfirmReaches(beh2_tbt,whichReach,behTriggerTimeWindow,behReadoutTimeWindow);
    postoutcome_reaches(reach_counter:reach_counter+length(reaches)-1)=reaches;
    reach_counter=reach_counter+length(reaches);
    fromwhichsess_reaches(reach_sess_counter:reach_sess_counter+length(reaches)-1)=ones(size(reaches))*i;
    reach_sess_counter=reach_sess_counter+length(reaches);
    % get duration of chewing
    chewendtimes=getChewEnd(beh2_tbt,whichReach,behTriggerTimeWindow,behReadoutTimeWindow);
    chewendings(chewendings_counter:chewendings_counter+length(chewendtimes)-1)=chewendtimes;
    chewendings_counter=chewendings_counter+length(chewendtimes);
    fromwhichsess_chews(chew_sess_counter:chew_sess_counter+length(chewendtimes)-1)=ones(size(chewendtimes))*i;
    chew_sess_counter=chew_sess_counter+length(chewendtimes);
end
reach_counter=reach_counter-1;
reach_sess_counter=reach_sess_counter-1;
chewendings_counter=chewendings_counter-1;
chew_sess_counter=chew_sess_counter-1;
chewendings=chewendings(1:chewendings_counter);
postoutcome_reaches=postoutcome_reaches(1:reach_counter);
fromwhichsess_reaches=fromwhichsess_reaches(1:reach_sess_counter);
fromwhichsess_chews=fromwhichsess_chews(1:chew_sess_counter);

end

function dp=rms_dprime(temp1,temp2)

% nonoptimal but simple
dp=(nanmean(temp1)-nanmean(temp2))./sqrt(nanstd(temp1,[],1).^2+nanstd(temp2,[],1).^2);

end

function chewendtimes=getChewEnd(beh2_tbt,reachType,reachWindow,behWindow)

[fout,dataout,n_events_in_av,alignmentCompanion,f_heatmap,plotBehFieldOut,phys_timepointsCompanion]=plotPhotometryResult(beh2_tbt,beh2_tbt,[],reachType,'chewingEnds','cueZone_onVoff','first',reachWindow,[]);
close all;
if isempty(dataout)
    chewendtimes=[];
    return
end
% get distribution of times when chewing ends
temp=nanmean(alignmentCompanion.y,1);
[~,ma]=nanmax(temp);
eventTime=alignmentCompanion.x(ma);
startAtInd=find(dataout.x>=eventTime+behWindow(1),1,'first');
temp=dataout.y(any(~isnan(alignmentCompanion.y),2),startAtInd:find(dataout.x<=eventTime+behWindow(2),1,'last'));
chewendtimes=nan(size(temp,1),1);
for i=1:size(temp,1)
    f=find(temp(i,:)>0.5,1,'first');
    if isempty(f)
        continue
    end
    chewendtimes(i)=dataout.x(startAtInd-1+f)-eventTime;
end

end

function reaches=getConfirmReaches(beh2_tbt,reachType,reachWindow,behWindow)

[fout,dataout,n_events_in_av,alignmentCompanion,f_heatmap,plotBehFieldOut,phys_timepointsCompanion]=plotPhotometryResult(beh2_tbt,beh2_tbt,[],reachType,'all_reachBatch','cueZone_onVoff','first',reachWindow,[]);
close all;
if isempty(dataout)
    reaches=[];
    return
end
% get # of reaches in behWindow, which is in seconds relative to peak of
% alignmentCompanion
temp=nanmean(alignmentCompanion.y,1);
[~,ma]=nanmax(temp);
eventTime=alignmentCompanion.x(ma);
reaches=nansum(dataout.y(any(~isnan(alignmentCompanion.y),2),find(dataout.x>=eventTime+behWindow(1),1,'first'):find(dataout.x<=eventTime+behWindow(2),1,'last')));

end