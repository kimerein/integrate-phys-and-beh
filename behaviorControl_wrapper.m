function behaviorControl_wrapper(data_loc_array)

cueOffset=-0.16; % match what used for physiology
behReadoutTimeWindow=[0 5]; % in sec from alignCompanion
maxTrialsPerSess=250; % cannot be more than this many trials of any type

dp_per_sess_cuedSucc_v_cuedFail=nan(size(data_loc_array,1),1);
dp_per_sess_uncuedSucc_v_uncuedFail=nan(size(data_loc_array,1),1);
dp_per_sess_cuedSucc_v_uncuedSucc=nan(size(data_loc_array,1),1);
dp_per_sess_cuedFail_v_uncuedFail=nan(size(data_loc_array,1),1);

% Cued success
event='success_fromPerchOrWheel';
timeWindow=[0+cueOffset 3]; % in seconds from cue onset
[chewendings,postoutcome_reaches,fromwhichsess_reaches,fromwhichsess_chews]=getChewsAndReachFromAllSess(data_loc_array,maxTrialsPerSess,event,cueOffset,timeWindow,behReadoutTimeWindow);

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
    % load behavior tbt
    load([data_loc_array{i,6} sep 'beh2_tbt.mat']);
    beh2_tbt=getChewingEnds(beh2_tbt);
    % for this session
    % get number of confirmatory reaches
    reaches=getConfirmReaches(beh2_tbt,whichReach,behTriggerTimeWindow,behReadoutTimeWindow);
    postoutcome_reaches(reach_counter:reach_counter+length(reaches)-1)=reaches;
    reach_counter=reach_counter+length(reaches);
    fromwhichsess_reaches(reach_sess_counter:reach_sess_counter+length(reaches)-1)=ones(size(reaches))*i;
    reach_sess_counter=reach_sess_counter+length(reaches);
    % get duration of chewing
    chewendtimes=getChewEnd(beh2_tbt,whichReach,[0+cueOffset 3],[0 5]);
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
% get distribution of times when chewing ends
temp=nanmean(alignmentCompanion.y,1);
[~,ma]=nanmax(temp);
eventTime=alignmentCompanion.x(ma);
startAtInd=find(dataout.x>=eventTime+behWindow(1),1,'first');
temp=dataout.y(any(~isnan(alignmentCompanion.y),2),startAtInd:find(dataout.x<=eventTime+behWindow(2),1,'last'));
chewendtimes=nan(size(temp,1),1);
for i=1:size(temp,1)
    f=find(temp(i,:)>0.5,1,'first');
    chewendtimes(i)=dataout.x(startAtInd-1+f)-eventTime;
end

end

function reaches=getConfirmReaches(beh2_tbt,reachType,reachWindow,behWindow)

[fout,dataout,n_events_in_av,alignmentCompanion,f_heatmap,plotBehFieldOut,phys_timepointsCompanion]=plotPhotometryResult(beh2_tbt,beh2_tbt,[],reachType,'all_reachBatch','cueZone_onVoff','first',reachWindow,[]);
close all;
% get # of reaches in behWindow, which is in seconds relative to peak of
% alignmentCompanion
temp=nanmean(alignmentCompanion.y,1);
[~,ma]=nanmax(temp);
eventTime=alignmentCompanion.x(ma);
reaches=nansum(dataout.y(any(~isnan(alignmentCompanion.y),2),find(dataout.x>=eventTime+behWindow(1),1,'first'):find(dataout.x<=eventTime+behWindow(2),1,'last')));

end