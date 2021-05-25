function alltbt=excludeExpressReaches(alltbt,trialTypes,timeFromCue,ledCond)

timesFromTrialStart=nanmean(alltbt.times_wrt_trial_start,1);
[~,ma]=nanmax(nanmean(alltbt.cueZone_onVoff,1));
cueTime=timesFromTrialStart(ma);
untilTime=cueTime+timeFromCue;
[~,untilInd]=nanmin(abs(timesFromTrialStart-untilTime));

alltbt.all_reachBatch(trialTypes.led==ledCond & any(alltbt.all_reachBatch(:,ma:untilInd)>0.05,2),ma:untilInd)=0;

end





