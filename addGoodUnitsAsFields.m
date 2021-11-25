function physiology_tbt=addGoodUnitsAsFields(physiology_tbt,spikes,goodUnitLabel,binsize,bsmooth)

% binsize in ms

maxTrialDuration=nanmax(physiology_tbt.cuetimes_wrt_trial_start(1:end));

temp=spikes.labels(:,1); 
useAssigns=unique(temp(ismember(spikes.labels(:,2),goodUnitLabel)));

figure();
for i=1:length(useAssigns)
    strunit=['unit' num2str(useAssigns(i))];
    [n,c,edges,x,y,ns]=psth_wStd_trialByTrial(filtspikes(spikes,0,'assigns',useAssigns(i)),binsize,bsmooth,maxTrialDuration,[],[]); % y is smoothed if 3rd arg is true
    physiology_tbt.(strunit)=ns;
    plot(x,y);
    hold all;
end
physiology_tbt.unitTimes=repmat(c,size(ns,1),1);
plot(nanmean(physiology_tbt.cuetimes_wrt_trial_start,1),nanmean(physiology_tbt.cue,1)*10,'Color','b');

end