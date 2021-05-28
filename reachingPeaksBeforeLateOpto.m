function trialTypes=reachingPeaksBeforeLateOpto(alltbt,metadata,trialTypes)

% assumes unique sessids for different days and mice
optoLength=0.5; % in seconds

trialTypes.reachingPeaksBeforeLateOpto=nan(size(trialTypes.led));
u=unique(metadata.sessid);
[~,cueInd]=nanmax(nanmean(alltbt.cueZone_onVoff,1));
timestep=mode(diff(nanmean(alltbt.times_wrt_trial_start,1)));
cueTime=timestep*cueInd;
for i=1:length(u)
    curr_sessid=u(i);
    figure(); 
    plot(nanmean(alltbt.times_wrt_trial_start,1),nanmean(alltbt.all_reachBatch(metadata.sessid==curr_sessid,:),1),'Color','k');
    hold on;
    plot(nanmean(alltbt.times_wrt_trial_start,1),nanmean(alltbt.cueZone_onVoff,1)/10,'Color','b');
    line([cueTime+optoLength cueTime+optoLength],[0 nanmax(nanmean(alltbt.all_reachBatch(metadata.sessid==curr_sessid,:),1))],'Color','r');
    line([cueTime cueTime],[0 nanmax(nanmean(alltbt.all_reachBatch(metadata.sessid==curr_sessid,:),1))],'Color','r');
    answer=questdlg('Is reaching peak before late opto?');
    switch answer
        case 'Yes'
            trialTypes.reachingPeaksBeforeLateOpto(metadata.sessid==curr_sessid)=1;
        case 'No'
            trialTypes.reachingPeaksBeforeLateOpto(metadata.sessid==curr_sessid)=0;
        otherwise
            return
    end
end

end