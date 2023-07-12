function plotByPosition(dataset, alltbt, trialTypes, metadata, usePositions, whichReach, cueWithinXsec)

seqAt=find(dataset.realDistributions.event_isSeq{1}); % index into alltbt dim 1

temp=nanmean(alltbt.times,1); timestep=mode(diff(temp(~isnan(temp)))); [~,cueAtInd]=nanmax(nanmean(alltbt.cueZone_onVoff,1));
disp(['Using timestep: ' num2str(timestep)]);

%usePositions=-5:5;

% iterate through each position in sequence
temp=alltbt.whichReach;
p_reach_preceded_by_cue=zeros(1,length(usePositions));
p_reach_followed_by_cue=zeros(1,length(usePositions));
ns_reach_preceded_by_cue=zeros(1,length(usePositions));
ns_reach_followed_by_cue=zeros(1,length(usePositions));
sample_std_for_binomial_precede=zeros(1,length(usePositions));
sample_std_for_binomial_follow=zeros(1,length(usePositions));
for i=1:length(usePositions)
    % get reaches
    reaches=temp(seqAt+usePositions(i),:);
    dataset.reaches=reaches;

    [was_reach_preceded_by_cue,trial_number,reach_time]=givenReach_probThatWasPrecededByCue(dataset, 'reaches', cueWithinXsec, timestep, cueAtInd);
    [was_reach_followed_by_cue,trial_number_follow,reach_time_follow]=givenReach_probThatWasPrecededByCue(dataset, 'reaches', -cueWithinXsec, timestep, cueAtInd);

    p_reach_preceded_by_cue(i)=nanmean(was_reach_preceded_by_cue);
    sample_std_for_binomial_precede(i)=sqrt((p_reach_preceded_by_cue(i)*(1-p_reach_preceded_by_cue(i)))/length(was_reach_preceded_by_cue));
    ns_reach_preceded_by_cue(i)=length(was_reach_preceded_by_cue);

    p_reach_followed_by_cue(i)=nanmean(was_reach_followed_by_cue);
    sample_std_for_binomial_follow(i)=sqrt((p_reach_followed_by_cue(i)*(1-p_reach_followed_by_cue(i)))/length(was_reach_followed_by_cue));
    ns_reach_followed_by_cue(i)=length(was_reach_followed_by_cue);
end


end