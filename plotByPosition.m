function out=plotByPosition(dataset, alltbt, trialTypes, metadata, usePositions, whichReach, cueWithinXsec)

% p_reach followed or preceded by cue for each position in a sequence

seqAt=find(dataset.realDistributions.event_isSeq{1}); % index into alltbt dim 1

temp=nanmean(alltbt.times,1); timestep=mode(diff(temp(~isnan(temp)))); [~,cueAtInd]=nanmax(nanmean(alltbt.cueZone_onVoff,1));
disp(['Using timestep: ' num2str(timestep)]);

%usePositions=-5:5;

% iterate through each position in sequence
temp=alltbt.(whichReach);
p_reach_preceded_by_cue=zeros(1,length(usePositions));
p_reach_followed_by_cue=zeros(1,length(usePositions));
ns_reach_preceded_by_cue=zeros(1,length(usePositions));
ns_reach_followed_by_cue=zeros(1,length(usePositions));
sample_std_for_binomial_precede=zeros(1,length(usePositions));
sample_std_for_binomial_follow=zeros(1,length(usePositions));
for i=1:length(usePositions)
    % get reaches
    thesepos=seqAt+usePositions(i);
    thesepos=thesepos(thesepos>0 & thesepos<=size(temp,1));
    reaches=temp(thesepos,:);
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

figure();
plot(usePositions,p_reach_preceded_by_cue,'Color','k'); hold on;
line([0 0],[min(p_reach_preceded_by_cue,[],'omitnan') max(p_reach_preceded_by_cue,[],'omitnan')],'Color','b');
for i=1:length(usePositions)
    line([usePositions(i) usePositions(i)],[p_reach_preceded_by_cue(i)-sample_std_for_binomial_precede(i) p_reach_preceded_by_cue(i)+sample_std_for_binomial_precede(i)],'Color','k');
end
title('p reach preceded by cue');
figure();
plot(usePositions,p_reach_followed_by_cue,'Color','r'); hold on;
line([0 0],[min(p_reach_followed_by_cue,[],'omitnan') max(p_reach_followed_by_cue,[],'omitnan')],'Color','b');
for i=1:length(usePositions)
    line([usePositions(i) usePositions(i)],[p_reach_followed_by_cue(i)-sample_std_for_binomial_follow(i) p_reach_followed_by_cue(i)+sample_std_for_binomial_follow(i)],'Color','r');
end
title('p reach followed by cue');

p_reach_preceded_by_cue(p_reach_preceded_by_cue==0)=1./(2.*ns_reach_preceded_by_cue(p_reach_preceded_by_cue==0));
p_reach_preceded_by_cue(p_reach_preceded_by_cue==1)=1-(1./(2.*ns_reach_preceded_by_cue(p_reach_preceded_by_cue==1)));
p_reach_followed_by_cue(p_reach_followed_by_cue==0)=1./(2.*ns_reach_followed_by_cue(p_reach_followed_by_cue==0));
p_reach_followed_by_cue(p_reach_followed_by_cue==1)=1-(1./(2.*ns_reach_followed_by_cue(p_reach_followed_by_cue==1)));
dprime_given_reach=norminv(p_reach_preceded_by_cue)-norminv(p_reach_followed_by_cue);
dprime_given_reachPLUSsd=norminv(p_reach_preceded_by_cue+sample_std_for_binomial_precede)-norminv(p_reach_followed_by_cue-sample_std_for_binomial_follow);
dprime_given_reachMINUSsd=norminv(p_reach_preceded_by_cue-sample_std_for_binomial_precede)-norminv(p_reach_followed_by_cue+sample_std_for_binomial_follow);
figure();
plot(usePositions,dprime_given_reach,'Color','g'); hold on;
plot(usePositions,dprime_given_reachPLUSsd,'Color','g');
plot(usePositions,dprime_given_reachMINUSsd,'Color','g');
line([0 0],[min(dprime_given_reachMINUSsd,[],'omitnan') max(dprime_given_reachPLUSsd,[],'omitnan')],'Color','b');

out.usePositions=usePositions;
out.dprime_given_reach=dprime_given_reach;
out.dprime_given_reachPLUSsd=dprime_given_reachPLUSsd;
out.dprime_given_reachMINUSsd=dprime_given_reachMINUSsd;
out.p_reach_preceded_by_cue=p_reach_preceded_by_cue;
out.sample_std_for_binomial_precede=sample_std_for_binomial_precede;
out.ns_reach_preceded_by_cue=ns_reach_preceded_by_cue;
out.p_reach_followed_by_cue=p_reach_followed_by_cue;
out.sample_std_for_binomial_follow=sample_std_for_binomial_follow;
out.ns_reach_followed_by_cue=ns_reach_followed_by_cue;

end