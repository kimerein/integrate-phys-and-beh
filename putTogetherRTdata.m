function putTogetherRTdata(alltbt,out,metadata)

bins=500;

unique_sess=unique(metadata.sessid);

f=fieldnames(alltbt);
fout=fieldnames(out);

allout.cued1back_touched1back_noPreempt.ledFalse=cell(1,length(unique_sess));
allout.cued1back_touched1back_noPreempt.ledTrue=cell(1,length(unique_sess));
allout.cued1back_touched1back.ledFalse=cell(1,length(unique_sess));
allout.cued1back_touched1back.ledTrue=cell(1,length(unique_sess));
allout.cued1back_failed1back.ledFalse=cell(1,length(unique_sess));
allout.cued1back_failed1back.ledTrue=cell(1,length(unique_sess));
allout.cued1back_touched1back_didnteat1back.ledFalse=cell(1,length(unique_sess));
allout.cued1back_touched1back_didnteat1back.ledTrue=cell(1,length(unique_sess));
allout.cued1back_touched1back_nochewTrialstart.ledFalse=cell(1,length(unique_sess));
allout.cued1back_touched1back_nochewTrialstart.ledTrue=cell(1,length(unique_sess));
allout.alltrials_nochewTrialstart.ledFalse=cell(1,length(unique_sess));
allout.alltrials_nochewTrialstart.ledTrue=cell(1,length(unique_sess));
allout.allreactiontimes.ledFalse=cell(1,length(unique_sess));
allout.allreactiontimes.ledTrue=cell(1,length(unique_sess));

for i=1:length(unique_sess)
    currsessid=unique_sess(i);
    for j=1:length(f)
        temp=alltbt.(f{j});
        currtbt.(f{j})=temp(metadata.sessid==currsessid,:);
    end
    for j=1:length(fout)
        temp=out.(fout{j});
        currout.(fout{j})=temp(metadata.sessid==currsessid);
    end
    output=getRTanalysis(currtbt,currout);
    if length(output.cued1back_touched1back.ledFalse)==length(output.cued1back_touched1back.ledTrue)
        allout.cued1back_touched1back_noPreempt.ledFalse{i}=output.cued1back_touched1back_noPreempt.ledFalse;
        allout.cued1back_touched1back_noPreempt.ledTrue{i}=output.cued1back_touched1back_noPreempt.ledTrue;
        allout.cued1back_touched1back.ledFalse{i}=output.cued1back_touched1back.ledFalse;
        allout.cued1back_touched1back.ledTrue{i}=output.cued1back_touched1back.ledTrue;
        allout.cued1back_failed1back.ledFalse{i}=output.cued1back_failed1back.ledFalse;
        allout.cued1back_failed1back.ledTrue{i}=output.cued1back_failed1back.ledTrue;
        allout.cued1back_touched1back_didnteat1back.ledFalse{i}=output.cued1back_touched1back_didnteat1back.ledFalse;
        allout.cued1back_touched1back_didnteat1back.ledTrue{i}=output.cued1back_touched1back_didnteat1back.ledTrue;
        allout.cued1back_touched1back_nochewTrialstart.ledFalse{i}=output.cued1back_touched1back_nochewTrialstart.ledFalse;
        allout.cued1back_touched1back_nochewTrialstart.ledTrue{i}=output.cued1back_touched1back_nochewTrialstart.ledTrue;
        allout.alltrials_nochewTrialstart.ledFalse{i}=output.alltrials_nochewTrialstart.ledFalse;
        allout.alltrials_nochewTrialstart.ledTrue{i}=output.alltrials_nochewTrialstart.ledTrue;
%         allout.allreactiontimes.ledFalse{i}=output.reactionTimes.ledFalse;
%         allout.allreactiontimes.ledTrue{i}=output.reactionTimes.ledTrue;
        fracThrough=metadata.fractionThroughSess(metadata.sessid==currsessid);
        allout.allreactiontimes.ledFalse{i}=output.reactionTimes.ledFalse(fracThrough>=0.8);
        allout.allreactiontimes.ledTrue{i}=output.reactionTimes.ledTrue(fracThrough<=0.2);
%         allout.allreactiontimes.ledFalse{i}=output.reactionTimes.ledFalse;
%         allout.allreactiontimes.ledTrue{i}=output.reactionTimes.ledTrue;
    else
        allout.cued1back_touched1back_noPreempt.ledFalse{i}=nan;
        allout.cued1back_touched1back_noPreempt.ledTrue{i}=nan;
        allout.cued1back_touched1back.ledFalse{i}=nan;
        allout.cued1back_touched1back.ledTrue{i}=nan;
        allout.cued1back_failed1back.ledFalse{i}=nan;
        allout.cued1back_failed1back.ledTrue{i}=nan;
        allout.cued1back_touched1back_didnteat1back.ledFalse{i}=nan;
        allout.cued1back_touched1back_didnteat1back.ledTrue{i}=nan;
        allout.cued1back_touched1back_nochewTrialstart.ledFalse{i}=nan;
        allout.cued1back_touched1back_nochewTrialstart.ledTrue{i}=nan;
        allout.alltrials_nochewTrialstart.ledFalse{i}=nan;
        allout.alltrials_nochewTrialstart.ledTrue{i}=nan;
        allout.allreactiontimes.ledFalse{i}=nan;
        allout.allreactiontimes.ledTrue{i}=nan;
    end
end

% Choose which to use
% success_ledFalse=allout.cued1back_touched1back.ledFalse;
% success_ledTrue=allout.cued1back_touched1back.ledTrue;
% fail_ledFalse=allout.cued1back_failed1back.ledFalse;
% fail_ledTrue=allout.cued1back_failed1back.ledTrue;
success_ledFalse=allout.cued1back_touched1back_nochewTrialstart.ledFalse;
success_ledTrue=allout.cued1back_touched1back_nochewTrialstart.ledTrue;
fail_ledFalse=allout.cued1back_failed1back.ledFalse;
fail_ledTrue=allout.cued1back_failed1back.ledTrue;
all_ledFalse=allout.alltrials_nochewTrialstart.ledFalse;
all_ledTrue=allout.alltrials_nochewTrialstart.ledTrue;
allRT_ledFalse=allout.allreactiontimes.ledFalse;
allRT_ledTrue=allout.allreactiontimes.ledTrue;

for i=1:length(unique_sess)
    success_med_ledFalse(i)=nanmedian(success_ledFalse{i});
    success_med_ledTrue(i)=nanmedian(success_ledTrue{i});
    fail_med_ledFalse(i)=nanmedian(fail_ledFalse{i});
    fail_med_ledTrue(i)=nanmedian(fail_ledTrue{i}); 
    all_med_ledFalse(i)=nanmean(all_ledFalse{i});
    all_med_ledTrue(i)=nanmean(all_ledTrue{i});
end

figure();
for i=1:length(unique_sess)
    plot([1 2],[success_med_ledFalse(i) success_med_ledTrue(i)]);
    hold on;
end
disp('p-value comparing medians');
if all(isnan(success_med_ledFalse) | isnan(success_med_ledTrue))
    disp('not enough data');
else
    disp(signrank(success_med_ledFalse, success_med_ledTrue));
end

% Put together all reaction times, but now subtract off median of fail
all_success_ledFalse=[];
all_success_ledTrue=[];
all_fails_ledFalse=[];
all_fails_ledTrue=[];
all_RT_ledFalse=[];
all_RT_ledTrue=[];
for i=1:length(unique_sess)
%     all_success_ledFalse=[all_success_ledFalse success_ledFalse{i}-nanmean([all_med_ledFalse(i)])];
%     all_success_ledTrue=[all_success_ledTrue success_ledTrue{i}-nanmean([all_med_ledFalse(i)])];
%     all_fails_ledFalse=[all_fails_ledFalse fail_ledFalse{i}-nanmean([all_med_ledFalse(i)])];
%     all_fails_ledTrue=[all_fails_ledTrue fail_ledTrue{i}-nanmean([all_med_ledFalse(i)])];
    all_success_ledFalse=[all_success_ledFalse success_ledFalse{i}];
    all_success_ledTrue=[all_success_ledTrue success_ledTrue{i}];
    all_fails_ledFalse=[all_fails_ledFalse fail_ledFalse{i}];
    all_fails_ledTrue=[all_fails_ledTrue fail_ledTrue{i}];
    all_RT_ledFalse=[all_RT_ledFalse allRT_ledFalse{i}];
    all_RT_ledTrue=[all_RT_ledTrue allRT_ledTrue{i}];
end


% Plot results
% [n,x]=histcounts(all_success_ledFalse,bins);
% [n,x]=histcounts(all_RT_ledFalse,0-0.015:0.03:15-0.03);
[n,x]=histcounts(all_RT_ledFalse,0-0.03:0.06:15-0.03);
% [n,x]=histcounts(all_RT_ledFalse,0-0.015:0.045:15-0.03);
x_backup=x;
[n,x]=cityscape_hist(n,x);
figure();
% plot(x,n,'Color',[0.8 0.8 0.8]);
plot(x,n./nansum(n),'Color',[0.8 0.8 0.8]);
hold on;
[n,x]=histcounts(all_RT_ledTrue,x_backup);
[n,x]=cityscape_hist(n,x);
% plot(x,n,'Color',[0.8 0.5 0.5]);
plot(x,n./nansum(n),'Color',[0.8 0.5 0.5]);
leg={'LED FALSE','LED TRUE'};
legend(leg);
xlabel('RT (sec)');
ylabel('Count');
disp('p-value of all together');
if all(isnan(all_RT_ledFalse)) | all(isnan(all_RT_ledTrue))
    disp('not enough data');
else
    disp(ranksum(all_RT_ledFalse,all_RT_ledTrue));
end
figure();
hax=gca;
parmhat=fitAndPlot(hax,all_RT_ledFalse,'lognfit',0-0.03:0.06:15-0.03,'k');
fitAndPlot(hax,all_RT_ledTrue,'lognfit',0-0.03:0.06:15-0.03,'r');

% Plot results
% [n,x]=histcounts(all_success_ledFalse,bins);
% [n,x]=histcounts(all_success_ledFalse,-100+0.03:0.06:100-0.06);
% [n,x]=histcounts(all_success_ledFalse,0:0.15:100-0.1);
% [n,x]=histcounts(all_success_ledFalse,-0.08:0.18:100-0.1);
% [n,x]=histcounts(all_success_ledFalse,-0.125:0.15:100-0.1);
% [n,x]=histcounts(all_success_ledFalse,-0.3:0.6:100-0.1);
[n,x]=histcounts(all_success_ledFalse,-0.1:0.2:100-0.1);
x_backup=x; 
[n,x]=cityscape_hist(n,x);
figure();
plot(x,n./nansum(n),'Color','k');
% plot(x,n,'Color','k');
hold on;
[n,x]=histcounts(all_success_ledTrue,x_backup);
[n,x]=cityscape_hist(n,x);
plot(x,n./nansum(n),'Color','r');
% plot(x,n,'Color','r');
leg={'LED FALSE','LED TRUE'};
legend(leg);
xlabel('Change in RT (sec)');
ylabel('Count');
disp('p-value of all together');
if all(isnan(success_med_ledFalse)) | all(isnan(success_med_ledTrue))
    disp('not enough data');
else
    % Note that if mouse gets more than 1.5 seconds faster on next trial,
    % animal is reaching BEFORE the cue, given that he did a cued reach on
    % previous trial
%     disp(ranksum(all_success_ledFalse(all_success_ledFalse<=1.5 & all_success_ledFalse>=-1.5),all_success_ledTrue(all_success_ledTrue<=1.5 & all_success_ledTrue>=-1.5)));
    disp(ranksum(all_success_ledFalse,all_success_ledTrue));
end
figure();
hax=gca;
% parmhat_red=fitAndPlot(hax,all_success_ledTrue,'lognfit',0:0.18:100-0.1,'r'); % parmhat_red should account for regression to the mean
% parmhat_red=fitAndPlot(hax,all_success_ledTrue,'lognfit',0:0.3:100-0.1,'r'); % parmhat_red should account for regression to the mean
% Note bins should be centered over 1!
parmhat_red=fitAndPlot(hax,all_success_ledTrue,'lognfit',-0.1:0.2:100-0.1,'r'); % parmhat_red should account for regression to the mean
disp(['parmhat_red mean is ' num2str(parmhat_red(1))]);
disp(['parmhat_red var is ' num2str(parmhat_red(2))]);
% fitAndPlot(hax,all_success_ledFalse,'sum_2_lognfit',0:0.18:100-0.1,'k',[parmhat_red(1) parmhat_red(2)]);
% fitAndPlot(hax,all_success_ledFalse,'sum_2_lognfit',0:0.3:100-0.1,'k',[parmhat_red(1) parmhat_red(2)]);
% Note bins should be centered over 1!
fitAndPlot(hax,all_success_ledFalse,'sum_2_lognfit',-0.1:0.2:100-0.1,'k',[parmhat_red(1) parmhat_red(2)]);

% [n,x]=histcounts(all_success_ledFalse,bins);
[n,x]=histcounts(all_fails_ledFalse,-100+0.03:0.06:100-0.06);
x_backup=x;
[n,x]=cityscape_hist(n,x);
figure();
plot(x,n./nansum(n),'Color','k');
hold on;
[n,x]=histcounts(all_fails_ledTrue,x_backup);
[n,x]=cityscape_hist(n,x);
plot(x,n./nansum(n),'Color','c');
leg={'LED FALSE','LED TRUE'};
legend(leg);
xlabel('Change in RT (sec)');
ylabel('Count');
disp('p-value of all together');
if all(isnan(success_med_ledFalse)) | all(isnan(success_med_ledTrue))
    disp('not enough data');
else
    disp(ranksum(all_fails_ledFalse,all_fails_ledTrue));
end


end

function output=getRTanalysis(alltbt,out)

nbins=[];

% Get reaction times for all trials where mouse reached after cue onset
[reactionTimes,alltbt]=plotOnlyFirstReach(alltbt,1,'reachStarts_noPawOnWheel','cueZone_onVoff',out,'led',0);
% output.reactionTimes.ledFalse=reactionTimes(out.led_1back==0 & out.cued_reach_1back==1 & out.touched_pellet_1back==1 & out.chewing_at_trial_start==0);
% output.reactionTimes.ledTrue=reactionTimes(out.led_1back==1 & out.cued_reach_1back==1 & out.touched_pellet_1back==1 & out.chewing_at_trial_start==0);
temp=reactionTimes;
% temp(out.led_1back==1 & out.cued_reach_1back==1 & out.touched_pellet_1back==1 & out.chewing_at_trial_start==0)=nan;
temp(~(out.led_1back==0 & out.led==0 & out.chewing_at_trial_start==0))=nan;
output.reactionTimes.ledFalse=temp;
temp=reactionTimes;
% temp(out.led_1back==0 & out.cued_reach_1back==1 & out.touched_pellet_1back==1 & out.chewing_at_trial_start==0)=nan;
temp(~(out.led_1back==0 & out.led==0 & out.chewing_at_trial_start==0))=nan;
output.reactionTimes.ledTrue=temp;
% temp(out.led_1back==1)=nan;
% output.reactionTimes.ledFalse=temp;
% temp=reactionTimes;
% temp(out.led_1back==0)=nan;
% output.reactionTimes.ledTrue=temp;

% Get trials where mouse paw was on pellet presenter wheel during wheel turn
% Note that enforcing this requires that the mouse not have extra
% information about pellet presentation / cue timing
dontUse=out.paw_during_wheel==1;
dontUse_1back=out.paw_during_wheel_1back==1;
dontUse_sequence=dontUse | dontUse_1back;

% Change in reaction time, no preemptive reach
% Test effects of LED on PREVIOUS trial 
% given that mouse successfully touched pellet and performed cued reach on
% previous trial
condition=out.cued_reach_1back==1 & out.touched_pellet_1back==1;
curr_rt=reactionTimes;
curr_rt(~(condition==1) | dontUse_sequence)=nan;
% plotCurrRT(curr_rt,'Cued reach & touched pellet 1 back & no preemptive reach for 2 trials in a row',out.led_1back==1);
rtchange=histRTchange_local(curr_rt,nbins,'Cued reach & touched pellet 1 back & no preemptive reach for 2 trials in a row',out.led_1back==1);
output.cued1back_touched1back_noPreempt.ledFalse=rtchange.rt_change_testcond0;
output.cued1back_touched1back_noPreempt.ledTrue=rtchange.rt_change_testcond1;
% [p,h]=ranksum(rtchange.rt_change_testcond0,rtchange.rt_change_testcond1);
% disp('p-value of cued reach & touched pellet 1 back & no preemptive reach for 2 trials in a row');
% disp(p);

% Change in reaction time
% Test effects of LED on PREVIOUS trial 
% given that mouse successfully touched pellet and performed cued reach on
% previous trial
condition=out.cued_reach_1back==1 & out.touched_pellet_1back==1;
curr_rt=reactionTimes;
curr_rt(~(condition==1))=nan;
% plotCurrRT(curr_rt,'Cued reach & touched pellet 1 back',out.led_1back==1);
% rtchange=histRTchange_local(curr_rt,nbins,'Cued reach & touched pellet 1 back',out.led_1back==1);
rtchange=histRTchange_local(curr_rt,nbins,'Cued reach & touched pellet 1 back',out.led_1back==1);
output.cued1back_touched1back.ledFalse=rtchange.rt_change_testcond0;
output.cued1back_touched1back.ledTrue=rtchange.rt_change_testcond1;
% [p,h]=ranksum(rtchange.rt_change_testcond0,rtchange.rt_change_testcond1);
% disp('p-value of cued reach & touched pellet 1 back');
% disp(p);

% Change in reaction time
% Test effects of LED on PREVIOUS trial 
% given that mouse failed to touch pellet despite performing a cued reach on
% previous trial
condition=out.cued_reach_1back==1 & out.touched_pellet_1back==0;
curr_rt=reactionTimes;
% curr_rt(~(condition==1) | dontUse_sequence)=nan;
curr_rt(~(condition==1))=nan;
% plotCurrRT(curr_rt,'Cued reach but failed to touch pellet 1 back & no preemptive reach for 2 trials in a row',out.led_1back==1);
% rtchange=histRTchange_local(curr_rt,nbins,'Cued reach but failed to touch pellet 1 back & no preemptive reach for 2 trials in a row',out.led_1back==1);
rtchange=histRTchange_local(curr_rt,nbins,'Cued reach but failed to touch pellet 1 back & no preemptive reach for 2 trials in a row',out.led_1back==1);
output.cued1back_failed1back.ledFalse=rtchange.rt_change_testcond0;
output.cued1back_failed1back.ledTrue=rtchange.rt_change_testcond1;

condition=out.cued_reach_1back==1 & out.touched_pellet_1back==1 & out.consumed_pellet_1back==0;
curr_rt=reactionTimes;
curr_rt(~(condition==1))=nan;
% plotCurrRT(curr_rt,'Playing 2',out.led_1back==1);
rtchange=histRTchange_local(curr_rt,nbins,'Playing 2',out.led_1back==1);
output.cued1back_touched1back_didnteat1back.ledFalse=rtchange.rt_change_testcond0;
output.cued1back_touched1back_didnteat1back.ledTrue=rtchange.rt_change_testcond1;

condition=out.cued_reach_1back==1 & out.touched_pellet_1back==1 & out.chewing_at_trial_start==0;
% condition=out.cued_reach_1back==1 & out.touched_pellet_1back==1;
curr_rt=reactionTimes;
curr_rt(~(condition==1))=nan;
% plotCurrRT(curr_rt,'Playing 3',out.led_1back==1);
% rtchange=histRTchange_local(curr_rt,nbins,'Playing 3',out.led_1back==1);
rtchange=histRTchange_local(curr_rt,nbins,'Playing 3',out.led_1back==1); % here
output.cued1back_touched1back_nochewTrialstart.ledFalse=rtchange.rt_change_testcond0;
output.cued1back_touched1back_nochewTrialstart.ledTrue=rtchange.rt_change_testcond1;

condition=out.chewing_at_trial_start==0;
% condition=ones(size(out.chewing_at_trial_start==0));
curr_rt=reactionTimes;
curr_rt(~(condition==1))=nan;
% plotCurrRT(curr_rt,'Playing 3',out.led_1back==1);
% rtchange=histRTchange_local(curr_rt,nbins,'Playing 3',out.led_1back==1);
rtchange=histRTchange_local(curr_rt,nbins,'Playing 4',out.led_1back==1); % here
output.alltrials_nochewTrialstart.ledFalse=rtchange.rt_change_testcond0;
output.alltrials_nochewTrialstart.ledTrue=rtchange.rt_change_testcond1;


end

function out=histRTchange_local(curr_rt,bins,tit,testcond)

out.rt_change_testcond0=nan;
out.rt_change_testcond1=nan;

if ~isempty(testcond)
    backup_rt=curr_rt;
    curr_rt(~(testcond==0))=nan; % test condition is false
end

if all(isnan(curr_rt(1:end-1)-curr_rt(2:end)))
    return
end
% [n,x]=histcounts(curr_rt(1:end-1)-curr_rt(2:end),bins);
% SUM
% out.rt_change_testcond0=curr_rt(1:end-1)-curr_rt(2:end);
% MULTIPLY
out.rt_change_testcond0=curr_rt(2:end)./curr_rt(1:end-1);
% x_backup=x;
% [n,x]=cityscape_hist(n,x);
% figure();
% plot(x,n./nansum(n),'Color','k');
% xlabel('Change in RT (sec)');
% ylabel('Count');
% title(tit);

if ~isempty(testcond)
%     hold on;
    curr_rt=backup_rt;
    curr_rt(~(testcond==1))=nan; % test condition is true
    if all(isnan(curr_rt(1:end-1)-curr_rt(2:end)))
        return
    end
%     [n,x]=histcounts(curr_rt(1:end-1)-curr_rt(2:end),x_backup);
    % SUM
%     out.rt_change_testcond1=curr_rt(1:end-1)-curr_rt(2:end);
    % MULTIPLY
    out.rt_change_testcond1=curr_rt(2:end)./curr_rt(1:end-1);
%     [n,x]=cityscape_hist(n,x);
%     plot(x,n./nansum(n),'Color','r');
%     leg={'testcond FALSE','testcond TRUE'};
%     legend(leg);
end

end