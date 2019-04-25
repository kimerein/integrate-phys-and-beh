function analyzeRTchanges(alltbt,out,metadata)

nbins=500; % for histograms

% Get reaction times for all trials where mouse reached after cue onset
[reactionTimes,alltbt]=plotOnlyFirstReach(alltbt,1,'reachStarts_noPawOnWheel','cueZone_onVoff',out,'led',0);

% Get trials where mouse paw was on pellet presenter wheel during wheel turn
% Note that enforcing this requires that the mouse not have extra
% information about pellet presentation / cue timing
dontUse=out.paw_during_wheel==1;
dontUse_1back=out.paw_during_wheel_1back==1;
dontUse_sequence=dontUse | dontUse_1back;

% Change in reaction time, all trials
condition=ones(length(out.led),1); % take all trials
curr_rt=reactionTimes;
curr_rt(~(condition==1))=nan;
plotCurrRT(curr_rt,'All trials',[]);
histRTchange(curr_rt,nbins,'All trials',[]);

% Change in reaction time, no preemptive reach
condition=ones(length(out.led),1); % take all trials
curr_rt=reactionTimes;
curr_rt(~(condition==1) | dontUse_sequence)=nan;
plotCurrRT(curr_rt,'All trials, no reach during wheel turn (2 trials in a row)',[]);
histRTchange(curr_rt,nbins,'All trials, no reach during wheel turn (2 trials in a row)',[]);

% Change in reaction time, no preemptive reach
% Test effects of LED on current trial
condition=ones(length(out.led),1); % take all trials
curr_rt=reactionTimes;
curr_rt(~(condition==1) | dontUse_sequence)=nan;
plotCurrRT(curr_rt,'Trials w and wout current trial led on, no reach during wheel turn (2 trials in a row)',out.led==1);
histRTchange(curr_rt,nbins,'Trials w and wout current trial led on, no reach during wheel turn (2 trials in a row)',out.led==1);

% Change in reaction time, no preemptive reach
% Test effects of LED on PREVIOUS trial
condition=ones(length(out.led),1); % take all trials
curr_rt=reactionTimes;
curr_rt(~(condition==1) | dontUse_sequence)=nan;
plotCurrRT(curr_rt,'Trials w and wout previous trial led on, no reach during wheel turn (2 trials in a row)',out.led_1back==1);
histRTchange(curr_rt,nbins,'Trials w and wout previous trial led on, no reach during wheel turn (2 trials in a row)',out.led_1back==1);

% Change in reaction time, no preemptive reach
% Test effects of LED on PREVIOUS trial 
% given that mouse successfully touched pellet and performed cued reach on
% previous trial
condition=out.cued_reach_1back==1 & out.touched_pellet_1back==1;
curr_rt=reactionTimes;
curr_rt(~(condition==1) | dontUse_sequence)=nan;
plotCurrRT(curr_rt,'Cued reach & touched pellet 1 back & no preemptive reach for 2 trials in a row',out.led_1back==1);
rtchange=histRTchange(curr_rt,nbins,'Cued reach & touched pellet 1 back & no preemptive reach for 2 trials in a row',out.led_1back==1);
try
    [p,h]=ranksum(rtchange.rt_change_testcond0,rtchange.rt_change_testcond1);
    disp('p-value of cued reach & touched pellet 1 back & no preemptive reach for 2 trials in a row');
    disp(p);
catch
end

% Change in reaction time, no preemptive reach
% Test effects of LED on PREVIOUS trial 
% given that mouse successfully touched pellet and performed cued reach on
% previous trial
condition=out.cued_reach_1back==1 & out.touched_pellet_1back==1;
curr_rt=reactionTimes;
curr_rt(~(condition==1))=nan;
plotCurrRT(curr_rt,'Cued reach & touched pellet 1 back',out.led_1back==1);
rtchange=histRTchange(curr_rt,nbins,'Cued reach & touched pellet 1 back',out.led_1back==1);
try
    [p,h]=ranksum(rtchange.rt_change_testcond0,rtchange.rt_change_testcond1);
    disp('p-value of cued reach & touched pellet 1 back');
    disp(p);
end

% Change in reaction time, no preemptive reach
% Test effects of LED on PREVIOUS trial 
% given that mouse failed to touch pellet despite performing a cued reach on
% previous trial
condition=out.cued_reach_1back==1 & out.touched_pellet_1back==0;
curr_rt=reactionTimes;
% curr_rt(~(condition==1) | dontUse_sequence)=nan;
curr_rt(~(condition==1))=nan;
plotCurrRT(curr_rt,'Cued reach but failed to touch pellet 1 back & no preemptive reach for 2 trials in a row',out.led_1back==1);
histRTchange(curr_rt,nbins,'Cued reach but failed to touch pellet 1 back & no preemptive reach for 2 trials in a row',out.led_1back==1);

% Playing
condition=out.cued_reach_1back==1 & out.touched_pellet_1back==1 & (out.paw_during_wheel==0 | out.paw_during_wheel_1back==0);
curr_rt=reactionTimes;
curr_rt(~(condition==1))=nan;
plotCurrRT(curr_rt,'Playing',out.led_1back==1);
rtchange=histRTchange(curr_rt,nbins,'Playing',out.led_1back==1);
try
    [p,h]=ranksum(rtchange.rt_change_testcond0,rtchange.rt_change_testcond1);
    disp('p-value of playing');
    disp(p);
catch
end

% Playing 2
condition=out.cued_reach_1back==1 & out.touched_pellet_1back==1 & out.consumed_pellet_1back==0;
% condition=out.dprime>0.1 & out.cued_reach_1back==1 & out.touched_pellet_1back==1;
curr_rt=reactionTimes;
curr_rt(~(condition==1))=nan;
plotCurrRT(curr_rt,'Playing 2',out.led_1back==1);
rtchange=histRTchange(curr_rt,nbins,'Playing 2',out.led_1back==1);
try
    [p,h]=ranksum(rtchange.rt_change_testcond0,rtchange.rt_change_testcond1);
    disp('p-value of playing 2');
    disp(p);
catch
end

% Playing 3
condition=out.touched_pellet_1back==1 & (out.consumed_pellet_1back==0 | out.chewing_at_trial_start==0);
% condition=out.cued_reach_1back==1 & out.touched_pellet_1back==1 & out.chewing_at_trial_start==0 & out.paw_during_wheel==0 & out.paw_during_wheel_1back==0;
% condition=out.chewing_at_trial_start==0;
curr_rt=reactionTimes;
curr_rt(~(condition==1))=nan;
plotCurrRT(curr_rt,'Playing 3',out.led_1back==1);
rtchange=histRTchange(curr_rt,nbins,'Playing 3',out.led_1back==1);
try
    [p,h]=ranksum(rtchange.rt_change_testcond0,rtchange.rt_change_testcond1);
    disp('p-value of playing 3');
    disp(p);
catch
end
% figure();
% plot(nanmean(alltbt.reachStarts_noPawOnWheel(condition==1 & out.led_1back==0,:),1),'Color','k');
% hold on;
% plot(nanmean(alltbt.reachStarts_noPawOnWheel(condition==1 & out.led_1back==1,:),1),'Color','r');
% plot(nanmean(alltbt.cueZone_onVoff,1),'Color','b');
figure();
plot(nanmean(alltbt.times,1),nanmean(alltbt.miss_reachStarts(out.chewing_at_trial_start==0 & out.led_1back==0,:),1),'Color','k');
hold on;
plot(nanmean(alltbt.times,1),nanmean(alltbt.miss_reachStarts(out.chewing_at_trial_start==0 & out.led_1back==1,:),1),'Color','r');
plot(nanmean(alltbt.times,1),nanmean(alltbt.cueZone_onVoff,1),'Color','b');

condition=out.cued_reach_1back==1 & out.touched_pellet_1back==1;
x=nanmean(alltbt.times,1);
% y1=nanmean(alltbt.reachStarts_noPawOnWheel(out.chewing_at_trial_start==0 & out.led_1back==0,:),1);
% y2=nanmean(alltbt.reachStarts_noPawOnWheel(out.chewing_at_trial_start==0 & out.led_1back==1,:),1);
x=0:0.025:max(reactionTimes);
[y1,x]=histcounts(reactionTimes(out.chewing_at_trial_start==0 & out.led_1back==0),x);
xbackup=x;
[y2,x]=histcounts(reactionTimes(out.chewing_at_trial_start==0 & out.led_1back==1),x);
% [y1,x]=cityscape_hist(y1,[xbackup-mode(diff(xbackup)) xbackup(end)+mode(diff(xbackup))]);
[y1,x]=cityscape_hist(y1,xbackup);
figure();
plot(x,y1./nansum(y1(x<=1.5)),'Color','k'); 
hold on;
% [y2,x]=cityscape_hist(y2,[xbackup-mode(diff(xbackup)) xbackup(end)+mode(diff(xbackup))]);
[y2,x]=cityscape_hist(y2,xbackup);
plot(x,y2./nansum(y2(x<=1.5)),'Color','r');

end

function plotCurrRT(curr_rt,tit,testcond)

if ~isempty(testcond)
    backup_rt=curr_rt;
    curr_rt(~(testcond==0))=nan; % test condition is false
end

if all(isnan(curr_rt(1:end-1)) | isnan(curr_rt(2:end)))
    return
end
figure();
jitter=0.03;
% jitter=0;
% disp('slope when passes through zero');
% yvals=curr_rt(2:end)+rand(size(curr_rt(2:end)));
% xvals=curr_rt(1:end-1)+rand(size(curr_rt(1:end-1)));
% yvals(yvals>1.5)=nan;
% xvals(xvals>1.5)=nan;
% disp(nanmean(yvals./xvals));
s=scatter(curr_rt(1:end-1)+rand(size(curr_rt(1:end-1))).*jitter,curr_rt(2:end)+rand(size(curr_rt(2:end))).*jitter,[],'k','filled');
m=get(s,'MarkerHandle');
% alpha=0.3;
% m.FaceColorData=uint8(255*[0;0;0;alpha]);
%m.EdgeColorData=uint8(255*[0;0;0;alpha]);
xlabel('RT previous trial in sec');
ylabel('RT current trial in sec');
title(tit);

if ~isempty(testcond)
    hold on;
    curr_rt=backup_rt;
    curr_rt(~(testcond==1))=nan; % test condition is true
    if all(isnan(curr_rt(1:end-1)) | isnan(curr_rt(2:end)))
        return
    end
%     disp('slope when passes through zero');
%     yvals=curr_rt(2:end)+rand(size(curr_rt(2:end)));
%     xvals=curr_rt(1:end-1)+rand(size(curr_rt(1:end-1)));
%     yvals(yvals>1.5)=nan;
%     xvals(xvals>1.5)=nan;
%     disp(nanmean(yvals./xvals));
    s=scatter(curr_rt(1:end-1)+rand(size(curr_rt(1:end-1))).*jitter,curr_rt(2:end)+rand(size(curr_rt(2:end))).*jitter,[],'r','filled');
    m=get(s,'MarkerHandle');
%     alpha=0.3;
%     m.FaceColorData=uint8(255*[1;0;0;alpha]);
    leg={'testcond FALSE','testcond TRUE'};
    legend(leg);
end

end

function out=histRTchange(curr_rt,bins,tit,testcond)

if ~isempty(testcond)
    backup_rt=curr_rt;
    curr_rt(~(testcond==0))=nan; % test condition is false
end

if all(isnan(curr_rt(1:end-1)-curr_rt(2:end)))
    return
end
[n,x]=histcounts(curr_rt(1:end-1)-curr_rt(2:end),bins);
out.rt_change_testcond0=curr_rt(1:end-1)-curr_rt(2:end);
x_backup=x;
[n,x]=cityscape_hist(n,x);
figure();
plot(x,n./nansum(n),'Color','k');
xlabel('Change in RT (sec)');
ylabel('Count');
title(tit);

if ~isempty(testcond)
    hold on;
    curr_rt=backup_rt;
    curr_rt(~(testcond==1))=nan; % test condition is true
    if all(isnan(curr_rt(1:end-1)-curr_rt(2:end)))
        return
    end
    [n,x]=histcounts(curr_rt(1:end-1)-curr_rt(2:end),x_backup);
    out.rt_change_testcond1=curr_rt(1:end-1)-curr_rt(2:end);
    [n,x]=cityscape_hist(n,x);
    plot(x,n./nansum(n),'Color','r');
    leg={'testcond FALSE','testcond TRUE'};
    legend(leg);
end

end