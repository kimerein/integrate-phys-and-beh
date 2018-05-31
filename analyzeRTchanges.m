function analyzeRTchanges(alltbt,out,metadata)

nbins=200; % for histograms

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
plotCurrRT(curr_rt,'All trials, no reach during wheel turn (2 trials in a row)',out.led==1);
histRTchange(curr_rt,nbins,'All trials, no reach during wheel turn (2 trials in a row)',out.led==1);

% Change in reaction time, no preemptive reach
% Test effects of LED on PREVIOUS trial
condition=ones(length(out.led),1); % take all trials
curr_rt=reactionTimes;
curr_rt(~(condition==1) | dontUse_sequence)=nan;
plotCurrRT(curr_rt,'All trials, no reach during wheel turn (2 trials in a row)',out.led_1back==1);
histRTchange(curr_rt,nbins,'All trials, no reach during wheel turn (2 trials in a row)',out.led_1back==1);

% Change in reaction time, no preemptive reach
% Test effects of LED on PREVIOUS trial 
% given that mouse successfully touched pellet and performed cued reach on
% previous trial
condition=out.cued_reach_1back==1 & out.touched_pellet_1back==1;
curr_rt=reactionTimes;
curr_rt(~(condition==1) | dontUse_sequence)=nan;
plotCurrRT(curr_rt,'Cued reach & touched pellet 1 back & no preemptive reach for 2 trials in a row',out.led_1back==1);
histRTchange(curr_rt,nbins,'Cued reach & touched pellet 1 back & no preemptive reach for 2 trials in a row',out.led_1back==1);

% Change in reaction time, no preemptive reach
% Test effects of LED on PREVIOUS trial 
% given that mouse successfully touched pellet and performed cued reach on
% previous trial
condition=out.cued_reach_1back==1 & out.touched_pellet_1back==1;
curr_rt=reactionTimes;
curr_rt(~(condition==1))=nan;
plotCurrRT(curr_rt,'Cued reach & touched pellet 1 back',out.led_1back==1);
histRTchange(curr_rt,nbins,'Cued reach & touched pellet 1 back',out.led_1back==1);

% Change in reaction time, no preemptive reach
% Test effects of LED on PREVIOUS trial 
% given that mouse failed to touch pellet despite performing a cued reach on
% previous trial
condition=out.cued_reach_1back==1 & out.touched_pellet_1back==0;
curr_rt=reactionTimes;
curr_rt(~(condition==1) | dontUse_sequence)=nan;
plotCurrRT(curr_rt,'Cued reach but failed to touch pellet 1 back & no preemptive reach for 2 trials in a row',out.led_1back==1);
histRTchange(curr_rt,nbins,'Cued reach but failed to touch pellet 1 back & no preemptive reach for 2 trials in a row',out.led_1back==1);

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
scatter(curr_rt(1:end-1)+rand(size(curr_rt(1:end-1))).*0.01,curr_rt(2:end)+rand(size(curr_rt(2:end))).*0.01,[],'k','filled');
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
    scatter(curr_rt(1:end-1)+rand(size(curr_rt(1:end-1))).*0.01,curr_rt(2:end)+rand(size(curr_rt(2:end))).*0.01,[],'r','filled');
    leg={'testcond FALSE','testcond TRUE'};
    legend(leg);
end

end

function histRTchange(curr_rt,bins,tit,testcond)

if ~isempty(testcond)
    backup_rt=curr_rt;
    curr_rt(~(testcond==0))=nan; % test condition is false
end

if all(isnan(curr_rt(1:end-1)-curr_rt(2:end)))
    return
end
[n,x]=histcounts(curr_rt(1:end-1)-curr_rt(2:end),bins);
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
    [n,x]=cityscape_hist(n,x);
    plot(x,n./nansum(n),'Color','r');
    leg={'testcond FALSE','testcond TRUE'};
    legend(leg);
end

end