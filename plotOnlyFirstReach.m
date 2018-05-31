function [reactionTimes,tbt]=plotOnlyFirstReach(tbt,ds,whichReach,useAsCue,out,outfield,doPlot)

settings=RTanalysis_settings();
excludePawOnWheelDuringCue=settings.excludePawOnWheelDuringCue;

% cue ind
avcue=nanmean(tbt.(useAsCue),1);
[~,maind]=max(avcue);

% which reaches to get
temp=tbt.(whichReach);
newtemp=nan(size(tbt.(whichReach)));

% get only first reaches
firstreaches=nan(1,size(temp,1));
% note that everything has already been cue-aligned
for i=1:size(temp,1)
    fi=find(temp(i,maind+1:end)>0,1,'first')+maind;
    if ~isempty(fi)
        fi=fi(1);
        firstreaches(i)=fi;
    end
    if ~isempty(fi)
        newtemp(i,:)=zeros(size(newtemp(i,:)));
        newtemp(i,fi)=1;
    end
end
reactionTimes=(firstreaches-maind).*mode(diff(nanmean(tbt.times,1)));
temp=newtemp;
tbt.firstReachesAfterCueOnset=newtemp;

if doPlot==1
    if excludePawOnWheelDuringCue==1
        useTrials=out.paw_during_wheel==0;
    else
        useTrials=ones(length(out.paw_during_wheel),1);
    end
    figure();
    plot(downSampAv(nanmean(tbt.times,1),ds),downSampAv(nanmean(tbt.(useAsCue),1),ds),'Color','b');
    hold on;
    noreachtrial=all(isnan(temp),2);
    if ischar(outfield)
        plot(downSampAv(nanmean(tbt.times,1),ds),downSampAv(nansum(temp(out.(outfield)==0 & useTrials==1 & noreachtrial==0,:),1)./sum(out.(outfield)==0 & useTrials==1 & noreachtrial==0),ds),'Color','k');
        plot(downSampAv(nanmean(tbt.times,1),ds),downSampAv(nansum(temp(out.(outfield)==1 & useTrials==1 & noreachtrial==0,:),1)./sum(out.(outfield)==1 & useTrials==1 & noreachtrial==0),ds),'Color','r');
        leg={'cue','out.field FALSE','out.field TRUE'};
    elseif islogical(outfield)
        plot(downSampAv(nanmean(tbt.times,1),ds),downSampAv(nansum(temp(outfield==false & useTrials==1 & noreachtrial==0,:),1)./sum(outfield==false & useTrials==1 & noreachtrial==0),ds),'Color','k');
        plot(downSampAv(nanmean(tbt.times,1),ds),downSampAv(nansum(temp(outfield==true & useTrials==1 & noreachtrial==0,:),1)./sum(outfield==true & useTrials==1 & noreachtrial==0),ds),'Color','r');
        leg={'cue','outfield FALSE','outfield TRUE'};
    elseif isnumeric(outfield)
        plot(downSampAv(nanmean(tbt.times,1),ds),downSampAv(nansum(temp(outfield==0 & useTrials==1 & noreachtrial==0,:),1)./sum(outfield==0 & useTrials==1 & noreachtrial==0),ds),'Color','k');
        plot(downSampAv(nanmean(tbt.times,1),ds),downSampAv(nansum(temp(outfield==1 & useTrials==1 & noreachtrial==0,:),1)./sum(outfield==1 & useTrials==1 & noreachtrial==0),ds),'Color','r');
        leg={'cue','outfield 0','outfield 1'};
    end
    xlabel('Time (seconds)');
    ylabel('# reaches / # trials');
    legend(leg);
end