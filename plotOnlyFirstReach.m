function [reactionTimes,tbt,reactionTimes_preCue]=plotOnlyFirstReach(tbt,ds,whichReach,useAsCue,out,outfield,doPlot)

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
if settings.returnPreCueRTs==1
    firstreaches_precue=nan(1,size(temp,1));
else
    reactionTimes_preCue=[];
end
baselineWindows=zeros(size(temp,1),maind);
% note that everything has already been cue-aligned
for i=1:size(temp,1)
    fi=find(temp(i,maind+1:end)>0,1,'first')+maind;
    if settings.subtractBaselineReaching==1 || settings.divideByBaseRate==1 || settings.returnPreCueRTs==1
        % subtract off baseline reaching
        % first, calculate baseline reaching
        fi_base=find(temp(i,maind:-1:1)>0,1,'first');
        if ~isempty(fi_base)
            baselineWindows(i,fi_base)=1;
        end
    end
    if settings.divideByBaseRate==1
        if fi_base<fi
            fi=[];
        end
    end
    if settings.returnPreCueRTs==1
        if ~isempty(fi_base)
            fi_base=fi_base(1);
            firstreaches_precue(i)=fi_base;
        end
    end
    if ~isempty(fi)
        fi=fi(1);
        firstreaches(i)=fi;
    end
    if settings.longRT_ifNoReach==1 && isempty(fi)
        fi=size(temp,2);
        firstreaches(i)=fi;
    end
    if ~isempty(fi)
        newtemp(i,:)=zeros(size(newtemp(i,:)));
        newtemp(i,fi)=1;
    end
end
reactionTimes=(firstreaches-maind).*mode(diff(nanmean(tbt.times,1)));
if settings.returnPreCueRTs==1
    reactionTimes_preCue=firstreaches_precue.*mode(diff(nanmean(tbt.times,1)));
end
binToRTmapping=nan(1,size(temp,2));
binToRTmapping(maind+1:end)=[1:length(binToRTmapping(maind+1:end))].*mode(diff(nanmean(tbt.times,1)));
temp=newtemp;
tbt.firstReachesAfterCueOnset=newtemp;
if settings.subtractBaselineReaching==1
    % find time bins with reach # below baseline reach rate
    % zero out these reaction times
%     baseThresh=nanmean(nansum(baselineWindows(:,end-floor(1./mode(diff(nanmean(tbt.times,1)))):end),1)./size(baselineWindows,1));
    baseThresh=nanmean(nansum(baselineWindows,1)./size(baselineWindows,1));
    RTsummary=smooth(nansum(newtemp,1)./size(newtemp,1),floor(0.3/mode(diff(nanmean(tbt.times,1)))));
    whichBinsBelow=RTsummary<=baseThresh;
    if settings.throwOutAllAfter1stBaseline==1
        firstAtBase=find(whichBinsBelow(maind+1:end)==1,1,'first');
        whichBinsBelow(maind+firstAtBase:end)=1;
    end
    tbt.firstReachesAfterCueOnset(:,whichBinsBelow)=0;
    zeroTheseRT=binToRTmapping(whichBinsBelow);
    reactionTimes(ismember(reactionTimes,zeroTheseRT))=nan;
end
if settings.noRTlessThan~=0
    reactionTimes(reactionTimes<settings.noRTlessThan)=nan;
    binFirstRT=find(binToRTmapping<settings.noRTlessThan,1,'last');
    tbt.firstReachesAfterCueOnset(:,1:binFirstRT)=0;
end

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