function reactionTimes=plotOnlyFirstReach(tbt,ds,whichReach,useAsCue,out,outfield,firstReachWindow)

excludePawOnWheelBeforeCue=1;
excludeReachesDuringCue=1;
cueDuration=0.2; % in seconds
cueDurationInds=floor(cueDuration/mode(diff(nanmean(tbt.times,1))));

useTheseTrials=ones(size(tbt.(whichReach),1),1);
if excludePawOnWheelBeforeCue==1
    useTheseTrials(out.pawOutDuringWheel==0)=0;
end

% cue ind
avcue=nanmean(tbt.(useAsCue),1);
[~,maind]=max(avcue);

% which reaches to plot
temp=tbt.(whichReach);
oldtemp=tbt.(whichReach);

% plot only first reaches?
firstreaches=nan(1,size(temp,1));
% note that everything has already been cue-aligned
for i=1:size(temp,1)
    % check if mouse reached during cue in this trial
    fi=find(temp(i,maind:maind+cueDurationInds)>0,1,'first');
    if excludeReachesDuringCue==1 && ~isempty(fi) % if yes
        % exclude this trial
        useTheseTrials(i)=0;
        continue
    end
    fi=find(temp(i,maind+1:end)>0,1,'first')+maind;
%     if length(fi)>1
    if ~isempty(fi)
        fi=fi(1);
        firstreaches(i)=fi;
    end
    temp(i,:)=zeros(size(temp(i,:)));
    if ~isempty(fi)
%         temp(i,fi)=oldtemp(i,fi);
        temp(i,fi)=1;
    end
end
[~,maind]=max(avcue);
reactionTimes=(firstreaches-maind).*mode(diff(nanmean(tbt.times,1)));

figure(); 
% plot(downSampAv(nanmean(tbt.times-0.2,1),ds),downSampAv(nanmean(tbt.(useAsCue),1),ds),'Color','b');
plot(downSampAv(nanmean(tbt.times,1),ds),downSampAv(nanmean(tbt.(useAsCue),1),ds),'Color','b');
hold on; 
plot(downSampAv(nanmean(tbt.times,1),ds),downSampAv(nansum(temp(out.(outfield)==0 & useTheseTrials==1,:),1)./sum(out.(outfield)==0 & useTheseTrials==1),ds),'Color','k');
plot(downSampAv(nanmean(tbt.times,1),ds),downSampAv(nansum(temp(out.(outfield)==1 & useTheseTrials==1,:),1)./sum(out.(outfield)==1 & useTheseTrials==1),ds),'Color','r');

disp(['Fraction of trials with first reach within first ' num2str(firstReachWindow) ' seconds after cue']);
u=unique(out.(outfield));
u=u(~isnan(u));
for i=1:length(u)
    disp(['For out.' outfield ' values equal to ' num2str(u(i))]);
    disp(sum(reactionTimes(out.(outfield)==u(i) & useTheseTrials==1)<firstReachWindow)./sum(out.(outfield)==u(i) & useTheseTrials==1));
end
