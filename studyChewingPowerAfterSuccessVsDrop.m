function [out1,out2]=studyChewingPowerAfterSuccessVsDrop(tbt,savehandles,zoneVals,eat,successReachName,dropReachName,overweightFP)

movieframesEarlyChews=savehandles.discardFirstNFrames+[1:length(eat.isChewing)];
settings=autoReachAnalysisSettings;
params.Fs=settings.movie_fps;
params.tapers=settings.chew.tapers;
params.fpass=settings.chew.fpass; % in Hz
tempie=zoneVals.eatZone;
eatzone=tempie(savehandles.discardFirstNFrames+1:end);
[S,t,f]=mtspecgramc(eatzone(~isnan(eatzone)),[5 0.25],params);
frameTimes=0:(1/settings.movie_fps):(length(movieframesEarlyChews(~isnan(movieframesEarlyChews)))-1)*(1/settings.movie_fps);

% ntimes_subsequent=1.667; % in sec
ntimes_subsequent=1.25; % in sec
ninds_subsequent=floor(ntimes_subsequent/(1/params.Fs));
[out1,bestThresh,bestThreshIntens]=getPowerDurationInWindow(tbt,successReachName,dropReachName,movieframesEarlyChews,t,frameTimes,eat,ninds_subsequent,[],[],[],[],overweightFP,eatzone,true);
if isempty(bestThresh) && isempty(bestThreshIntens)
    out1.threshold=['chewingPower>' num2str(-100000) ' & rawIntensity>' num2str(-100000)];
elseif isempty(bestThresh) && ~isempty(bestThreshIntens)
    out1.threshold=['chewingPower>' num2str(-100000) ' & rawIntensity>' num2str(bestThreshIntens)];
elseif ~isempty(bestThresh) && isempty(bestThreshIntens)
    out1.threshold=['chewingPower>' num2str(bestThresh) ' & rawIntensity>' num2str(-100000)];
else
    out1.threshold=['chewingPower>' num2str(bestThresh) ' & rawIntensity>' num2str(bestThreshIntens)];
end

% ntimes_subsequent=6; % in sec
% ninds_subsequent=floor(ntimes_subsequent/(1/params.Fs));
% out2=getPowerDurationInWindow(tbt,successReachName,dropReachName,movieframesEarlyChews,t,frameTimes,eat,ninds_subsequent,[],(6/20)*ninds_subsequent,[],[]);
% out2.threshold=['chewingDuration>' num2str((6/20)*ninds_subsequent)];

ntimes_subsequent=4; % in sec
delaytime=2; % in sec
ninds_subsequent=floor(ntimes_subsequent/(1/params.Fs));
delay=floor(delaytime/(1/params.Fs));
[out2,~,~,bestThreshDur]=getPowerDurationInWindow(tbt,successReachName,dropReachName,movieframesEarlyChews,t,frameTimes,eat,ninds_subsequent,[],[],delay,delaytime,overweightFP,eatzone,true);
if ~isempty(bestThreshDur)
    out2.threshold=['chewingDuration>' num2str(bestThreshDur)];
else
    out2.threshold=['chewingDuration>' num2str(-100000)];
end
% out2=getPowerDurationInWindow(tbt,successReachName,dropReachName,movieframesEarlyChews,t,frameTimes,eat,ninds_subsequent,[],(6/20)*ninds_subsequent,delay,delaytime,overweightFP,eatzone,false);
% out2.threshold=['chewingDuration>' num2str((6/20)*ninds_subsequent)];

end

function [out,bestThresh,bestThreshIntens,bestThreshDur]=getPowerDurationInWindow(tbt,successReachName,dropReachName,movieframesEarlyChews,t,frameTimes,eat,ninds_subsequent,x_thresh,y_thresh,delay,delayTime,overweightFP,eatZone,showROCthresh)

temp=tbt.(successReachName);
temp=temp(1:end);
temp=temp';
fi=find(temp>0.5);
temp=tbt.movieframeinds(1:end);
linearmovieframes=temp';
subsequent_chewingPower=nan(1,length(fi));
subsequent_chewingDuration=nan(1,length(fi));
rawIntensity=nan(1,length(fi));
% figure(); plot(t,eat.chewingpower); title('chewing power');
out.chewingPower=[];
out.chewingDuration=[];
out.movieFrameInds=[];
out.rawIntensity=[];
for i=1:length(fi)
    currmovieind=linearmovieframes(fi(i));
    [~,mi]=nanmin(abs(movieframesEarlyChews-currmovieind));
    [~,closest_t]=nanmin(abs(t-frameTimes(mi)));
    if ~isempty(delay)
        [~,closest_t]=nanmin(abs(t-(frameTimes(mi)+delayTime)));
        mi=mi+delay;
    end
    if closest_t+ninds_subsequent>length(eat.chewingpower)
        subsequent_chewingPower(i)=nanmean(eat.chewingpower(closest_t:end));
    else
        subsequent_chewingPower(i)=nanmean(eat.chewingpower(closest_t:closest_t+ninds_subsequent));
    end
    if mi+ninds_subsequent>length(eat.isChewing)
        subsequent_chewingDuration(i)=nansum(eat.isChewing(mi:end))*(ninds_subsequent/length(eat.isChewing(mi:end)));
        rawIntensity(i)=nanmean(eatZone(mi:end));
    else
        subsequent_chewingDuration(i)=nansum(eat.isChewing(mi:mi+ninds_subsequent));
        rawIntensity(i)=nanmean(eatZone(mi:mi+ninds_subsequent));
    end
    if currmovieind>34130-5 && currmovieind<34130+5
        disp(subsequent_chewingDuration(i));
        disp(subsequent_chewingPower(i));
    end
    out.chewingPower(i)=subsequent_chewingPower(i);
    out.chewingDuration(i)=subsequent_chewingDuration(i);
    out.rawIntensity(i)=rawIntensity(i);
    out.movieFrameInds(i)=currmovieind;
end
f=figure(); 
scatter(subsequent_chewingPower,subsequent_chewingDuration,[],'g');
if ~isempty(x_thresh)
    mi=nanmin(subsequent_chewingDuration);
    ma=nanmax(subsequent_chewingDuration);
elseif ~isempty(y_thresh)
    mi=nanmin(subsequent_chewingPower);
    ma=nanmax(subsequent_chewingPower);
end
hold on;

temp=tbt.(dropReachName);
temp=temp(1:end);
temp=temp';
fi=find(temp>0.5);
subsequent_chewingPower=nan(1,length(fi));
subsequent_chewingDuration=nan(1,length(fi));
rawIntensity=nan(1,length(fi));
for i=1:length(fi)
    currmovieind=linearmovieframes(fi(i));
    [~,mi]=nanmin(abs(movieframesEarlyChews-currmovieind));
    [~,closest_t]=nanmin(abs(t-frameTimes(mi)));
    if ~isempty(delay)
        [~,closest_t]=nanmin(abs(t-(frameTimes(mi)+delayTime)));
        mi=mi+delay;
    end
    if closest_t+ninds_subsequent>length(eat.chewingpower)
        subsequent_chewingPower(i)=nanmean(eat.chewingpower(closest_t:end));
    else
        subsequent_chewingPower(i)=nanmean(eat.chewingpower(closest_t:closest_t+ninds_subsequent));
    end
    if mi+ninds_subsequent>length(eat.isChewing)
        subsequent_chewingDuration(i)=nansum(eat.isChewing(mi:end))*(ninds_subsequent/length(eat.isChewing(mi:end)));
        rawIntensity(i)=nanmean(eatZone(mi:end));
    else
        subsequent_chewingDuration(i)=nansum(eat.isChewing(mi:mi+ninds_subsequent));
        rawIntensity(i)=nanmean(eatZone(mi:mi+ninds_subsequent));
    end
    if currmovieind>34130-10 && currmovieind<34130+10
        disp(subsequent_chewingDuration(i));
        disp(subsequent_chewingPower(i));
    end
end
scatter(subsequent_chewingPower,subsequent_chewingDuration,[],'r');
xlabel('Chewing power');
ylabel('Chewing duration');
if ~isempty(x_thresh)
    line([x_thresh x_thresh],[nanmin([subsequent_chewingDuration mi]) nanmax([subsequent_chewingDuration ma])],'Color','k');
elseif ~isempty(y_thresh)
    line([nanmin([subsequent_chewingPower mi]) nanmax([subsequent_chewingPower ma])],[y_thresh y_thresh],'Color','k');
end
bestThresh=buildROC([out.chewingPower subsequent_chewingPower],[ones(size(out.chewingPower)) zeros(size(subsequent_chewingPower))],overweightFP);
if showROCthresh==true && ~isempty(bestThresh)
    set(0,'CurrentFigure',f);
    temp=[out.chewingDuration subsequent_chewingDuration];
    line([bestThresh bestThresh],[nanmin(temp)-1 nanmax(temp)+1],'Color','k');
end
bestThreshDur=buildROC([out.chewingDuration subsequent_chewingDuration],[ones(size(out.chewingDuration)) zeros(size(subsequent_chewingDuration))],overweightFP);
if showROCthresh==true && ~isempty(bestThreshDur)
    set(0,'CurrentFigure',f);
    temp=[out.chewingPower subsequent_chewingPower];
    line([nanmin(temp)-1 nanmax(temp)+1],[bestThreshDur bestThreshDur],'Color','k');
end

f=figure(); 
scatter(out.chewingPower,out.rawIntensity,[],'g');
hold on; scatter(subsequent_chewingPower,rawIntensity,[],'r');
xlabel('Chewing power');
ylabel('Raw intensity in eat zone');
bestThreshIntens=buildROC([out.rawIntensity rawIntensity],[ones(size(out.chewingPower)) zeros(size(subsequent_chewingPower))],overweightFP);
if showROCthresh==true && ~isempty(bestThreshIntens)
    set(0,'CurrentFigure',f);
    temp=[out.chewingPower subsequent_chewingPower];
    line([nanmin(temp)-1 nanmax(temp)+1],[bestThreshIntens bestThreshIntens],'Color','k');
end



end

function bestThresh=buildROC(pred,isSuccess,overweightFP)

% try thresholds
trythresh=nanmin(pred)-1:0.1:nanmax(pred)+1;
tpr=nan(1,length(trythresh));
fpr=nan(1,length(trythresh));
for i=1:length(trythresh)
    guessSuccess=pred>trythresh(i);
    tpr(i)=sum(isSuccess==1 & guessSuccess==1,'all','omitnan');
    fpr(i)=sum(isSuccess==0 & guessSuccess==1,'all','omitnan');
end
figure();
scatter(fpr,tpr,[],'k');
hold on;
xlabel('False positive rate');
ylabel('True positive rate');
% maximize TPR - FPR
[~,ma]=max(tpr-overweightFP*fpr,[],2,'omitnan');
bestThresh=trythresh(ma);
scatter(fpr(ma),tpr(ma),[],'r');
% check whether this threshold actually provides any value
if tpr(ma)/nansum(isSuccess==1) < 0.8 || fpr(ma)/nansum(isSuccess==0) > 0.2
    % not helpful
    bestThresh=[];
end

end