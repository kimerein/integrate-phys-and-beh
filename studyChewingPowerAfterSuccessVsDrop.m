function [out1,out2]=studyChewingPowerAfterSuccessVsDrop(tbt,savehandles,zoneVals,eat,successReachName,dropReachName)

movieframesEarlyChews=savehandles.discardFirstNFrames+[1:length(eat.isChewing)];
settings=autoReachAnalysisSettings;
params.Fs=settings.movie_fps;
params.tapers=settings.chew.tapers;
params.fpass=settings.chew.fpass; % in Hz
tempie=zoneVals.eatZone;
tempie=tempie(savehandles.discardFirstNFrames+1:end);
[S,t,f]=mtspecgramc(tempie(~isnan(tempie)),[5 0.25],params);
frameTimes=0:(1/settings.movie_fps):(length(movieframesEarlyChews(~isnan(movieframesEarlyChews)))-1)*(1/settings.movie_fps);

ninds_subsequent=50;
out1=getPowerDurationInWindow(tbt,successReachName,dropReachName,movieframesEarlyChews,t,frameTimes,eat,ninds_subsequent,3,[]);
out1.threshold='chewingPower>3';

ninds_subsequent=180;
out2=getPowerDurationInWindow(tbt,successReachName,dropReachName,movieframesEarlyChews,t,frameTimes,eat,ninds_subsequent,[],(6/20)*ninds_subsequent);
out2.threshold=['chewingDuration>' num2str((6/20)*ninds_subsequent)];


function out=getPowerDurationInWindow(tbt,successReachName,dropReachName,movieframesEarlyChews,t,frameTimes,eat,ninds_subsequent,x_thresh,y_thresh)

temp=tbt.(successReachName);
temp=temp(1:end);
temp=temp';
fi=find(temp>0.5);
temp=tbt.movieframeinds(1:end);
linearmovieframes=temp';
subsequent_chewingPower=nan(1,length(fi));
subsequent_chewingDuration=nan(1,length(fi));
% figure(); plot(t,eat.chewingpower); title('chewing power');
out.chewingPower=[];
out.chewingDuration=[];
out.movieFrameInds=[];
for i=1:length(fi)
    currmovieind=linearmovieframes(fi(i));
    [~,mi]=nanmin(abs(movieframesEarlyChews-currmovieind));
    [~,closest_t]=nanmin(abs(t-frameTimes(mi)));
    if closest_t+ninds_subsequent>length(eat.chewingpower)
        subsequent_chewingPower(i)=nanmean(eat.chewingpower(closest_t:end));
    else
        subsequent_chewingPower(i)=nanmean(eat.chewingpower(closest_t:closest_t+ninds_subsequent));
    end
    if mi+ninds_subsequent>length(eat.isChewing)
        subsequent_chewingDuration(i)=nansum(eat.isChewing(mi:end))*(ninds_subsequent/length(eat.isChewing(mi:end)));
    else
        subsequent_chewingDuration(i)=nansum(eat.isChewing(mi:mi+ninds_subsequent));
    end
    out.chewingPower(i)=subsequent_chewingPower(i);
    out.chewingDuration(i)=subsequent_chewingDuration(i);
    out.movieFrameInds(i)=currmovieind;
end
figure(); 
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
for i=1:length(fi)
    currmovieind=linearmovieframes(fi(i));
    [~,mi]=nanmin(abs(movieframesEarlyChews-currmovieind));
    [~,closest_t]=nanmin(abs(t-frameTimes(mi)));
    if closest_t+ninds_subsequent>length(eat.chewingpower)
        subsequent_chewingPower(i)=nanmean(eat.chewingpower(closest_t:end));
    else
        subsequent_chewingPower(i)=nanmean(eat.chewingpower(closest_t:closest_t+ninds_subsequent));
    end
    if mi+ninds_subsequent>length(eat.isChewing)
        subsequent_chewingDuration(i)=nansum(eat.isChewing(mi:end))*(ninds_subsequent/length(eat.isChewing(mi:end)));
    else
        subsequent_chewingDuration(i)=nansum(eat.isChewing(mi:mi+ninds_subsequent));
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