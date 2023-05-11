function [lowspeed_tbt,highspeed_tbt]=alignHighSpeedEventsToLowSpeedEvents(lowspeed_tbt,location_of_rig_events,fps,cue_duration,distractor_duration,alignment,nHighSpeedVidFiles,framesPerHighSpeedVid,nLowSpeedFrames)

dsby=5;
minTrialLength=9; % in seconds
distractorLength=0.25; % in seconds

% timestep on Kim's rig is usually 0.0039, Fs = 256 fps
timestep=1/fps;
ds_timestep=timestep*dsby;

load(location_of_rig_events);

% Guess at initial alignment
durationOfHighSpeedAcq=nHighSpeedVidFiles*framesPerHighSpeedVid*timestep;
temp=lowspeed_tbt.times_wrt_trial_start; temp=temp'; timestepls=mode(diff(temp(1:end)));
durationOfLowSpeedAcq=nLowSpeedFrames*timestepls;
disp(['duration of low speed acquisition: ' num2str(durationOfLowSpeedAcq) ' seconds']);
disp(['duration of high speed acquisition: ' num2str(durationOfHighSpeedAcq) ' seconds']);

% Make high speed tbt, downsampling 10 times
% This is just for initial alignment
% Then will precisely align to each cue
distractor=zeros(1,ceil(nanmax(distractorDiffEvs_minus/dsby)));
cue=zeros(1,ceil(nanmax(distractorDiffEvs_minus/dsby)));
wheel=zeros(1,ceil(nanmax(distractorDiffEvs_minus/dsby)));
reach=zeros(1,ceil(nanmax(distractorDiffEvs_minus/dsby)));

cuediff=floor(cueDiffEvs_plus./dsby);
distractordiff=floor(distractorDiffEvs_plus./dsby);
wheeldiff=floor(wheelDiffEvs_plus./dsby);
reachdiff=floor(reachDiffEvs_plus./dsby);

cue(cuediff)=1;
distractor(distractordiff)=1;
wheel(wheeldiff)=1;
reach(reachdiff)=1;

t=nanmean(lowspeed_tbt.times_wrt_trial_start,1);
[~,f]=nanmax(nanmean(lowspeed_tbt.cueZone_onVoff,1));
cuedelay=t(f);
cueindsbefore=floor(cuedelay/ds_timestep);
trialLength=nanmax(nanmean(lowspeed_tbt.times_wrt_trial_start,1));
% Use trialLength to clean up cues
f=find(cue>0.5);
triallengthinds=floor((minTrialLength/2)/ds_timestep);
for i=1:length(f)
    if cue(f)<0.5 
        continue
    end
    cue(f(i)+1:f(i)+triallengthinds)=0;
end
figure();
plot(cue); title('Cleaned up cue from high speed video');
cuediff=find(cue>0.5);

% Use distractorLength to clean up distractor
f=find(distractor>0.5);
distractorLengthinds=floor((2*distractorLength/3)/ds_timestep);
for i=1:length(f)
    if distractor(f)<0.5 
        continue
    end
    distractor(f(i)+1:f(i)+distractorLengthinds)=0;
end
figure();
plot(distractor); title('Cleaned up distractor from high speed video');
distractordiff=find(distractor>0.5);

% Start with big alignment
% fractionThrough=0.3;
% isInBackHalf=false;
% distract_thresh_movie=0.5;
% distract_thresh_paw=0.5;
% settings.photo_fs=floor(fps/dsby); 
% ard_ts=diff(alignment.timesfromarduino);
% ard_ts(ard_ts==0)=nan;
% ard_ts=mode(ard_ts); % is in ms
% settings.movie_fs=1/(ard_ts/1000); % movie data sampling rate in Hz
% ard_ts=ard_ts./1000;
% disp(['Check this. Movie data sampling rate is ' num2str(settings.movie_fs)]);
% pause;
% settings.scale_factor=floor(settings.photo_fs/settings.movie_fs);
% settings.photo_dec=1;
% settings.movie_dec=1;
% % discard this many frames from beginning
% settings.minlagForInitialAlign=[];
% settings.maxlagForInitialAlign=[]; % [] is don't want to constrain alignment
% settings.try_delay1=0;
% settings.try_delay2=0;
% if ~isempty(settings.minlagForInitialAlign) || ~isempty(settings.maxlagForInitialAlign)
%     questdlg('Preset min and max lag. Continue?');
% end

distractor=distractor-cue;
distractor(distractor<0)=0;

% %alignDistractors(alignment.movie_distractor,distractor,distract_thresh_movie,distract_thresh_paw,0:ard_ts:(length(alignment.movie_distractor)-1)*ard_ts,0:ds_timestep:(length(distractor)-1)*ds_timestep,settings,isInBackHalf,fractionThrough);
% alignDistractors(distractor+cue,alignment.movie_distractor+alignment.cueZone_onVoff,distract_thresh_paw,distract_thresh_movie,0:ds_timestep:(length(distractor)-1)*ds_timestep,0:ard_ts:(length(alignment.movie_distractor)-1)*ard_ts,settings,isInBackHalf,fractionThrough);
% % alignDistractors(distractor,alignment.movie_distractor,distract_thresh_paw,distract_thresh_movie,0:ds_timestep:(length(distractor)-1)*ds_timestep,0:ard_ts:(length(alignment.movie_distractor)-1)*ard_ts,settings,isInBackHalf,fractionThrough);

highspeed_tbt.cue=zeros(length(cuediff),floor(trialLength/ds_timestep));
highspeed_tbt.distractor=zeros(length(cuediff),floor(trialLength/ds_timestep));
highspeed_tbt.wheel=zeros(length(cuediff),floor(trialLength/ds_timestep));
highspeed_tbt.times=zeros(length(cuediff),floor(trialLength/ds_timestep));
highspeed_tbt.reach=zeros(length(cuediff),floor(trialLength/ds_timestep));
for i=1:length(cuediff)
    currcueind=cuediff(i);
    tempinds=currcueind-cueindsbefore:currcueind-cueindsbefore+size(highspeed_tbt.cue,2)-1;
    addnanatfront=0;
    addnanatback=0;
    if tempinds(1)<1
        addnanatfront=nansum(tempinds<1);
        tempinds=tempinds(find(tempinds>=1,1,'first'):end);
    end
    if tempinds(end)>length(cue)
        addnanatback=nansum(tempinds>length(cue));
        tempinds=tempinds(1:find(tempinds<=length(cue),1,'last'));
    end
    highspeed_tbt.cue(i,:)=[nan(1,addnanatfront) cue(tempinds) nan(1,addnanatback)];
    highspeed_tbt.distractor(i,:)=[nan(1,addnanatfront) distractor(tempinds) nan(1,addnanatback)];
    highspeed_tbt.wheel(i,:)=[nan(1,addnanatfront) wheel(tempinds) nan(1,addnanatback)]; 
    highspeed_tbt.reach(i,:)=[nan(1,addnanatfront) reach(tempinds) nan(1,addnanatback)]; 
    highspeed_tbt.times(i,:)=0:ds_timestep:(size(highspeed_tbt.times,2)-1)*ds_timestep;
end

% Fill out distractor and cue with true durations
highspeed_tbt=fleshOutDuration(highspeed_tbt,'cue',floor(cue_duration./ds_timestep));
highspeed_tbt=fleshOutDuration(highspeed_tbt,'distractor',floor(distractor_duration./ds_timestep));

%%%%%%%%%%% USER SETS THIS
lowspeed_tbt=moveRedForwardOrBack(lowspeed_tbt,[5:size(lowspeed_tbt.cue,1)],[],'drop');

% Plot behavior from high speed vs low speed video
plotBehFromLowSpeedMovie(lowspeed_tbt);
plotBehFromHighSpeedMovie(highspeed_tbt);

% useDistractorAlignment(lowspeed_tbt,'times_wrt_trial_start','movie_distractor',highspeed_tbt,'times','distractor','data2',false,'cueZone_onVoff','cue');

end

function tbt=fleshOutDuration(tbt,whichField,indsToFill)

temp=tbt.(whichField);
for i=1:size(temp,1)
    temprow=temp(i,:);
    f=find(temprow>0.5);
    for j=1:length(f)
        if f(j)+indsToFill>length(temprow)
            temprow(f(j):end)=1;
        else
            temprow(f(j):f(j)+indsToFill)=1;
        end
    end
    temp(i,:)=temprow;
end
tbt.(whichField)=temp;

end

function plotBehFromHighSpeedMovie(alltbt)

settings=plotCueTriggered_settings();
settings.plotfields={'cue','distractor','wheel','reach'};
settings.plotevents=settings.plotfields;
settings.eventOutlines={'b','y','k','g'};
settings.eventThresh={[0.5],[0.5],[0.5],[0.5]};
settings.eventColors={'b','y','k','g'};
settings.firstN={1,10,1,10};
settings.histoplotfields={'cue','reach'};
settings.shading_type=[];
plotBehavior(alltbt,'cue',false,1:size(alltbt.cue,1),settings);

end

function plotBehFromLowSpeedMovie(alltbt)

settings=plotCueTriggered_settings();
settings.plotfields={'movie_distractor','cueZone_onVoff','optoZone','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts','reachStarts'};
settings.plotevents=settings.plotfields;
settings.eventOutlines={'y','b','m','g','g','g','g','g'};
figure(); plot(alltbt.optoZone(1:20,:)');
temp=input('Thresh for optoZone: ',"s");
if strcmp(temp,'optoOn')
    temp=0.5;
    settings.plotfields{2}='optoOn';
    settings.plotevents{2}=settings.plotfields{2};
else
    temp=eval(temp);
end
settings.eventThresh={[0.5],[0.5],[temp],[0.5],[0.5],[0.5],[0.5],[0.5]};
settings.eventColors={'y','b','none','g','g','g','g','g'};
settings.firstN={'all',1,[5],'all','all','all','all','all'};
settings.histoplotfields={'cueZone_onVoff','all_reachBatch'};
settings.shading_type=[];
plotBehavior(alltbt,'cueZone_onVoff',false,1:size(alltbt.cueZone_onVoff,1),settings);

end

function [discardedPhotoFrames_time,frontShift_time,scaleBy,movie_LED,movie_times,scaleMovieTimes,addToMovieTimes,padPhotoTimesAtFront]=alignDistractors(movie_distract,photo_distract,distract_thresh_movie,distract_thresh_photometry,movie_times,photo_times,settings,isInBackHalf,fractionThrough)

% Get when LED was on in movie vs off
movie_LED=single(movie_distract>distract_thresh_movie);

% Find best alignment of distractor LED in movie and photometry output -- note
% different sampling rates
photo_LED=single(photo_distract>distract_thresh_photometry);

% Find alignment
% First down-sample photometry LED
photo_dec=settings.photo_dec;
movie_dec=settings.movie_dec;

photo_LED=decimate(double(photo_LED),photo_dec);
photo_times=decimate(double(photo_times),photo_dec);

movie_LED=decimate(double(movie_LED),movie_dec);
movie_times=decimate(double(movie_times),movie_dec);

% Do an initial alignment
temp=photo_LED;
photo_LED(temp>=0.5)=1;
photo_LED(temp<0.5)=0;
temp=movie_LED;
movie_LED(temp>=0.5)=1;
movie_LED(temp<0.5)=0;

% for purposes of alignment, discard photometry before first distractor
% onset
[pks_arduino,locs_arduino]=findpeaks(photo_LED);
discardedPhotoFrames=locs_arduino(1)-1;
discardedPhotoFrames_time=mode(photo_times(2:end)-photo_times(1:end-1))*discardedPhotoFrames; % in seconds
photo_LED=photo_LED(locs_arduino(1)-1:end);
arduino_LED_ITIs=diff(photo_times(locs_arduino));
[pks,locs]=findpeaks(movie_LED);
movie_LED_ITIs=diff(movie_times(locs)); 

photo_times=photo_times(locs_arduino(1)-1:end);

if isInBackHalf==true
    backup_arduino_LED_ITIs=arduino_LED_ITIs;
    midLength=ceil(fractionThrough*length(arduino_LED_ITIs));
    arduino_LED_ITIs=arduino_LED_ITIs(midLength+1:end);
end

% temp1=arduino_LED_ITIs./max(arduino_LED_ITIs);
% temp2=movie_LED_ITIs./max(movie_LED_ITIs);
temp1=arduino_LED_ITIs;
temp2=movie_LED_ITIs;

if ~isempty(settings.minlagForInitialAlign)
    temp2=[nan(1,settings.minlagForInitialAlign) temp2];
else
    settings.minlagForInitialAlign=0;
end

if isempty(settings.maxlagForInitialAlign) 
    [X,Y,D]=alignsignals(temp1,temp2); 
%     [X,Y,D]=alignsignals(photo_LED,movie_LED); 
else
    [X,Y,D]=alignsignals(temp1,temp2,settings.maxlagForInitialAlign); 
%     [X,Y,D]=alignsignals(photo_LED,movie_LED,settings.maxlagForInitialAlign); 
end

if isInBackHalf==true
    D=D-midLength;
    arduino_LED_ITIs=backup_arduino_LED_ITIs;
    if D>0
%         X=[zeros(1,D) arduino_LED_ITIs./max(arduino_LED_ITIs)];
        X=[zeros(1,D) arduino_LED_ITIs];
%         Y=movie_LED_ITIs./max(movie_LED_ITIs);
        Y=movie_LED_ITIs;
    elseif D<0
%         Y=[zeros(1,-D) movie_LED_ITIs./max(movie_LED_ITIs)];
        Y=[zeros(1,-D) movie_LED_ITIs];
%         X=arduino_LED_ITIs./max(arduino_LED_ITIs);
        X=arduino_LED_ITIs;
    end
end

if D>0
    % photometry starts after movie
    % pad photometry at front
    figure();
    plot(X,'Color','b'); title('Blue is paw distractor, red is movie distractor');
    hold on;
    plot(Y,'Color','r');
    anchor=input('Choose index for alignment anchor. ');
    alignAtPeak_arduino=locs_arduino(anchor-D);
    alignAtPeak_movie=locs(anchor-settings.minlagForInitialAlign);
    photo_LED=[nan(1,(-alignAtPeak_arduino+alignAtPeak_movie)) photo_LED];
    photo_times=[nan(1,(-alignAtPeak_arduino+alignAtPeak_movie)) photo_times];
    tempphot=photo_times(2:end)-photo_times(1:end-1);
    padPhotoTimesAtFront=mode(tempphot(~isnan(tempphot)))*D;
    guess_best_delay=0;
else
    padPhotoTimesAtFront=0;
    figure();
    plot(X,'Color','b');
    hold on;
    plot(Y,'Color','r');
    title('Preliminary alignment of photometry distractor intervals onto arduino distractor intervals');
    legend({'Photometry distractor intervals','Movie distractor intervals'});
    
    anchor=input('Choose index for alignment anchor. ');
    % align at peak
    alignAtPeak_arduino=locs_arduino(anchor);
    alignAtPeak_movie=locs(anchor+D-settings.minlagForInitialAlign);
    % pad movie to align at this peak
    figure();
    plot([nan(1,(alignAtPeak_arduino-alignAtPeak_movie)) movie_LED],'Color','r');
    hold on;
    plot(photo_LED,'Color','b');

    guess_best_delay=alignAtPeak_arduino-alignAtPeak_movie;
    if guess_best_delay<0
        guess_best_delay=0;
        disp('WARNING: Movie appears to start before photometry');
    end
end
trydelays=guess_best_delay+settings.try_delay1:guess_best_delay+settings.try_delay2;
% Note that fixed, so now best scale is 1
guess_best_scale=1;


figure();
plot([photo_LED(guess_best_delay+1:end) nan(1,guess_best_delay)],'Color','b');
hold on;
plot(movie_LED,'Color','r');
title('Preliminary alignment of movie distractor onto photometry distractor');
legend({'Photometry distractor','Movie distractor'});

if isInBackHalf==true
    trydelays=guess_best_delay;
end

% Test signal alignment and scaling
trydelays=trydelays(trydelays>=0);
disp('Now refining alignment ...');
sumdiffs=nan(1,length(trydelays));
for i=1:length(trydelays)
    currdelay=trydelays(i);
    if mod(i,500)==0
        % disp(i);
    end
    temp_movie=movie_LED;
    temp_arduino=[photo_LED(currdelay+1:end-currdelay) nan(1,currdelay)];
    if length(temp_arduino)>length(temp_movie)
        temp_movie=[temp_movie nan(1,length(temp_arduino)-length(temp_movie))];
        sumdiffs(i)=nansum(abs(temp_movie-temp_arduino));
    elseif length(temp_movie)>length(temp_arduino)
        temp_movie=temp_movie(1:length(temp_arduino));
    else
        sumdiffs(i)=nansum(abs(temp_movie-temp_arduino))./nansum(~isnan(temp_arduino));
    end
end
sumdiffs(isinf(sumdiffs))=nan;
ma=max(sumdiffs(:));
sumdiffs(isnan(sumdiffs))=3*ma;
[minval,mi]=min(sumdiffs(:));

figure(); 
imagesc(sumdiffs);
title('Finding best alignment');
xlabel('Trying different delays');
disp('Best delay');
disp(trydelays(mi));
guess_best_delay=trydelays(mi);

frontShift=guess_best_delay;
tempphot=photo_times(2:end)-photo_times(1:end-1);
frontShift_time=mode(tempphot(~isnan(tempphot)))*frontShift; % shift movie back by this amount, i.e., shift photometry forward by this amount
scaleBy=1;
photo_LED=[photo_LED(frontShift+1:end) nan(1,frontShift)];
photo_times=[photo_times(frontShift+1:end) nan(1,frontShift)];
movie_times=movie_times/1000; % convert to secs instead of ms
fir=find(~isnan(photo_times),1,'first');
figure(); plot(photo_times,photo_LED,'Color','k'); hold on; plot(movie_times+(photo_times(fir)-movie_times(fir)),movie_LED,'Color','r'); title('Shifted alignment movie to photometry distractors');
figure(); plot(photo_LED,'Color','k'); hold on; plot(movie_LED,'Color','r');

%disp(['scale photometry by ' num2str(scaleBy)]);
% scale movie_times subtly to match photometry
g1=input('Enter starting index / alignment anchor for black (photometry) distractors. ');
g2=input('Enter ending index for black (photometry) distractors. ');
g3=input('Enter starting index / alignment anchor for red (movie) distractors. ');
g4=input('Enter ending index for red (movie) distractors. ');

timestep_photometry=photo_times(g2)-photo_times(g1);
timestep_movie=movie_times(g4)-movie_times(g3); % convert to s
movie_times=movie_times*(timestep_photometry/timestep_movie);
scaleMovieTimes=(timestep_photometry/timestep_movie);
movie_times=movie_times+(photo_times(g1)-movie_times(g3));
addToMovieTimes=photo_times(g1)-movie_times(g3);

figure(); plot(photo_times,photo_LED,'Color','k'); hold on; plot(movie_times+(photo_times(fir)-movie_times(fir)),movie_LED,'Color','r'); 
title('Alignment after movie time scaling');
xlabel('Time (s)');
ylabel('Distractor');

end