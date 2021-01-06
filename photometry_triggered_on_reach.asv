function data=photometry_triggered_on_reach(datadir_forPhotometry,datadir_forBehavior)

% rig-specific settings
% distractor threshold
distract_thresh_photometry=0.5;

% other settings
distract_thresh_movie=0.5; % from alignment

% load and process photometry data
data=processPhotometry(datadir_forPhotometry);

% load already processed behavior data


end

function alignDistractors(movie_distract,photo_distract,distract_thresh_movie,distract_thresh_photometry,movie_times,photo_times)

% Get when LED was on in movie vs off
movie_LED=single(movie_distract>distract_thresh_movie);
movie_LED_for_finalalignment=movie_LED;

% Find best alignment of distractor LED in movie and photometry output -- note
% different sampling rates
temptimes=[];
for i=1:size(out.allTrialTimes,1)
    temptimes=[temptimes out.allTrialTimes(i,:)];
end
temp=temp(1:end);
photo_LED=single(photo_distract>distract_thresh_photometry);
arduino_times=temptimes(~isnan(temptimes));
backup_arduino_times=arduino_times;

% Find alignment
% First down-sample arduino LED
arduino_dec=settings.arduino_dec;
movie_dec=settings.movie_dec;

arduino_LED=decimate(arduino_LED,arduino_dec);
arduino_times=decimate(arduino_times,arduino_dec);

testRun_movieLED=double(movie_LED);

movie_LED=decimate(double(movie_LED),movie_dec);
movie_times=decimate(movie_times,movie_dec);

% Do an initial alignment
temp=arduino_LED;
arduino_LED(temp>=0.5)=1;
arduino_LED(temp<0.5)=0;
temp=movie_LED;
movie_LED(temp>=0.5)=1;
movie_LED(temp<0.5)=0;

% Discard beginning of Arduino LED
arduino_LED(arduino_times<settings.discardTimeArduinoLED*1000)=nan; % convert discardTimeArduinoLED to ms

% Throw out LED distractor on intervals less than settings.useDistractorThresh
% This deals with skipping of low frame rate DVR
allEvents_arduino_LED=arduino_LED;
allEvents_movie_LED=movie_LED;
arduino_LED=throwOutOnStretches(arduino_LED,arduino_times);
[movie_LED,throwOutMovie]=throwOutOnStretches(movie_LED,movie_times);

[pks_arduino,locs_arduino]=findpeaks(arduino_LED);
arduino_LED_ITIs=diff(arduino_times(locs_arduino));
[pks,locs]=findpeaks(movie_LED);
movie_LED_ITIs=diff(movie_times(locs)); 

if settings.isInSecondHalf==true
    backup_arduino_LED_ITIs=arduino_LED_ITIs;
    midLength=floor(settings.fractionThroughArduino*length(arduino_LED_ITIs))+20;
    arduino_LED_ITIs=arduino_LED_ITIs(midLength+1:end);
end
temp1=arduino_LED_ITIs./max(arduino_LED_ITIs);
temp2=movie_LED_ITIs./max(movie_LED_ITIs);
% temp1=temp1-nanmean(temp1);
% temp2=temp2-nanmean(temp2);
if isempty(settings.maxlagForInitialAlign) 
    [X,Y,D]=alignsignals(temp1,temp2); 
else
    [X,Y,D]=alignsignals(temp1,temp2,settings.maxlagForInitialAlign); 
end
if settings.isInSecondHalf==true
    D=D-midLength;
    arduino_LED_ITIs=backup_arduino_LED_ITIs;
    if D>0
        X=[zeros(1,D) arduino_LED_ITIs./max(arduino_LED_ITIs)];
        Y=movie_LED_ITIs./max(movie_LED_ITIs);
    elseif D<0
        Y=[zeros(1,-D) movie_LED_ITIs./max(movie_LED_ITIs)];
        X=arduino_LED_ITIs./max(arduino_LED_ITIs);
    end
end

tryinc=settings.tryinc; % this is the increment for trying different scalings of movie onto arduino data
if D>0
    error('Why does movie start before Arduino?');
else
    movie_peakloc=1;
    arduino_peakloc=abs(D)+1;
    movie_peak_indexIntoMovie=locs(movie_peakloc);
    arduino_peak_indexIntoArduino=locs_arduino(arduino_peakloc);
    if (-D)+1+length(movie_LED_ITIs)>length(locs_arduino)
        size_of_arduino=length(arduino_LED(locs_arduino((-D)+1):end));
    else
        size_of_arduino=length(arduino_LED(locs_arduino((-D)+1):locs_arduino((-D)+1+length(movie_LED_ITIs))));
    end
    size_of_movie=length(movie_LED(locs(1):locs(end)));
%     size_of_movie=length(movie_LED(locs(383-312):locs(403-312)));
%     size_of_arduino=length(arduino_LED(locs_arduino(383):locs_arduino(403)));
%     movie_peak_indexIntoMovie=locs(383-312);
%     arduino_peak_indexIntoArduino=locs_arduino(383);
    guess_best_scale=size_of_arduino/size_of_movie;
    guess_best_scale1=guess_best_scale;
    % Adjust according to guess_best_scale
    movledinds=1:length(movie_LED);
    movie_LED=resample(movie_LED,floor(mod(size_of_arduino/size_of_movie,1)*100)+floor((guess_best_scale*100)/100)*100,100);
    movledinds=resample(movledinds,floor(mod(size_of_arduino/size_of_movie,1)*100)+floor((guess_best_scale*100)/100)*100,100);
    [~,mi]=min(abs(movledinds-movie_peak_indexIntoMovie));
    movie_peak_indexIntoMovie=mi;
    guess_best_delay=arduino_peak_indexIntoArduino-movie_peak_indexIntoMovie;
    trydelays=guess_best_delay+settings.try_delay1:guess_best_delay+settings.try_delay2;
    % Note that fixed, so now best scale is 1
    guess_best_scale=1;
    tryscales=guess_best_scale+settings.try_scale1:tryinc:guess_best_scale+settings.try_scale2;
end

figz(2)=figure();
plot(X,'Color','b');
hold on;
plot(Y,'Color','r');
title('Preliminary alignment of movie distractor intervals onto arduino distractor intervals');
legend({'Arduino distractor intervals','Movie distractor intervals'});

figz(3)=figure();
plot(arduino_LED,'Color','b');
hold on;
plot([nan(1,guess_best_delay) movie_LED],'Color','r');
title('Preliminary alignment of movie distractor onto arduino distractor');
legend({'Arduino distractor','Movie distractor'});

% Wait for user to confirm preliminary alignment
if settings.isOrchestra~=1
    pause;
end
 
% Test signal alignment and scaling
disp('Now refining alignment ...');
sumdiffs=nan(length(tryscales),length(trydelays));
if settings.alignWithAllEvents==1
    backup_movie_LED=allEvents_movie_LED;
    forSecondaryAlignment=backup_movie_LED;
    backup_arduino_LED=allEvents_arduino_LED;
    % Remove LED distractor on intervals that are too short (i.e., may have
    % been missed in movie, shorter than 3 movie frames)
    [backup_movie_LED,throwOutMovie]=throwOutOnStretches(backup_movie_LED,1:length(backup_movie_LED),3);
    backup_movie_LED=resample(backup_movie_LED,floor(mod(size_of_arduino/size_of_movie,1)*100)+floor((guess_best_scale*100)/100)*100,100);
    throwOutMovie=interp(throwOutMovie,movie_dec);
    movie_LED_for_finalalignment(throwOutMovie>0.5)=0;
    figz(4)=figure(); 
    plot(allEvents_movie_LED,'Color','b'); 
    hold on; 
    plot(movie_LED_for_finalalignment,'Color','r');
    legend({'before removal','after removal'});
    title('Remove short, skipped frame LED distractors before alignment');
else
    figz(4)=figure;
    backup_movie_LED=movie_LED;
    forSecondaryAlignment=backup_movie_LED;
    backup_arduino_LED=arduino_LED;
end
for j=1:length(tryscales)
    if mod(j,10)==0
        disp('Processing ...');
%         disp(j);
    end
    currscale=tryscales(j);
    movie_LED=resample(backup_movie_LED,floor(currscale*(1/tryinc)),floor((1/tryinc)));
    for i=1:length(trydelays)
        currdelay=trydelays(i);
        if mod(i,500)==0
%             disp(i);
        end 
        temp_movie=[nan(1,currdelay) movie_LED];
        temp_arduino=[arduino_LED nan(1,length(temp_movie)-length(arduino_LED))];
        if length(temp_arduino)>length(temp_movie)
            temp_movie=[temp_movie nan(1,length(temp_arduino)-length(temp_movie))];
            sumdiffs(j,i)=nansum(abs(temp_movie-temp_arduino));
        elseif length(temp_movie)>length(temp_arduino)
            % This should not happen
            % Arduino should include movie
            sumdiffs(j,i)=inf;
            error('Why is movie longer than Arduino?');
        else
            sumdiffs(j,i)=nansum(abs(temp_movie-temp_arduino)); 
        end
    end
end
sumdiffs(isinf(sumdiffs))=nan;
ma=max(sumdiffs(:));
sumdiffs(isnan(sumdiffs))=3*ma;
[minval,mi]=min(sumdiffs(:));
[mi_row,mi_col]=ind2sub(size(sumdiffs),mi);

figz(5)=figure(); 
imagesc(sumdiffs);
title('Finding best alignment');
xlabel('Trying different delays');
ylabel('Trying different scales');
disp('Best delay column');
disp(mi_col);
disp('Best scale row');
disp(mi_row);

end
