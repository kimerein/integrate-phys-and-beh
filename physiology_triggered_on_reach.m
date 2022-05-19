function [physiology_tbt,behavior_tbt,data]=physiology_triggered_on_reach(datadir_forBehavior,auxData)

% rig-specific settings
% distractor threshold
distract_thresh_photometry=0.5;
minTimeBetweenCues=2; % in seconds, include some slop in alignment
useCue='cueZone_onVoff';
chronuxPathWFindpeaks='C:\Users\sabatini\Documents\MATLAB\chronux_2_11\spectral_analysis\continuous';

rmpath(chronuxPathWFindpeaks);

% other settings
distract_thresh_movie=0.5; % from alignment

data.Fs=auxData.distractorData.ADFreq;
data.distractor=auxData.distractorData.Values;
data.phys_timepoints=0:1/data.Fs:(length(data.distractor)-1)*(1/data.Fs);
data.phys_inds=1:length(data.distractor);
data.cue=auxData.cueData.Values;
data.opto=auxData.optoData.Values;
% If also photometry_times as a way to keep track of indexing into
% photometry data, rename photometry_times as green_time
if isfield(auxData,'photometry_times')
    data.green_time=auxData.photometry_times;
end
% make all 1D fields row vectors
f=fieldnames(data);
for i=1:length(f)
    if size(data.(f{i}),2)==1 && size(data.(f{i}),1)>1
        % column vector
        temp=data.(f{i});
        data.(f{i})=temp';
    end
end

% load already processed behavior data
if iscell(datadir_forBehavior)
    a=cell(1,length(datadir_forBehavior));
    t=cell(1,length(datadir_forBehavior));
    for i=1:length(datadir_forBehavior)
        currdir=datadir_forBehavior{i};
        al=load([currdir '\final_aligned_data.mat']);
        a{i}=al.alignment;
        tbtl=load([currdir '\tbt.mat']);
        t{i}=tbtl.tbt;
    end
    
    if length(a)>1
        % concat alignments
        totalalignment=concatAlignmnents(a{1},a{2});
        for i=3:length(datadir_forBehavior)
            totalalignment=concatAlignmnents(totalalignment,a{i});
        end
    
        alltbt=concatTbt(t{1},t{2});
        for i=3:length(datadir_forBehavior)
            alltbt=concatTbt(alltbt,t{i});
        end
    else
        totalalignment=a{1};
        totalalignment.from_first_video=ones(size(totalalignment.cue));
        alltbt=t{1};
        alltbt.from_first_video=ones(size(alltbt.cue));
    end
else
    a=load([datadir_forBehavior '\final_aligned_data.mat']);
    totalalignment=a.alignment;
end

[totalalignment,~,cueInds]=findMatchedCues(totalalignment,useCue);
disp(['found ' num2str(length(cueInds)) ' cues in totalalignment']);

% Need to resample movie and arduino data so that indices represent
% matching times
settings.photo_fs=data.Fs; % photometry data sampling rate in Hz
tempalign1=a{1};
mo=mode(diff(tempalign1.timesfromarduino));
if mo==0
    temp=diff(tempalign1.timesfromarduino);
    temp(temp==0)=nan;
    mo=mode(temp);
end
settings.movie_fs=1/(mo/1000); % movie data sampling rate in Hz
disp(['Check this. Movie data sampling rate is ' num2str(settings.movie_fs)]);
pause;
settings.scale_factor=floor(settings.photo_fs/settings.movie_fs);
settings.photo_dec=settings.scale_factor-15;
settings.movie_dec=1;
settings.minlagForInitialAlign=[];
settings.maxlagForInitialAlign=[]; % [] is don't want to constrain alignment
settings.try_delay1=0;
settings.try_delay2=0;

% will adjust time subtly (<1%) in the behavior data to match the photometry
% will not alter the photometry data itself
data.distractor_times=0:(1/data.Fs):(length(data.distractor)-1)*(1/data.Fs);
data.cue_times=0:(1/data.Fs):(length(data.distractor)-1)*(1/data.Fs);
% because there is a break between the videos, and user can throw out first
% frames in each video, align each video separately
[discardedPhotoFrames_time,frontShift_time,scaleBy,movie_LED,movie_times,scaleMovieTimes,addToMovieTimes,padPhotoTimesAtFront]=alignDistractors(totalalignment.movie_distractor(totalalignment.from_first_video==1),data.distractor,distract_thresh_movie,distract_thresh_photometry,totalalignment.timesfromarduino(totalalignment.from_first_video==1),data.distractor_times,settings,[],[]);
alltbt1=scaleMovTimes(alltbt,scaleMovieTimes,addToMovieTimes,alltbt.from_first_video(:,1));
[tbt_data_vid1,shifted_data_vid1,alltbt1]=shiftPhotometryToBehavior(data,discardedPhotoFrames_time,frontShift_time,movie_LED,movie_times,totalalignment.(useCue),totalalignment,alltbt1,minTimeBetweenCues,totalalignment.from_first_video==1,alltbt.from_first_video(:,1),padPhotoTimesAtFront);
fromvid=alltbt.from_first_video(:,1);
if isfield(alltbt1,'threwOutCuesAtEnd')
    fromvid(end-alltbt1.threwOutCuesAtEnd+1:end)=0;
end
if isfield(alltbt1,'threwOutCuesAtFront')
    fromvid(find(fromvid==1,1,'first'):find(fromvid==1,1,'first')+(alltbt1.threwOutCuesAtFront-1))=0;
end
alltbt1=selectRows(alltbt,fromvid);

if length(a)>1
%     isInBackHalf=true;
%     fractionThru=0.85;
    isInBackHalf=[];
    fractionThru=[];
    [discardedPhotoFrames_time,frontShift_time,scaleBy,movie_LED,movie_times,scaleMovieTimes,addToMovieTimes,padPhotoTimesAtFront]=alignDistractors(totalalignment.movie_distractor(totalalignment.from_second_video==1),data.distractor,distract_thresh_movie,distract_thresh_photometry,totalalignment.timesfromarduino(totalalignment.from_second_video==1),data.distractor_times,settings,isInBackHalf,fractionThru);
    alltbt2=scaleMovTimes(alltbt,scaleMovieTimes,addToMovieTimes,alltbt.from_second_video(:,1));
    [tbt_data_vid2,shifted_data_vid2,alltbt2]=shiftPhotometryToBehavior(data,discardedPhotoFrames_time,frontShift_time,movie_LED,movie_times,totalalignment.(useCue),totalalignment,alltbt2,minTimeBetweenCues,totalalignment.from_second_video==1,alltbt.from_second_video(:,1),padPhotoTimesAtFront);
    fromvid=alltbt.from_second_video(:,1);
    if isfield(alltbt2,'threwOutCuesAtEnd')
        fromvid(end-alltbt2.threwOutCuesAtEnd+1:end)=0;
    end
    if isfield(alltbt2,'threwOutCuesAtFront')
        fromvid(find(fromvid==1,1,'first'):find(fromvid==1,1,'first')+(alltbt2.threwOutCuesAtFront-1))=0;
    end
    alltbt2=selectRows(alltbt,fromvid);
    
    if isfield(totalalignment,'from_third_video')
%         isInBackHalf=true;
%         fractionThru=0.95;
        isInBackHalf=[];
        fractionThru=[];
        [discardedPhotoFrames_time,frontShift_time,scaleBy,movie_LED,movie_times,scaleMovieTimes,addToMovieTimes,padPhotoTimesAtFront]=alignDistractors(totalalignment.movie_distractor(totalalignment.from_third_video==1),data.distractor,distract_thresh_movie,distract_thresh_photometry,totalalignment.timesfromarduino(totalalignment.from_third_video==1),data.distractor_times,settings,isInBackHalf,fractionThru);
        alltbt3=scaleMovTimes(alltbt,scaleMovieTimes,addToMovieTimes,alltbt.from_third_video(:,1));
        [tbt_data_vid3,shifted_data_vid3,alltbt3]=shiftPhotometryToBehavior(data,discardedPhotoFrames_time,frontShift_time,movie_LED,movie_times,totalalignment.(useCue),totalalignment,alltbt3,minTimeBetweenCues,totalalignment.from_third_video==1,alltbt.from_third_video(:,1),padPhotoTimesAtFront);
        fromvid=alltbt.from_third_video(:,1);
        if isfield(alltbt3,'threwOutCuesAtEnd')
            fromvid(end-alltbt3.threwOutCuesAtEnd+1:end)=0;
        end
        if isfield(alltbt3,'threwOutCuesAtFront')
            fromvid(find(fromvid==1,1,'first'):find(fromvid==1,1,'first')+(alltbt3.threwOutCuesAtFront-1))=0;
        end
        alltbt3=selectRows(alltbt,fromvid);
    end
end

% ds=100; % use this downsample factor for the fields in downsampfields
% downsampfields={'green_mod','red_mod','opto','cue','distractor','cue_times'};
% if ds~=1
%     tbt_data_vid1=downSampleData(tbt_data_vid1,ds,downsampfields);
%     if length(a)>1
%         tbt_data_vid2=downSampleData(tbt_data_vid2,ds,downsampfields);
%         if isfield(totalalignment,'from_third_video')
%             tbt_data_vid3=downSampleData(tbt_data_vid3,ds,downsampfields);
%         end
%     end
% end

% concatenate tbts from two videos, should match cues in behavior tbt
if length(a)>1
    [photometry_tbt]=concatTbt(tbt_data_vid1,tbt_data_vid2);
    [behavior_tbt]=concatTbt(alltbt1,alltbt2);
    if isfield(totalalignment,'from_third_video')
        [photometry_tbt]=concatTbt(photometry_tbt,tbt_data_vid3);
        [behavior_tbt]=concatTbt(behavior_tbt,alltbt3);
    end
else
    photometry_tbt=tbt_data_vid1;
    behavior_tbt=alltbt1;
end

% add times wrt trial start for photometry
if isfield(data,'green_time')
    photometry_tbt.trialDurations=nanmax(photometry_tbt.green_time,[],2)-photometry_tbt.green_time(:,1);
    photometry_tbt.times_wrt_trial_start=photometry_tbt.green_time-repmat(photometry_tbt.green_time(:,1),1,size(photometry_tbt.green_time,2));
    photometry_tbt.times_wrt_trial_start(photometry_tbt.times_wrt_trial_start<0)=nan;
end
photometry_tbt.cuetimes_wrt_trial_start=photometry_tbt.cue_times-repmat(photometry_tbt.cue_times(:,1),1,size(photometry_tbt.cue_times,2));
photometry_tbt.cuetimes_wrt_trial_start(photometry_tbt.cuetimes_wrt_trial_start<0)=nan;

% Actually this is physiology
physiology_tbt=photometry_tbt;

end

function data=downSampleData(data,ds,downsampfields)

for i=1:length(downsampfields)
    data.(downsampfields{i})=downSampMatrix(data.(downsampfields{i}),ds);
end

end

function alltbt=selectRows(alltbt,takeTheseRows)

f=fieldnames(alltbt);
for i=1:length(f)
    if strcmp(f{i},'threwOutCuesAtEnd') || strcmp(f{i},'threwOutCuesAtFront')
        continue
    end
    temp=alltbt.(f{i});
    temp=temp(takeTheseRows==1,:);
    alltbt.(f{i})=temp;
end

end

function alltbt=scaleMovTimes(alltbt,scaleMovieTimes,addToMovieTimes,takeTheseRows)

fields_like_times={'timesfromarduino','times','times_wrt_trial_start'};
for i=1:length(fields_like_times)
    temp=alltbt.(fields_like_times{i});
    temp(takeTheseRows==1,:)=(temp(takeTheseRows==1,:)*scaleMovieTimes)+addToMovieTimes;
    alltbt.(fields_like_times{i})=temp;
end

end

function [alltbt]=concatTbt(tbt1,tbt2)

f=fieldnames(tbt1);
for i=1:length(f)
    if ~isfield(tbt2,f{i})
        continue
    end
    if size(tbt1.(f{i}),2)>size(tbt2.(f{i}),2)
        tbt2.(f{i})=[tbt2.(f{i}) nan(size(tbt2.(f{i}),1),size(tbt1.(f{i}),2)-size(tbt2.(f{i}),2))]; % pad w nans at the end
    elseif size(tbt1.(f{i}),2)<size(tbt2.(f{i}),2)
        tbt1.(f{i})=[tbt1.(f{i}) nan(size(tbt1.(f{i}),1),size(tbt2.(f{i}),2)-size(tbt1.(f{i}),2))];
    end
    alltbt.(f{i})=[tbt1.(f{i}); tbt2.(f{i})];
end

if ~isfield(tbt1,'from_second_video') 
    alltbt.from_first_video=zeros(size(alltbt.(f{1})));
    alltbt.from_first_video(1:size(tbt1.(f{1}),1),:)=1;
    alltbt.from_second_video=zeros(size(alltbt.(f{1})));
    alltbt.from_second_video(size(tbt1.(f{1}),1)+1:end,:)=1;
elseif ~isfield(tbt1,'from_third_video') && isfield(tbt1,'from_second_video')
    alltbt.from_first_video=zeros(size(alltbt.(f{1})));
    alltbt.from_first_video(1:nansum(tbt1.from_first_video(:,1)),:)=tbt1.from_first_video(1:nansum(tbt1.from_first_video(:,1)),:);
    alltbt.from_second_video=zeros(size(alltbt.(f{1})));
    alltbt.from_second_video(nansum(tbt1.from_first_video(:,1))+1:nansum(tbt1.from_first_video(:,1))+nansum(tbt1.from_second_video(:,1)),:)=tbt1.from_second_video(nansum(tbt1.from_first_video(:,1))+1:nansum(tbt1.from_first_video(:,1))+nansum(tbt1.from_second_video(:,1)),:);
    alltbt.from_third_video=zeros(size(alltbt.(f{1})));
    alltbt.from_third_video(size(tbt1.(f{1}),1)+1:end,:)=1;
else
end

end

function outalignment=concatAlignmnents(alignment1,alignment2)

f=fieldnames(alignment1);
for i=1:length(f)
    if strcmp(f{i},'timesfromarduino')
        timestep=mode(alignment1.timesfromarduino(2:end)-alignment1.timesfromarduino(1:end-1));
        nm=nanmin(alignment2.timesfromarduino);
        outalignment.(f{i})=[alignment1.(f{i}) alignment1.timesfromarduino(end)+timestep-nm+alignment2.(f{i})];
    else
        if ~isfield(alignment2,f{i})
            continue
        end
        outalignment.(f{i})=[alignment1.(f{i}) alignment2.(f{i})];
    end
end

if ~isfield(alignment1,'from_first_video')
    outalignment.from_first_video=zeros(size(outalignment.(f{1})));
    outalignment.from_first_video(1:length(alignment1.(f{1})))=1;
    outalignment.from_second_video=zeros(size(outalignment.(f{1})));
    outalignment.from_second_video(length(alignment1.(f{1}))+1:end)=1;
elseif ~isfield(alignment1,'from_third_video')
    outalignment.from_first_video=zeros(size(outalignment.(f{1})));
    outalignment.from_first_video(1:nansum(alignment1.from_first_video))=alignment1.from_first_video(1:nansum(alignment1.from_first_video));
    outalignment.from_second_video=zeros(size(outalignment.(f{1})));
    outalignment.from_second_video(nansum(alignment1.from_first_video)+1:nansum(alignment1.from_first_video)+nansum(alignment1.from_second_video))=alignment1.from_second_video(nansum(alignment1.from_first_video)+1:nansum(alignment1.from_first_video)+nansum(alignment1.from_second_video));
    outalignment.from_third_video=zeros(size(outalignment.(f{1})));
    outalignment.from_third_video(length(alignment1.(f{1}))+1:end)=1;
else
end

end

function [tbt_data,shifted_data,alltbt]=shiftPhotometryToBehavior(data,discardedPhotoFrames_time,frontShift_time,movie_LED,movie_times,movie_cue,totalalignment,alltbt,minTimeBetweenCues,fromCurrVid,fromCurrVid_tbt,padPhotoTimesAtFront)

% note that only makes trial-by-trial (tbt) for the current video

% will not resample photometry or physiology, might introduce artifacts
% only shift
% then make trial-by-trial and align at cues

fields_like_distractor={'cue','opto','cue_times','phys_timepoints','phys_inds'};
fields_like_photometry={'green_time'};
if ~isfield(data,'green_time')
    fields_like_photometry={};
end

shifted_data=data;
% start with distractor
% discard time at front
padPhotoFrames_inds=round(padPhotoTimesAtFront./mode(data.distractor_times(2:end)-data.distractor_times(1:end-1)));
discardedPhotoFrames_inds=round(discardedPhotoFrames_time./mode(data.distractor_times(2:end)-data.distractor_times(1:end-1)));
shiftByInds=round(frontShift_time./mode(data.distractor_times(2:end)-data.distractor_times(1:end-1)));
shifted_data.distractor=[nan(1,padPhotoFrames_inds) shifted_data.distractor]; 
shifted_data.distractor=shifted_data.distractor(discardedPhotoFrames_inds+1:end); 
shifted_data.distractor=shifted_data.distractor(shiftByInds+1:end);
shifted_data.distractor_times=[nan(1,padPhotoFrames_inds) shifted_data.distractor_times]; 
shifted_data.distractor_times=shifted_data.distractor_times(discardedPhotoFrames_inds+1:end);
shifted_data.distractor_times=shifted_data.distractor_times(shiftByInds+1:end);

figure(); 
plot(shifted_data.distractor_times,shifted_data.distractor./nanmax(shifted_data.distractor),'Color','b'); 
hold on;
plot(movie_times,movie_LED,'Color','r');
title('Shifted distractors');

% then discard same times from fields like distractor
for i=1:length(fields_like_distractor)
    % discard time at front
    temp=shifted_data.(fields_like_distractor{i});
    temp=[nan(1,padPhotoFrames_inds) temp]; 
    temp=temp(discardedPhotoFrames_inds+1:end);
    shifted_data.(fields_like_distractor{i})=temp(shiftByInds+1:end);
end
disp(['min cue times in shifted data is ' num2str(nanmin(shifted_data.cue_times))]);

% then discard same times from fields like photometry
if ~isempty(fields_like_photometry)
    if isfield(data,'green_time')
        padPhotoFrames_inds_photo=floor(padPhotoTimesAtFront./mode(data.green_time(2:end)-data.green_time(1:end-1)));
        discardedPhotoFrames_inds_photo=floor(discardedPhotoFrames_time./mode(data.green_time(2:end)-data.green_time(1:end-1)));
        shiftByInds_photo=floor(frontShift_time./mode(data.green_time(2:end)-data.green_time(1:end-1)));
        for i=1:length(fields_like_photometry)
            % discard time at front
            temp=shifted_data.(fields_like_photometry{i});
            temp=[nan(1,padPhotoFrames_inds_photo) temp];
            temp=temp(discardedPhotoFrames_inds_photo+1:end);
            shifted_data.(fields_like_photometry{i})=temp(shiftByInds_photo+1:end);
        end
        disp(['min photometry times in shifted data is ' num2str(nanmin(shifted_data.green_time))]);
    end
end

movie_cue=movie_cue(fromCurrVid);

figure(); 
plot(shifted_data.distractor_times,shifted_data.cue./nanmax(shifted_data.cue),'Color','b'); 
hold on;
plot(movie_times,movie_cue,'Color','r');
title('All shifted cues');
xlabel('Time (sec)');

figure(); 
plot(shifted_data.distractor_times,shifted_data.cue./nanmax(shifted_data.cue),'Color','b'); 
hold on;
plot(movie_times,movie_cue,'Color','r');
title('Shifted cues - current video');
xlabel('Time (sec)');
anchor=input('Choose index for alignment anchor. ');
minSpacingBetweenCues_sig1=floor(minTimeBetweenCues/mode(shifted_data.distractor_times(2:end)-shifted_data.distractor_times(1:end-1)));
temp=movie_times(2:end)-movie_times(1:end-1);
temp=temp(temp~=0);
minSpacingBetweenCues_sig2=floor(minTimeBetweenCues/mode(temp));
[~,anchor_index_signal1]=nanmin(abs(shifted_data.distractor_times-anchor));
[~,anchor_index_signal2]=nanmin(abs(movie_times-anchor));
zo=input('Do you want to zero out any of the blue cues? If yes, enter time start and time end as [start end] to zero out. ');
if ~isempty(zo)
    shifted_data.cue=zeroOutPartOfOne(shifted_data.cue,shifted_data.distractor_times,zo); disp('zeroing out');
end
[mapping_signal1,mapping_signal2]=countCues(shifted_data.cue./nanmax(shifted_data.cue),shifted_data.distractor_times,movie_cue,movie_times,0.5,anchor_index_signal1,anchor_index_signal2,minSpacingBetweenCues_sig1,minSpacingBetweenCues_sig2);
% signal2 is movie, signal1 is data.distractor
% sort cues from first to last
[~,si]=sort(mapping_signal2);
mapping_signal2=mapping_signal2(si);
mapping_signal1=mapping_signal1(si);
if isfield(data,'green_time')
    mapping_signal3=getPhotoCueMapping(mapping_signal1,shifted_data.distractor_times,shifted_data.green_time);
    mapping_signal3=mapping_signal3(si);
else
    mapping_signal3=[];
end

fields_like_distractor{length(fields_like_distractor)+1}='distractor';
[tbt_data,alltbt]=makeTbtData(shifted_data,alltbt,fromCurrVid_tbt,mapping_signal1,mapping_signal2,mapping_signal3,'cueZone_onVoff',fields_like_distractor,fields_like_photometry);

end

function data=zeroOutPartOfOne(data,datatimes,zeroOutTime)

data(datatimes>zeroOutTime(1) & datatimes<zeroOutTime(2))=0;

end

function [data,cue,cueInds,cueIndITIs]=findMatchedCues(data,nameOfCue)
% data is alignment

% Get cue/data type for triggering trial-by-trial data
cue=data.(nameOfCue);

settings=plotCueTriggered_settings();
% In case of issues with aliasing of instantaneous cue
maxITI=settings.maxITI; % in seconds, maximal ITI
minITI=settings.minITI; % in seconds, minimal ITI

% Get time delay
timeIncs=diff(data.timesfromarduino(data.timesfromarduino~=0));
mo=mode(timeIncs);
bettermode=mo;
if mo==0
    timeIncs(timeIncs==mo)=nan;
    bettermode=mode(timeIncs); % in ms
end
bettermode=bettermode/1000; % in seconds

% Fix aliasing issues with resampled data
% if strcmp(nameOfCue,'cueZone_onVoff')
if strcmp(nameOfCue,'cue') | strcmp(nameOfCue,'cueZone_onVoff') | strcmp(nameOfCue,'falseCueOn') | strcmp(nameOfCue,'movie_distractor') | strcmp(nameOfCue,'optoOnly')
    [cue,cueInds,cueIndITIs]=fixAlias_forThreshCue(cue,maxITI,minITI,bettermode);
else
    [cue,cueInds,cueIndITIs]=fixAliasing(cue,maxITI,minITI,bettermode);
end
data.(nameOfCue)=cue;
[data.pelletPresented,presentedInds]=fixAliasing(data.pelletPresented,maxITI,minITI,bettermode);

end

function new_mapping=getPhotoCueMapping(currcuemapping,currtimes,newtimes)

% gets cue inds into down-sampled photometry data
new_mapping=nan(1,length(currcuemapping));
for i=1:length(currcuemapping)
    [~,new_mapping(i)]=nanmin(abs(newtimes-currtimes(currcuemapping(i))));
end

end

function [tbt_data,alltbt]=makeTbtData(data,alltbt,fromCurrVid,cue_mapping_data_distract,cue_mapping_movie_distract,cue_mapping_data_photo,useCue,fields_like_distractor,fields_like_photometry)

temp=alltbt.(useCue);
disp([num2str(size(temp(fromCurrVid==1,:),1)) ' cues in behavior tbt for this vid']);
disp([num2str(length(cue_mapping_data_distract)) ' cues in data tbt']);
alltbt.threwOutCuesAtFront=0;
alltbt.threwOutCuesAtEnd=0;
if size(temp(fromCurrVid==1,:),1)~=length(cue_mapping_data_distract)
    n=input('Number of cues in data does not match number of cues in movie. Could be a problem. To throw out n cues at end of movie (blue), enter positive number (negative number for front of movie). Or [neg (front) pos (end)]. Or quit. n: ');
    if length(n)==1
        if n>0
            n=[0 n];
        elseif n<0
            n=[n 0];
        else
            n=[0 0];
        end
    end
    if n(2)>0
        % throw out n cues at end of video
        f=fieldnames(alltbt);
        for i=1:length(f)
            if strcmp(f{i},'threwOutCuesAtEnd') || strcmp(f{i},'threwOutCuesAtFront')
                continue
            end
            temp=alltbt.(f{i});
            temp=temp(1:end-n(2),:);
            alltbt.(f{i})=temp;
        end
        alltbt.threwOutCuesAtEnd=n(2);
    end
    if n(1)<0
        n(1)=abs(n(1));
        % throw out n cues at beginning of video
        f=fieldnames(alltbt);
        for i=1:length(f)
            if strcmp(f{i},'threwOutCuesAtEnd') || strcmp(f{i},'threwOutCuesAtFront')
                continue
            end
            temp=alltbt.(f{i});
            temp=temp(n(1)+1:end,:);
            alltbt.(f{i})=temp;
        end
        alltbt.threwOutCuesAtFront=n(1);
    end
end

% get size of baseline in behavior tbt
cueInd=find(nanmean(alltbt.(useCue),1)>0.5,1,'first');
secBeforeCue=nanmean(alltbt.times_wrt_trial_start(:,cueInd));
distractorIndsBeforeCue=floor(secBeforeCue/mode(data.distractor_times(2:end)-data.distractor_times(1:end-1)));
maxDistractorIndsAfterCue=nanmax(diff(cue_mapping_data_distract));
maxTimeAfterCue=maxDistractorIndsAfterCue*mode(data.distractor_times(2:end)-data.distractor_times(1:end-1));
if isfield(data,'green_time')
    photoIndsBeforeCue=floor(secBeforeCue/mode(data.green_time(2:end)-data.green_time(1:end-1)));
    maxPhotoIndsAfterCue=floor(maxTimeAfterCue/mode(data.green_time(2:end)-data.green_time(1:end-1)));
end

% line up each cue in data to get trial-by-trial (tbt)
% make each row include the pre-cue baseline and all points up to but not
% including the next cue (same as in plotCueTriggeredBehavior.m)

for i=1:length(fields_like_distractor)
    % align distractor-like fields
    tempdata=data.(fields_like_distractor{i});
    temptbt=nan(length(cue_mapping_movie_distract),maxDistractorIndsAfterCue+distractorIndsBeforeCue);
    for j=1:length(cue_mapping_movie_distract)
        % for each movie cue, make a data cue
        if cue_mapping_data_distract(j)-distractorIndsBeforeCue<1
            firstInd=1;
        else
            firstInd=cue_mapping_data_distract(j)-distractorIndsBeforeCue;
        end
        if j==length(cue_mapping_movie_distract)
            if cue_mapping_data_distract(j)+maxDistractorIndsAfterCue+distractorIndsBeforeCue-1>length(tempdata)
                secondInd=length(tempdata);
            else
                secondInd=cue_mapping_data_distract(j)+maxDistractorIndsAfterCue+distractorIndsBeforeCue-1;
            end
        else
            secondInd=cue_mapping_data_distract(j+1)-1;
        end
        if cue_mapping_data_distract(j)-distractorIndsBeforeCue<1
            tempie=[nan(1,abs(cue_mapping_data_distract(j)-distractorIndsBeforeCue)) tempdata(firstInd:secondInd)];
        else
            tempie=tempdata(firstInd:secondInd);
        end
        temptbt(j,1:length(tempie))=tempie;
    end
    tbt_data.(fields_like_distractor{i})=temptbt;
end

for i=1:length(fields_like_photometry)
    % align photometry-like fields
    tempdata=data.(fields_like_photometry{i});
    temp=nan(length(cue_mapping_movie_distract),photoIndsBeforeCue+maxPhotoIndsAfterCue);
    for j=1:length(cue_mapping_movie_distract)
        % for each movie cue, make a data cue
        if cue_mapping_data_photo(j)-photoIndsBeforeCue<1
            firstInd=1;
        else
            firstInd=cue_mapping_data_photo(j)-photoIndsBeforeCue;
        end
        if j==length(cue_mapping_data_photo)
            if cue_mapping_data_photo(j)+maxPhotoIndsAfterCue+photoIndsBeforeCue-1>length(tempdata)
                secondInd=length(tempdata);
            else
                secondInd=cue_mapping_data_photo(j)+maxPhotoIndsAfterCue+photoIndsBeforeCue-1;
            end
        else
            secondInd=cue_mapping_data_photo(j+1)-1;
        end
        if cue_mapping_data_photo(j)-photoIndsBeforeCue<1
            tempie=[nan(1,abs(cue_mapping_data_photo(j)-photoIndsBeforeCue)) tempdata(firstInd:secondInd)];
        else
            tempie=tempdata(firstInd:secondInd);
        end
        temp(j,1:length(tempie))=tempie;
    end
    tbt_data.(fields_like_photometry{i})=temp;
end

end

function [mapping_signal1,mapping_signal2]=countCues(signal1,sig1_x,signal2,sig2_x,thresh,anchor1,anchor2,minSpacingBetweenCues_sig1,minSpacingBetweenCues_sig2)

countedcue_signal1=signal1;
countedcue_signal2=signal2;

% find signal1>thresh closest to index anchor
f1=find(signal1>thresh);
[~,mi]=nanmin(abs(f1-anchor1));
closest1=f1(mi);

% find signal2>thresh closest to anchor2
f2=find(signal2>thresh);
[~,mi]=nanmin(abs(f2-anchor2));
closest2=f2(mi);

mapping_signal1(1)=closest1;
mapping_signal2(1)=closest2;

% moving to the right, count cues after this anchor
backup_countedcue_signal1=countedcue_signal1;
backup_countedcue_signal2=countedcue_signal2;
countedcue_signal1(1:closest1)=0;
countedcue_signal2(1:closest2)=0;
if closest1+minSpacingBetweenCues_sig1>length(countedcue_signal1)
    takeInds=closest1:length(countedcue_signal1);
else
    takeInds=closest1:closest1+minSpacingBetweenCues_sig1;
end
countedcue_signal1(takeInds)=0;
if closest2+minSpacingBetweenCues_sig2>length(countedcue_signal2)
    takeInds=closest2:length(countedcue_signal2);
else
    takeInds=closest2:closest2+minSpacingBetweenCues_sig2;
end
countedcue_signal2(takeInds)=0;
remainingInds=length(signal1(closest1:end));
j=2;
for i=1:remainingInds
    % if there are no more points > thresh, break
    if nansum(countedcue_signal1(closest1:end)>thresh)==0 || nansum(countedcue_signal2(closest2:end)>thresh)==0
        break
    end
    f1=find(countedcue_signal1>thresh,1,'first');
    mapping_signal1(j)=f1;
    if f1+minSpacingBetweenCues_sig1>length(countedcue_signal1)
        takeInds=f1:length(countedcue_signal1);
    else
        takeInds=f1:f1+minSpacingBetweenCues_sig1;
    end
    countedcue_signal1(takeInds)=0;
    f2=find(countedcue_signal2>thresh,1,'first');
    mapping_signal2(j)=f2;
    if f2+minSpacingBetweenCues_sig2>length(countedcue_signal2)
        takeInds=f2:length(countedcue_signal2);
    else
        takeInds=f2:f2+minSpacingBetweenCues_sig2;
    end
    countedcue_signal2(takeInds)=0;
    j=j+1;
end

% moving to the left, count cues before this anchor
countedcue_signal1(1:closest1)=backup_countedcue_signal1(1:closest1);
countedcue_signal1(closest1:end)=0;
countedcue_signal2(1:closest2)=backup_countedcue_signal2(1:closest2);
countedcue_signal2(closest2:end)=0;
if closest1-minSpacingBetweenCues_sig1<1
    takeInds=1:closest1;
else
    takeInds=closest1-minSpacingBetweenCues_sig1:closest1;
end
countedcue_signal1(takeInds)=0;
if closest2-minSpacingBetweenCues_sig2<1
    takeInds=1:closest2;
else
    takeInds=closest2-minSpacingBetweenCues_sig2:closest2;
end
countedcue_signal2(takeInds)=0;
remainingInds=length(signal1(1:closest1));
for i=1:remainingInds
    % if there are no more points > thresh, break
    if nansum(countedcue_signal1(1:closest1)>thresh)==0 || nansum(countedcue_signal2(1:closest2)>thresh)==0
        break
    end
    f1=find(countedcue_signal1>thresh,1,'last');
    mapping_signal1(j)=f1;
    if f1-minSpacingBetweenCues_sig1<1
        takeInds=1:f1;
    else
        takeInds=f1-minSpacingBetweenCues_sig1:f1;
    end
    countedcue_signal1(takeInds)=0;
    f2=find(countedcue_signal2>thresh,1,'last');
    mapping_signal2(j)=f2;
    if f2-minSpacingBetweenCues_sig2<1
        takeInds=1:f2;
    else
        takeInds=f2-minSpacingBetweenCues_sig2:f2;
    end
    countedcue_signal2(takeInds)=0;
    j=j+1;
end

% Plot mapping of cues
figure();
plot(sig1_x,signal1,'Color','k');
hold on;
plot(sig2_x,signal2,'Color','b');
for i=1:length(mapping_signal1)
    c=rand(1,3);
    scatter(sig1_x(mapping_signal1(i)),nanmax(signal1),10,c,'filled');
    scatter(sig2_x(mapping_signal2(i)),nanmax(signal2),20,c);
end

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

temp1=arduino_LED_ITIs./max(arduino_LED_ITIs);
temp2=movie_LED_ITIs./max(movie_LED_ITIs);

if ~isempty(settings.minlagForInitialAlign)
    temp2=[nan(1,settings.minlagForInitialAlign) temp2];
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
        X=[zeros(1,D) arduino_LED_ITIs./max(arduino_LED_ITIs)];
        Y=movie_LED_ITIs./max(movie_LED_ITIs);
    elseif D<0
        Y=[zeros(1,-D) movie_LED_ITIs./max(movie_LED_ITIs)];
        X=arduino_LED_ITIs./max(arduino_LED_ITIs);
    end
end

if D>0
    % photometry starts after movie
    % pad photometry at front
    figure();
    plot(X,'Color','b');
    hold on;
    plot(Y,'Color','r');
    anchor=input('Choose index for alignment anchor. ');
    alignAtPeak_arduino=locs_arduino(anchor);
    alignAtPeak_movie=locs(anchor+D);
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
    alignAtPeak_movie=locs(anchor+D);
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

function [cue,cueInds,cueIndITIs]=fixAlias_forThreshCue(cue,maxITI,minITI,bettermode)

settings=plotCueTriggered_settings();
peakHeight=0.5;

[pks,locs]=findpeaks(cue);
cueInds=locs(pks>peakHeight);
% cueInds=[1 cueInds length(cue)]; % in case aliasing problem is at edges
cueIndITIs=diff(cueInds);
checkTheseIntervals=find(cueIndITIs*bettermode>(maxITI*1.5));
for i=1:length(checkTheseIntervals)
    indsIntoCue=cueInds(checkTheseIntervals(i))+floor((maxITI/2)./bettermode):cueInds(checkTheseIntervals(i)+1)-floor((maxITI/2)./bettermode);
    if any(cue(indsIntoCue)>0.001)
        [~,ma]=max(cue(indsIntoCue));
        cue(indsIntoCue(ma))=max(cue);
    end
end

% [pks,locs]=findpeaks(cue);
[pks,locs]=findpeaks(cue,'MinPeakDistance',floor((minITI*0.75)/bettermode),'MinPeakProminence',peakHeight);
cueInds=locs(pks>peakHeight);
cueIndITIs=diff(cueInds);
checkTheseIntervals=find(cueIndITIs*bettermode<(minITI*0.75));
if ~isempty(checkTheseIntervals)
    for i=1:length(checkTheseIntervals)
        cue(cueInds(checkTheseIntervals(i)))=0;
        cueInds(checkTheseIntervals(i))=nan;
    end
end
cueInds=cueInds(~isnan(cueInds));
cueIndITIs=diff(cueInds);

cue=cue./nanmax(cue);

end

function [cue,cueInds,cueIndITIs]=fixAliasing(cue,maxITI,minITI,bettermode)

cue=nonparamZscore(cue); % non-parametric Z score

settings=plotCueTriggered_settings();
peakHeight=nanmean(cue)+settings.nStdDevs*nanstd(cue);
relativePeakHeight=settings.nStdDevs*nanstd(cue);

[pks,locs]=findpeaks(cue);
cueInds=locs(pks>peakHeight);
% cueInds=[1 cueInds length(cue)]; % in case aliasing problem is at edges
cueIndITIs=diff(cueInds);
checkTheseIntervals=find(cueIndITIs*bettermode>(maxITI*1.5));
for i=1:length(checkTheseIntervals)
    indsIntoCue=cueInds(checkTheseIntervals(i))+floor((maxITI/2)./bettermode):cueInds(checkTheseIntervals(i)+1)-floor((maxITI/2)./bettermode);
    if any(cue(indsIntoCue)>0.001)
        [~,ma]=max(cue(indsIntoCue));
        cue(indsIntoCue(ma))=max(cue);
    end
end

% [pks,locs]=findpeaks(cue);
[pks,locs]=findpeaks(cue,'MinPeakDistance',floor((minITI*0.75)/bettermode),'MinPeakProminence',relativePeakHeight);
cueInds=locs(pks>peakHeight);
cueIndITIs=diff(cueInds);
checkTheseIntervals=find(cueIndITIs*bettermode<(minITI*0.75));
if ~isempty(checkTheseIntervals)
    for i=1:length(checkTheseIntervals)
        cue(cueInds(checkTheseIntervals(i)))=0;
        cueInds(checkTheseIntervals(i))=nan;
    end
end
cueInds=cueInds(~isnan(cueInds));
cueIndITIs=diff(cueInds);

cue=cue./nanmax(cue);

end
