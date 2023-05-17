function alignment=alignContinuous(alignment,data)

% Fs of photometry 2000
% Fs of movie 30 fps
% small videos each have ~7513 frames
nframessmallvid=7513;

% use distractors to align
photo_distract=data.distractor;
movie_distract=alignment.movie_distractor;

% resample
frames=alignment.movieframeinds;
movie_distract=resample(movie_distract,66,1); % from relative frame rates
frames=resample(frames,66,1);

% overall big alignment
[xa,ya,delay_big]=alignsignals(movie_distract,photo_distract);
movie_distract=xa;
photo_distract=ya;
% align onsets
movie_distract=double(movie_distract>0.5);
movie_distract=[diff(movie_distract)>0.5 nan];
photo_distract=double(photo_distract>0.5);
photo_distract=[diff(photo_distract)>0.5 nan];

% shift smaller segments bcz frames dropped at each video transition
maxvids=20;
currchunk_firstframe=1;
delay_small=nan(1,length(maxvids));
framechunks=nan(length(maxvids),2);
frames_notyet_aligned=frames;
for i=1:maxvids
    if currchunk_firstframe>length(movie_distract)
        endedat=i-1;
        break
    end
    f=find(frames_notyet_aligned>nframessmallvid*i,1,'first');
    framechunks(i,:)=[currchunk_firstframe f];
    if f>length(photo_distract)
        f=length(photo_distract);
        framechunks(i,:)=[currchunk_firstframe f];
        temp1=movie_distract(currchunk_firstframe:f); temp1(isnan(temp1))=0;
        temp2=photo_distract(currchunk_firstframe:f); temp2(isnan(temp2))=0;
        [xa,ya,delay_small(i)]=alignsignals(temp1,temp2);
        endedat=i;
        break
    elseif f>length(movie_distract)
        f=length(movie_distract);
        framechunks(i,:)=[currchunk_firstframe f];
        temp1=movie_distract(currchunk_firstframe:f); temp1(isnan(temp1))=0;
        temp2=photo_distract(currchunk_firstframe:f); temp2(isnan(temp2))=0;
        [xa,ya,delay_small(i)]=alignsignals(temp1,temp2);
        endedat=i;
        break
    end
    temp1=movie_distract(currchunk_firstframe:f); temp1(isnan(temp1))=0;
    temp2=photo_distract(currchunk_firstframe:f); temp2(isnan(temp2))=0;
    [xa,ya,delay_small(i)]=alignsignals(temp1,temp2);
    frames_notyet_aligned(currchunk_firstframe:f)=nan;
    currchunk_firstframe=f+1;
    endedat=i;
end
framechunks=framechunks(1:endedat,:);
delay_small=delay_small(1:endedat);

% movie_distract=alignByDelaySmall(movie_distract,delay_small,framechunks);
% plot output
% figure();
% plot(movie_distract,'Color','k'); hold on;
% plot(photo_distract,'Color','r');

% Align everything
f=fieldnames(alignment);
zeroOrOneFields={'all_reachBatch','reachBatch_drop_reachStarts','reachBatch_success_reachStarts','reachBatch_miss_reachStarts_pawOnWheel','reachBatch_success_reachStarts_pawOnWheel',...
    'reachBatch_miss_reachStarts','drop_reachStarts_pawOnWheel_backup','success_reachStarts_pawOnWheel_backup','drop_reachStarts_backup','success_reachStarts_backup',...
    'reachFidgetBegins','lickStarts','reach_ongoing','miss_reachStarts_pawOnWheel','drop_reachStarts_pawOnWheel','success_reachStarts_pawOnWheel','pawOnWheel',...
    'eating','pelletmissingreach_reachStarts','miss_reachStarts','drop_reachStarts','success_reachStarts','reachStarts_pelletPresent','reachEnds','reachStarts','isHold',...
    'isChewing'};
for i=1:length(f)
    temp=alignment.(f{i});
    temp=resample(temp,66,1);
    temp=alignByDelayBig(temp,delay_big);
    temp=alignByDelaySmall(temp,delay_small,framechunks);
    if ismember(f{i},zeroOrOneFields) % clean up 0 or 1 fields
        temp=double(temp>0.5);
    end
    alignment.(f{i})=temp;
end

% Plot distractors as sanity check
figure();
plot(alignment.movie_distractor,'Color','k'); hold on;
plot(data.distractor,'Color','r');

end

function movie_distract=alignByDelayBig(movie_distract,delay_big)

temp=movie_distract(1:end);
D=delay_big;
if D>0
    movie_distract=[nan(1,abs(D)-1) temp(1:end-abs(D)+1)];
elseif D<0
    movie_distract=[temp(abs(D):end) nan(1,abs(D)-1)];
else
    movie_distract=temp;
end

end

function movie_distract=alignByDelaySmall(movie_distract,delay_small,framechunks)

% align distractor according to delay_small
for i=1:length(delay_small)
    temp=movie_distract(framechunks(i,1):framechunks(i,2));
    D=delay_small(i);
    if D>0
        movie_distract(framechunks(i,1):framechunks(i,2))=[nan(1,abs(D)-1) temp(1:end-abs(D)+1)];
    elseif D<0
        movie_distract(framechunks(i,1):framechunks(i,2))=[temp(abs(D):end) nan(1,abs(D)-1)];
    else
        movie_distract(framechunks(i,1):framechunks(i,2))=temp;
    end
end

end
