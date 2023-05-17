function alignContinuous(alignment,data)

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
        [xa,ya,delay_small(i)]=alignsignals(movie_distract(currchunk_firstframe:f),photo_distract(currchunk_firstframe:f));
        frames_notyet_aligned(currchunk_firstframe:f)=nan;
        currchunk_firstframe=f+1;
        endedat=i;
        break
    elseif f>length(movie_distract)
        f=length(movie_distract);
        [xa,ya,delay_small(i)]=alignsignals(movie_distract(currchunk_firstframe:f),photo_distract(currchunk_firstframe:f));
        frames_notyet_aligned(currchunk_firstframe:f)=nan;
        currchunk_firstframe=f+1;
        endedat=i;
        break
    end
    [xa,ya,delay_small(i)]=alignsignals(movie_distract(currchunk_firstframe:f),photo_distract(currchunk_firstframe:f));
    frames_notyet_aligned(currchunk_firstframe:f)=nan;
    currchunk_firstframe=f+1;
    endedat=i;
end
framechunks=framechunks(1:endedat,:);
delay_small=delay_small(1:endedat);

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

% plot output
figure();
plot(movie_distract,'Color','k'); hold on;
plot(photo_distract,'Color','r');



end
