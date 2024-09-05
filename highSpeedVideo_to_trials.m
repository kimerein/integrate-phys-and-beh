function tbt=highSpeedVideo_to_trials(filename,cueVarName,distractorVarName,rawVarName,wheelVarName,settings)

% filename is file that contains output of getRigEvents.py
% i.e., the cue zone, wheel zone, distractor zone and raw differences
% between frames, in order to sort into trials
% 
% settings.fps = frames per second of video
% settings.distractorDuration = in seconds
% settings.cueDuration = in seconds
% settings.wheelSecBeforeCue = when wheel starts in seconds before cue
% e.g.,
%                    fps: 256
%     distractorDuration: 0.2500
%            cueDuration: 0.2500
%      wheelSecBeforeCue: 0.9600

% have user set these
useRawDiffsAsWheelTurn=true; % if true, will take the raw pixel change as indication of the wheel turning, else will use the wheel zone and wheel thresh
fillInWheelStarts=true;
cueThresh=100;
maxFramesToPlot=30000;
rawThresh=1; % for where video skips, this is after the wheel has been turning for a little bit, about 0.96 s in Kim's rig
distractorThresh=100;
% will likely want to downsample this data for alignment, because high
% speed video is huge
ds=10; % downsample by this amount
maxTrialLength=19; % no trial is longer than this in seconds
maxTrialInds=floor(maxTrialLength * settings.fps);

disp('Loading rig events file');

a=load(filename);
if isfield(a,'cueDiffs')
    % haven't yet found difference events
    foundDiffEvs=false;
    cueDiffs=a.(cueVarName);
    distractorDiffs=a.(distractorVarName);
    rawDiffs=a.(rawVarName);
    wheelDiffs=a.(wheelVarName);
else
    foundDiffEvs=true;
    cueDiffEvs_plus=double(a.([cueVarName '_plus']));
    cueDiffEvs_minus=double(a.([cueVarName '_minus']));
    rawDiffEvs_plus=double(a.([rawVarName '_plus']));
    rawDiffEvs_minus=double(a.([rawVarName '_minus']));
    distractorDiffEvs_plus=double(a.([distractorVarName '_plus']));
    distractorDiffEvs_minus=double(a.([distractorVarName '_minus']));
    wheelDiffEvs_plus=double(a.([wheelVarName '_plus']));
    wheelDiffEvs_minus=double(a.([wheelVarName '_minus']));
    howmanyframes=double(a.('howmanyframes'));
end    
    
cueInds=floor(settings.cueDuration * settings.fps);
distractorInds=floor(settings.distractorDuration * settings.fps);
wheelBeforeCueInds=floor(settings.wheelSecBeforeCue * settings.fps);
indsSlop=floor(0.05 * settings.fps);
% indsSlop=floor(0.1 * settings.fps);

% get possible cues first
% then look for wheel turns at a fixed time before cue onset
figure();
if foundDiffEvs==false
    framesmax=length(cueDiffs);
    if framesmax>maxFramesToPlot
        framesmax=maxFramesToPlot;
    end
    plot(cueDiffs(1:framesmax),'Color','b'); hold on;
    line([1 framesmax],[cueThresh cueThresh],'Color','b');
    line([1 framesmax],[-cueThresh -cueThresh],'Color','b');
else
    for i=1:length(cueDiffEvs_plus)
        scatter(cueDiffEvs_plus(i),cueThresh,[],'b'); hold on;
    end
    for i=1:length(cueDiffEvs_minus)
        scatter(cueDiffEvs_minus(i),-cueThresh,[],'b'); hold on;
    end
end
ylabel('Difference frame-to-frame in cue zone');
xlabel('Frames');
% find cues of specific duration
if foundDiffEvs==false
    [cueStarts,cueEnds]=findEventsOfSpecificDuration(cueDiffs,cueThresh,cueInds,indsSlop);
else
    pause;
    [cueStarts,cueEnds]=findDiffsOfSpecificDuration(cueDiffEvs_plus,cueDiffEvs_minus,cueInds,indsSlop);
end
% plot
for i=1:length(cueStarts)
    line([cueStarts(i) cueEnds(i)],[cueThresh-cueThresh*0.25 cueThresh-cueThresh*0.25],'Color','c','LineWidth',2);
end

% get wheel turns beginning at a specific time before cue
if foundDiffEvs==false
    wheelOnsets=find(rawDiffs>rawThresh);
    comboDiffs=zeros(size(cueDiffs));
    comboDiffs(cueStarts)=-1;
else
    if useRawDiffsAsWheelTurn==true
        wheelOnsets=rawDiffEvs_plus;
    else
        wheelOnsets=[wheelDiffEvs_plus wheelDiffEvs_minus];
    end
end
% only keep wheelOnsets closest to cue onset, because cue is clearer
% signal
wheelStarts=zeros(size(cueStarts));
wheelOnsets=[0 wheelOnsets];
for i=1:length(cueStarts)
    wheelStarts(i)=wheelOnsets(find(wheelOnsets<cueStarts(i)-10,1,'last'));
end
if fillInWheelStarts==true
    wheelStarts=cueStarts-wheelBeforeCueInds;
    wheelStarts(wheelStarts<0)=0;
end
if foundDiffEvs==false
    comboDiffs(wheelStarts)=1;
    [wheelStarts,wheelEnds]=findEventsOfSpecificDuration(comboDiffs,0.5,wheelBeforeCueInds,indsSlop);
    plot(rawDiffs(1:framesmax),'Color','g'); hold on;
    line([1 framesmax],[rawThresh rawThresh],'Color','g');
else
    for i=1:length(rawDiffEvs_plus)
        scatter(rawDiffEvs_plus(i),cueThresh*0.75,[],'g'); hold on;
    end
%     [wheelStarts,wheelEnds]=findDiffsOfSpecificDuration(wheelStarts,cueDiffEvs_plus,wheelBeforeCueInds,indsSlop);
    [wheelStarts,wheelEnds]=findDiffsOfSpecificDuration(wheelStarts,cueStarts,wheelBeforeCueInds,indsSlop);
end
for i=1:length(wheelStarts)
    line([wheelStarts(i) wheelEnds(i)],[cueThresh-cueThresh*0.3 cueThresh-cueThresh*0.3],'Color',[0 0.75 0],'LineWidth',2);
end

% get distractors
if foundDiffEvs==false
    plot(distractorDiffs(1:framesmax),'Color','r'); hold on;
    line([1 framesmax],[distractorThresh distractorThresh],'Color','r');
    line([1 framesmax],[-distractorThresh -distractorThresh],'Color','r');
    [distractorStarts,distractorEnds]=findEventsOfSpecificDuration(distractorDiffs,distractorThresh,distractorInds,indsSlop);
else
    for i=1:length(distractorDiffEvs_plus)
        scatter(distractorDiffEvs_plus(i),0.5*cueThresh,[],'r'); hold on;
    end
    for i=1:length(distractorDiffEvs_minus)
        scatter(distractorDiffEvs_minus(i),-0.5*cueThresh,[],'r'); hold on;
    end
    [distractorStarts,distractorEnds]=findDiffsOfSpecificDuration(distractorDiffEvs_plus,distractorDiffEvs_minus,distractorInds,indsSlop);
end
for i=1:length(distractorStarts)
    line([distractorStarts(i) distractorEnds(i)],[cueThresh-cueThresh*0.15 cueThresh-cueThresh*0.15],'Color','m','LineWidth',2);
end

if foundDiffEvs==false
    frameinds=1+[1:length(cueDiffs)];
else
    frameinds=1+[1:howmanyframes];
end
cueOn=zeros(size(frameinds));
for i=1:length(cueStarts)
    cueOn(cueStarts(i):cueEnds(i))=1;
end
distractorOn=zeros(size(frameinds));
for i=1:length(distractorStarts)
    distractorOn(distractorStarts(i):distractorEnds(i))=1;
end
wheelOn=zeros(size(frameinds));
for i=1:length(wheelStarts)
    wheelOn(wheelStarts(i):wheelEnds(i))=1;
end
% convert to tbt (i.e., "trial by trial")
tbt.highspeed_frameinds=nan(length(wheelStarts),maxTrialInds);
tbt.cueZone_onVoff=nan(length(wheelStarts),maxTrialInds);
tbt.movie_distractor=nan(length(wheelStarts),maxTrialInds);
tbt.wheelTurning=nan(length(wheelStarts),maxTrialInds);
% wheel is less reliable than cue zone -- trigger on cue, not wheel
for i=1:length(wheelStarts)
    if i+1>length(wheelStarts)
        if length(frameinds)>wheelStarts(i)+maxTrialInds
            takeInds=wheelStarts(i):wheelStarts(i)+maxTrialInds;
        else
            takeInds=wheelStarts(i):length(frameinds);
        end
    else
        takeInds=wheelStarts(i):wheelStarts(i+1)-1;
    end
    tbt.highspeed_frameinds(i,1:length(downSampAv(frameinds(takeInds),ds)))=downSampAv(frameinds(takeInds),ds);
    tbt.cueZone_onVoff(i,1:length(downSampAv(cueOn(takeInds),ds)))=downSampAv(cueOn(takeInds),ds);
    tbt.movie_distractor(i,1:length(downSampAv(distractorOn(takeInds),ds)))=downSampAv(distractorOn(takeInds),ds);
    tbt.wheelTurning(i,1:length(downSampAv(wheelOn(takeInds),ds)))=downSampAv(wheelOn(takeInds),ds);
end

end

function [cueStarts,cueEnds]=findDiffsOfSpecificDuration(cueDiffEvs_plus,cueDiffEvs_minus,cueInds,indsSlop)

fcue=cueDiffEvs_plus;
fcueoff=cueDiffEvs_minus;
if length(fcueoff)>length(fcue)
    fcueoff=fcueoff(1:length(fcue)); % truncate and drop the last events
elseif length(fcue)>length(fcueoff)
    fcue=fcue(1:length(fcueoff));
end
diff1=fcueoff-fcue; % if starts with cue on, will be positive
diff2=fcueoff(2:end)-fcue(1:end-1); % if starts with cue off
if mode(sign(diff1))<0
    % started with cue on to off
    % use diff2
    % find differences approx. matching cue duration
    currfcue=fcue(1:end-1);
    currfcueoff=fcueoff(2:end);
    cueStarts=currfcue(diff2>cueInds-indsSlop & diff2<cueInds+indsSlop);
    cueEnds=currfcueoff(diff2>cueInds-indsSlop & diff2<cueInds+indsSlop);
else
    % use diff1
    currfcue=fcue;
    currfcueoff=fcueoff;
    cueStarts=currfcue(diff1>cueInds-indsSlop & diff1<cueInds+indsSlop);
    cueEnds=currfcueoff(diff1>cueInds-indsSlop & diff1<cueInds+indsSlop);
end

end

function [cueStarts,cueEnds]=findEventsOfSpecificDuration(cueDiffs,cueThresh,cueInds,indsSlop)

fcue=find(cueDiffs>cueThresh);
fcueoff=find(cueDiffs<-cueThresh);
if length(fcueoff)>length(fcue)
    fcueoff=fcueoff(1:length(fcue)); % truncate and drop the last events
elseif length(fcue)>length(fcueoff)
    fcue=fcue(1:length(fcueoff));
end
diff1=fcueoff-fcue; % if starts with cue on, will be positive
diff2=fcueoff(2:end)-fcue(1:end-1); % if starts with cue off
if mode(sign(diff1))<0
    % started with cue on to off
    % use diff2
    % find differences approx. matching cue duration
    currfcue=fcue(1:end-1);
    currfcueoff=fcueoff(2:end);
    cueStarts=currfcue(diff2>cueInds-indsSlop & diff2<cueInds+indsSlop);
    cueEnds=currfcueoff(diff2>cueInds-indsSlop & diff2<cueInds+indsSlop);
else
    % use diff1
    currfcue=fcue;
    currfcueoff=fcueoff;
    cueStarts=currfcue(diff1>cueInds-indsSlop & diff1<cueInds+indsSlop);
    cueEnds=currfcueoff(diff1>cueInds-indsSlop & diff1<cueInds+indsSlop);
end

end