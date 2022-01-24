function highSpeedVideo_to_trials(filename,cueVarName,distractorVarName,rawVarName,wheelVarName,settings)

% filename is file that contains output of getRigEvents.py
% i.e., the cue zone, wheel zone, distractor zone and raw differences
% between frames, in order to sort into trials
% 
% settings.fps = frames per second of video
% settings.distractorDuration = in seconds
% settings.cueDuration = in seconds
% settings.wheelSecBeforeCue = when wheel starts in seconds before cue

% have user set these
cueThresh=100;


a=load(filename);
cueDiffs=a.(cueVarName);
distractorDiffs=a.(distractorVarName);
rawDiffs=a.(rawVarName);
wheelDiffs=a.(wheelVarName);
cueInds=floor(settings.cueDuration * settings.fps);
distractorInds=floor(settings.distractorDuration * settings.fps);
wheelBeforeCueInds=floor(settings.wheelSecBeforeCue * settings.fps);
indsSlop=floor(0.05 * settings.fps);

% get possible cues first
% then look for wheel turns at a fixed time before cue onset
figure();
framesmax=length(cueDiffs);
if framesmax>30000
    framesmax=30000;
end
plot(cueDiffs(1:framesmax),'Color','b'); hold on;
line([1 framesmax],[cueThresh cueThresh],'Color','r');
line([1 framesmax],[-cueThresh -cueThresh],'Color','r');
ylabel('Difference frame-to-frame in cue zone');
xlabel('Frames');
% find cues of specific duration
[cueStarts,cueEnds]=findEventsOfSpecificDuration(cueDiffs,cueThresh,cueInds,indsSlop);
% plot
for i=1:length(cueStarts)
    line([cueStarts(i) cueEnds(i)],[cueThresh-cueThresh*0.25 cueThresh-cueThresh*0.25],'Color','c','LineWidth',2);
end


end

function [cueStarts,cueEnds]=findEventsOfSpecificDuration(cueDiffs,cueThresh,cueInds,indsSlop)

fcue=find(cueDiffs>cueThresh);
fcueoff=find(cueDiffs<-cueThresh);
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