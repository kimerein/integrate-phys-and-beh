function highSpeedVideo_to_trials(filename,cueVarName,distractorVarName,rawVarName,wheelVarName,settings)

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
cueThresh=100;
rawThresh=1; % for where video skips, this is after the wheel has been turning for a little bit, about 0.96 s in Kim's rig
distractorThresh=100;

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
line([1 framesmax],[cueThresh cueThresh],'Color','b');
line([1 framesmax],[-cueThresh -cueThresh],'Color','b');
ylabel('Difference frame-to-frame in cue zone');
xlabel('Frames');
% find cues of specific duration
[cueStarts,cueEnds]=findEventsOfSpecificDuration(cueDiffs,cueThresh,cueInds,indsSlop);
% plot
for i=1:length(cueStarts)
    line([cueStarts(i) cueEnds(i)],[cueThresh-cueThresh*0.25 cueThresh-cueThresh*0.25],'Color','c','LineWidth',2);
end

% get wheel turns beginning at a specific time before cue
wheelOnsets=find(rawDiffs>rawThresh);
comboDiffs=zeros(size(cueDiffs));
comboDiffs(cueStarts)=-1;
% only keep wheelOnsets closest to cue onset, because cue is clearer
% signal
wheelStarts=zeros(size(cueStarts));
for i=1:length(cueStarts)
    % find closest wheel onset preceding cue
    wheelStarts(i)=wheelOnsets(find(wheelOnsets<cueStarts(i),1,'last'));
end
comboDiffs(wheelStarts)=1;
[wheelStarts,wheelEnds]=findEventsOfSpecificDuration(comboDiffs,0.5,wheelBeforeCueInds,indsSlop);
plot(rawDiffs(1:framesmax),'Color','g'); hold on;
line([1 framesmax],[rawThresh rawThresh],'Color','g');
for i=1:length(wheelStarts)
    line([wheelStarts(i) wheelEnds(i)],[cueThresh-cueThresh*0.3 cueThresh-cueThresh*0.3],'Color',[0 0.75 0],'LineWidth',2);
end

% get distractors
plot(distractorDiffs(1:framesmax),'Color','r'); hold on;
line([1 framesmax],[distractorThresh distractorThresh],'Color','r');
line([1 framesmax],[-distractorThresh -distractorThresh],'Color','r');
[distractorStarts,distractorEnds]=findEventsOfSpecificDuration(distractorDiffs,distractorThresh,distractorInds,indsSlop);
% plot
for i=1:length(distractorStarts)
    line([distractorStarts(i) distractorEnds(i)],[cueThresh-cueThresh*0.15 cueThresh-cueThresh*0.15],'Color','m','LineWidth',2);
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