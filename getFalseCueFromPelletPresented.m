function alltbt=getFalseCueFromPelletPresented(alltbt,trialTypes,metadata)

% false cue if
% pellet is presented again after trial duration
% trial is long ITI, i.e., double length
% cue does not turn on

% note that each trial starts before cue and then goes to the onset of the
% next cue, nearly exactly

% pelletPresent second time goes up from inds 250 to 350 (this is second
% pellet present)

% If pelletPresent LOW before 250 and then pelletPresent goes HIGH after
% 250, for a time duration longer than pellet beginning to enter the cue
% zone, and this is a long ITI trial, then this is false cue

% find time duration from start of pelletPresent going up to movie cue
% (i.e., cueZone_onVoff) -- this is due to the pellet moving into position
% and thus entering the pellet zone before the cue turns on
% from figure(); plot(nanmean(alltbt.cueZone_onVoff,1),'Color','b'); hold on; plot(nanmean(alltbt.pelletPresent,1),'Color','g');
% this time duration is 
% 1.085 sec
% And this exactly matches the Arduino code

falseCueHere=all(alltbt.pelletPresent(:,249)<0.1,2) & nansum(alltbt.pelletPresent(:,250:end),2)>29 & trialTypes.isLongITI==1;
% this is a slightly stricter cutoff but ensures falsecueon

% cue would have been at index 275
alltbt.falseCueOn=zeros(size(alltbt.cue));
alltbt.falseCueOn(falseCueHere==1,275)=1;

end