function settings=settingsForDprimes(alltbt,nameOfCue,display)

settings.preCueWindow_start1=0.529; %0; % define start of time window from trial onset, in seconds -- for first window, assuming that trial onset is 0 sec
settings.preCueWindow_end1=1.029; %1; % define end of time window from trial onset, in seconds -- for first window
% preCueWindow_start2=3.81; % define start of time window from trial onset, in seconds -- for second window
% preCueWindow_end2=5.31; % define end of time window from trial onset, in seconds -- for second window
[~,ma]=nanmax(nanmean(alltbt.(nameOfCue),1));
cuetimeat=mode(diff(nanmean(alltbt.times,1)))*ma;
settings.preCueWindow_start2=cuetimeat+9; %3; %8.5;  
settings.preCueWindow_end2=cuetimeat+9.5; %4; %9.5; 
settings.reachAfterCueWindow_start=0; %-0.25; % after cue window in seconds, for hits, i.e., real cued reaching
settings.reachAfterCueWindow_end=0.5; %0.75; % after cue window in seconds, for hits, i.e., real cued reaching

% PREVIOUSLY USED
% settings.preCueWindow_start1=0; % define start of time window from trial onset, in seconds -- for first window, assuming that trial onset is 0 sec
% settings.preCueWindow_end1=1.5; % define end of time window from trial onset, in seconds -- for first window
% % preCueWindow_start2=3.81; % define start of time window from trial onset, in seconds -- for second window
% % preCueWindow_end2=5.31; % define end of time window from trial onset, in seconds -- for second window
% [~,ma]=nanmax(nanmean(alltbt.(nameOfCue),1));
% cuetimeat=mode(diff(nanmean(alltbt.times,1)))*ma;
% settings.preCueWindow_start2=cuetimeat+1.5; 
% settings.preCueWindow_end2=cuetimeat+3;
% settings.reachAfterCueWindow_start=0; % after cue window in seconds, for hits, i.e., real cued reaching
% settings.reachAfterCueWindow_end=1.5; % after cue window in seconds, for hits, i.e., real cued reaching

if display==true
    disp(settings);
    pause;
end

end