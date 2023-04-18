function settings=settingsForDprimes(alltbt,nameOfCue,display)

% when does cue turn on?
timestep=mode(diff(nanmean(alltbt.times,1)));
tims=0:timestep:(size(alltbt.times,2)-1)*timestep;
[~,ma]=nanmax(nanmean(alltbt.cueZone_onVoff,1));
% figure(); plot(tims,nanmean(alltbt.cueZone_onVoff,1)); xlabel('Time (sec)'); ylabel('Cue');
disp(['Cue turns on at ' num2str(tims(ma)) ' seconds']);
cuetimeat=tims(ma);

settings.preCueWindow_start1=cuetimeat-0.14-1.4; %0; % define start of time window from trial onset, in seconds -- for first window, assuming that trial onset is 0 sec
settings.preCueWindow_end1=cuetimeat-0.14-0.7; %1; % define end of time window from trial onset, in seconds -- for first window
% preCueWindow_start2=3.81; % define start of time window from trial onset, in seconds -- for second window
% preCueWindow_end2=5.31; % define end of time window from trial onset, in seconds -- for second window
settings.preCueWindow_start2=cuetimeat+8.3; %3; %8.5;  
settings.preCueWindow_end2=cuetimeat+9; %4; %9.5; 
settings.reachAfterCueWindow_start=-0.2; %-0.25; % after cue window in seconds, for hits, i.e., real cued reaching
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