function settings=settingsForReachRates(alltbt,nameOfCue,display)

% when does cue turn on?
timestep=mode(diff(nanmean(alltbt.times,1)));
tims=0:timestep:(size(alltbt.times,2)-1)*timestep;
[~,ma]=nanmax(nanmean(alltbt.(nameOfCue),1));
settings.maxIndNameOfCue=ma;
% figure(); plot(tims,nanmean(alltbt.cueZone_onVoff,1)); xlabel('Time (sec)'); ylabel('Cue');
disp(['Cue turns on at ' num2str(tims(ma)) ' seconds']);
cuetimeat=tims(ma);

% for external vis cue only, ramp slows down reach, and cueZone_onVoff is
% detected mid-ramp
% Note that preCueWindows are not used
% settings.preCueWindow_start1=cuetimeat-1; %+0.25; %-0.45; %0; % define start of time window from trial onset, in seconds -- for first window, assuming that trial onset is 0 sec
% settings.preCueWindow_end1=cuetimeat-0.25; %1; % define end of time window from trial onset, in seconds -- for first window
% % preCueWindow_start2=3.81; % define start of time window from trial onset, in seconds -- for second window
% % preCueWindow_end2=5.31; % define end of time window from trial onset, in seconds -- for second window
% settings.preCueWindow_start2=cuetimeat-1; % if these are same time, will not use second window %0.75; %0.5; %0.5; %8.5;  
% settings.preCueWindow_end2=cuetimeat-0.25; %1.25; %-0.25; %1; %0.95; %9.5; 
% settings.reachAfterCueWindow_start=-0.05; %-0.05; %-0.25; % after cue window in seconds, for hits, i.e., real cued reaching
% settings.reachAfterCueWindow_end=0.4; %0.4; %0.75; % after cue window in seconds, for hits, i.e., real cued reaching

% for opto cue
% settings.preCueWindow_start1=cuetimeat-2; %+0.25; %-0.45; %0; % define start of time window from trial onset, in seconds -- for first window, assuming that trial onset is 0 sec
% settings.preCueWindow_end1=cuetimeat; %1; % define end of time window from trial onset, in seconds -- for first window
% % preCueWindow_start2=3.81; % define start of time window from trial onset, in seconds -- for second window
% % preCueWindow_end2=5.31; % define end of time window from trial onset, in seconds -- for second window
% settings.preCueWindow_start2=cuetimeat+1.5; % if these are same time, will not use second window %0.75; %0.5; %0.5; %8.5;  
% settings.preCueWindow_end2=cuetimeat+9; %1.25; %-0.25; %1; %0.95; %9.5; 
% settings.reachAfterCueWindow_start=-0.05; %-0.05; %-0.25; % after cue window in seconds, for hits, i.e., real cued reaching
% settings.reachAfterCueWindow_end=0.35; %0.4; %0.75; % after cue window in seconds, for hits, i.e., real cued reaching

% for SC
settings.preCueWindow_start1=cuetimeat-0.4; %+0.25; %-0.45; %0; % define start of time window from trial onset, in seconds -- for first window, assuming that trial onset is 0 sec
settings.preCueWindow_end1=cuetimeat; %1; % define end of time window from trial onset, in seconds -- for first window
% preCueWindow_start2=3.81; % define start of time window from trial onset, in seconds -- for second window
% preCueWindow_end2=5.31; % define end of time window from trial onset, in seconds -- for second window
settings.preCueWindow_start2=cuetimeat-0.4; % if these are same time, will not use second window %0.75; %0.5; %0.5; %8.5;  
settings.preCueWindow_end2=cuetimeat; %1.25; %-0.25; %1; %0.95; %9.5; 
settings.reachAfterCueWindow_start=0; %-0.05; %-0.25; % after cue window in seconds, for hits, i.e., real cued reaching
settings.reachAfterCueWindow_end=0.4; %0.75; % after cue window in seconds, for hits, i.e., real cued reaching

if display==true
    disp('SETTINGS FOR REACH RATES');
    disp(settings);
    pause;
end

end