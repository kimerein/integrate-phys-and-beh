function settings=settingsForBiasTerm(alltbt,nameOfCue,display)

% when does cue turn on?
timestep=mode(diff(nanmean(alltbt.times,1)));
tims=0:timestep:(size(alltbt.times,2)-1)*timestep;
[~,ma]=nanmax(nanmean(alltbt.(nameOfCue),1));
% figure(); plot(tims,nanmean(alltbt.cueZone_onVoff,1)); xlabel('Time (sec)'); ylabel('Cue');
disp(['Cue turns on at ' num2str(tims(ma)) ' seconds']);
cuetimeat=tims(ma);

settings.preCueWindow_start1=cuetimeat-0.98-0.4; % define start of time window from trial onset, in seconds -- for first window, assuming that trial onset is 0 sec
settings.preCueWindow_end1=cuetimeat-0.98; % define end of time window from trial onset, in seconds -- for first window
settings.preCueWindow_start2=cuetimeat+9; % if these are same time, will not use second window 
settings.preCueWindow_end2=cuetimeat+9;
settings.reachAfterCueWindow_start=-0.98; 
settings.reachAfterCueWindow_end=-0.98+0.4;

if display==true
    disp('SETTINGS FOR BIAS TERM')
    disp(settings);
    pause;
end

end