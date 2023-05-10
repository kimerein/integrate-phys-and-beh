function alignHighSpeedEventsToLowSpeedEvents(lowspeed_tbt,location_of_rig_events,fps)

dsby=10;

% timestep on Kim's rig is usually 0.0039, Fs = 256 fps
timestep=1/fps;
ds_timestep=timestep*dsby;

load(location_of_rig_events);

% Make high speed tbt, downsampling 10 times
% This is just for initial alignment
% Then will precisely align to each cue
distractor=zeros(1,ceil(nanmax(cueDiffEvs_minus/dsby)));
cue=zeros(1,ceil(nanmax(cueDiffEvs_minus/dsby)));
wheel=zeros(1,ceil(nanmax(cueDiffEvs_minus/dsby)));

cuediff=floor(cueDiffEvs_plus./dsby);
distractordiff=floor(distractorDiffEvs_plus./dsby);
wheeldiff=floor(wheelDiffEvs_plus./dsby);

cue(cuediff)=1;
distractor(distractordiff)=1;
wheel(wheeldiff)=1;

t=nanmean(lowspeed_tbt.times_wrt_trial_start,1);
[~,f]=nanmax(nanmean(lowspeed_tbt.cueZone_onVoff,1));
cuedelay=t(f);
cueindsbefore=floor(cuedelay/ds_timestep);
trialLength=nanmax(nanmean(lowspeed_tbt.times_wrt_trial_start,1));
highspeed_tbt.cue=zeros(length(cuediff),floor(trialLength/ds_timestep));
highspeed_tbt.distractor=zeros(length(cuediff),floor(trialLength/ds_timestep));
highspeed_tbt.wheel=zeros(length(cuediff),floor(trialLength/ds_timestep));
for i=1:length(cuediff)
    currcueind=cuediff(i);
    tempinds=currcueind-cueindsbefore:currcueind-cueindsbefore+size(highspeed_tbt.cue,2)-1;
    addnanatfront=0;
    addnanatback=0;
    if tempinds(1)<1
        addnanatfront=nansum(tempinds<1);
        tempinds=tempinds(find(tempinds>=1,1,'first'):end);
    end
    if tempinds(end)>length(cue)
        addnanatback=nansum(tempinds>length(cue));
        tempinds=tempinds(1:find(tempinds<=length(cue),1,'last'));
    end
    highspeed_tbt.cue(i,:)=[nan(1,addnanatfront) cue(tempinds) nan(1,addnanatback)];
    highspeed_tbt.distractor(i,:)=[nan(1,addnanatfront) distractor(tempinds) nan(1,addnanatback)];
    highspeed_tbt.wheel(i,:)=[nan(1,addnanatfront) wheel(tempinds) nan(1,addnanatback)]; 
end

end