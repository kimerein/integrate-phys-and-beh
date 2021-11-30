function [spikes,inTrial_spikes]=organizeSpikesToMatch_physiology_tbt(spikes,physiology_tbt)

spikes.trials=nan(size(spikes.trials));
spikes.spiketimes=nan(size(spikes.trials));
for i=1:size(physiology_tbt.phys_timepoints,1)-1
    % exclude all zeros
    temp=physiology_tbt.phys_timepoints(i,:);
    temp_nexttrial=physiology_tbt.phys_timepoints(i+1,:);
    temp(temp==0)=nan;
    temp_nexttrial(temp_nexttrial==0)=nan;
    % label all spikes that belong to this trial
    spikes.trials(spikes.unwrapped_times>=nanmin(temp) & spikes.unwrapped_times<nanmin(temp_nexttrial))=i;
    % get spikes.spiketimes relative to each trial start
    spikes.spiketimes(spikes.unwrapped_times>=nanmin(temp) & spikes.unwrapped_times<nanmin(temp_nexttrial))=spikes.unwrapped_times(spikes.unwrapped_times>=nanmin(temp) & spikes.unwrapped_times<nanmin(temp_nexttrial))-nanmin(temp);
end
temp=physiology_tbt.phys_timepoints(end,:);
temp(temp==0)=nan;
spikes.trials(spikes.unwrapped_times>=nanmin(temp) & spikes.unwrapped_times<nanmax(temp))=size(physiology_tbt.phys_timepoints,1);
spikes.spiketimes(spikes.unwrapped_times>=nanmin(temp) & spikes.unwrapped_times<nanmax(temp))=spikes.unwrapped_times(spikes.unwrapped_times>=nanmin(temp) & spikes.unwrapped_times<nanmax(temp))-nanmin(temp);

spikes.sweeps.trials=1:size(physiology_tbt.phys_timepoints,1);

% filter out spikes that occur outside of physiology_tbt
f=fieldnames(spikes);
for i=1:length(f)
    if strcmp(f{i},'params') || strcmp(f{i},'info') || strcmp(f{i},'labels') || strcmp(f{i},'sweeps')
        inTrial_spikes.(f{i})=spikes.(f{i});
    elseif strcmp(f{i},'waveforms')
        temp=spikes.(f{i});
        inTrial_spikes.(f{i})=temp(~isnan(spikes.trials),:,:);
    else
        temp=spikes.(f{i});
        inTrial_spikes.(f{i})=temp(~isnan(spikes.trials));
    end
end