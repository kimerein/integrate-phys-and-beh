function bernardoTrialTable(tbt, bin1, bin2, bin3, signalSpikeFields, spikeFieldsNames, tbt_photo, signalPhotoFields, photoFieldsNames)

% Time delay Arduino to cue onset


% Names of columns:
varNames={'Trials', 'TrialStart', 'TrialStart_s', 'Cue_onset', 'Opto_onset', 'Opto_offset', ...
          'First_reach_timing', 'First_reach_type', 'Bin1_start', 'Bin1_end', 'N_reaches_bin1', ...
          'Bin2_start', 'Bin2_end', 'N_reaches_bin2', 'Bin3_start', 'Bin3_end', 'N_reaches_bin3'};

% Re-expand alignment into single row for signals
% Using alignment of photometry and physiology
% and times in tbt's
backup_tbt=[];
if ~isempty(signalPhotoFields) && isempty(signalSpikeFields)
    % this experiment had photometry but not physiology
    backup_tbt=tbt;
    tbt=tbt_photo; % use behavior tbt associated with photometry
elseif isempty(signalPhotoFields) && ~isempty(signalSpikeFields)
    % this experiment had physiology but not photometry
    % tbt is behavior tbt associated with physiology
elseif ~isempty(signalPhotoFields) && ~isempty(signalSpikeFields)
    % this experiment had both
    % only populate bernardoTrialTable for trials with both phys and photo
    % throw out trials that don't exist in both
    [tbt, tbt_photo, tookForPhys, tookForPhoto]=trimTbtsToMatch(tbt, tbt_photo);
    signalSpikeFields=takeOnlyTheseTrials(signalSpikeFields,tookForPhys);
    signalPhotoFields=takeOnlyTheseTrials(signalPhotoFields,tookForPhoto);
else
    % no neural data, only behavior data
end
[tbt, tbt_photo, signalSpikeFields, signalPhotoFields]=addTimes(tbt, tbt_photo, signalSpikeFields, signalPhotoFields);

% Add signals and names of signals
[signal_cue,signal_cue_times,signal_cue_times_wrt_thisTrialCueStart,signal_which_trial]=getLinedUpCueSignal(tbt, signalSpikeFields);
% Now fill in other signals
nextSignalTbt=tbt.all_reachBatch;
nextSignalTimes=tbt.times_wrt_trial_start;
nextSignalCue=tbt.cueZone_onVoff;
sig=fillInOtherSignals(nextSignalTbt,nextSignalTimes,nextSignalCue,signal_cue_times_wrt_thisTrialCueStart,signal_which_trial);
signal_ind=1;
signals{signal_ind}=sig;
signal_names{signal_ind}='all_reachBatch';
% then go through phys and photo signals, adding them
% single units
for i=1:length(spikeFieldsNames)
    nextSignalTbt=signalSpikeFields.(spikeFieldsNames{i});
    nextSignalTimes=signalSpikeFields.times_wrt_trial_start;
    nextSignalCue=signalSpikeFields.cue;
    sig=fillInOtherSignals(nextSignalTbt,nextSignalTimes,nextSignalCue,signal_cue_times_wrt_thisTrialCueStart,signal_which_trial);
    signal_ind=signal_ind+1;
    signals{signal_ind}=sig;
    signal_names{signal_ind}=spikeFieldsNames{i};
end
% photometry channels
for i=1:length(photoFieldsNames)
    nextSignalTbt=signalPhotoFields.(photoFieldsNames{i});
    nextSignalTimes=signalPhotoFields.times_wrt_trial_start;
    nextSignalCue=signalPhotoFields.cue;
    sig=fillInOtherSignals(nextSignalTbt,nextSignalTimes,nextSignalCue,signal_cue_times_wrt_thisTrialCueStart,signal_which_trial);
    signal_ind=signal_ind+1;
    signals{signal_ind}=sig;
    signal_names{signal_ind}=photoFieldsNames{i};
end


% Trials
Trials=(1:size(tbt.cue,1))';
% Cue onset
%Cue_onset=
% Pellet presentation wheel begins to move
%TrialStart=
% Trial start seconds


% Opto manipulation onset
% Opto manipulation end
% First reach batch timing
% First reach batch type
% Bin 1 start
% Bin 1 end
% Number of reaches bin 1 (before the cue)
% Bin 2 start
% Bin 2 end
% Number of reaches bin 2 (after the cue)
% Bin 3 start
% Bin 3 end
% Number of reaches bin 3 (long time after the cue)

% Convert times to indices into signals




end

function [tbt, tbt_photo, signalSpikeFields, signalPhotoFields]=addTimes(tbt, tbt_photo, signalSpikeFields, signalPhotoFields)

% have already trimmed all tbt's to have same number of rows in each field
if ~isempty(tbt)
    if ~isfield(tbt,'times_wrt_trial_start') 
        tbt=addTimesWrtTrialStarts(tbt,'times','cueZone_onVoff');
    end
end
if ~isempty(tbt_photo)
    if ~isfield(tbt_photo,'times_wrt_trial_start')
        tbt_photo=addTimesWrtTrialStarts(tbt_photo,'times','cueZone_onVoff');
    end
end
if ~isempty(signalSpikeFields)
    if ~isfield(signalSpikeFields,'times_wrt_trial_start')
        signalSpikeFields=addTimesWrtTrialStarts(signalSpikeFields,'phys_timepoints','cue');
    end
end
if ~isempty(signalPhotoFields)
    if ~isfield(signalPhotoFields,'times_wrt_trial_start')  
        signalPhotoFields=addTimesWrtTrialStarts(signalPhotoFields,'green_time','cue;);
    end
end

end

function [signal_cue,signal_cue_times,signal_cue_times_wrt_thisTrialCueStart,signal_which_trial]=getLinedUpCueSignal(tbt, signalSpikeFields)

% could use timesfromarduino to line trials back up 
% or movieframeinds*timestep sec
temp=tbt.times_wrt_trial_start; 
temp=temp'; 
temp2=diff(temp(1:end)); 
timestep=mode(temp2(temp2>0));
% each signal has its own time
% use cues to align across signals
% Plus need to resample each signal
% I think the best way to do this is going to be to put the cues at fixed
% times, then just fill in the other time points using the value from the
% closest time in the existing data structure -- yes

% If is physiology, use cue times from physiology
% Else use cue times from behavior tbt
if ~isempty(signalSpikeFields)
    [signal_cue,signal_cue_times,signal_cue_times_wrt_thisTrialCueStart,signal_which_trial]=lineUpCueChunks(signalSpikeFields, 'cue', 'cue_times');
else
    [signal_cue,signal_cue_times,signal_cue_times_wrt_thisTrialCueStart,signal_which_trial]=lineUpCueChunks(tbt, 'cueZone_onVoff', 'times_wrt_trial_start');
    signal_cue_times=0:timestep:(length(signal_cue_times)-1)*timestep;
end


end

function sig=fillInOtherSignals(nextSignalTbt,nextSignalTimes,nextSignalCue,signal_cue_times_wrt_thisTrialCueStart,signal_which_trial)

sig=nan(size(signal_cue_times_wrt_thisTrialCueStart)); % must be the same size as other signals
% just use nearest neighbor to interpolate or downsample
for i=1:size(nextSignalTbt,1)
    temp=nextSignalTbt(i,:);
    times=nextSignalTimes(i,:);
    tempcue=nextSignalCue(i,:);
    % relative to the timing of the cue in this trial
    % resample temp
    f=find(signal_which_trial==i);
    timeswrtcuestart=signal_cue_times_wrt_thisTrialCueStart(f);
    for j=1:length(timeswrtcuestart)
        fcue=find(tempcue>0.5,1,'first');
        timeswrtcuestart_nextsignal=times-times(fcue);
        % find nearest neighbor
        [~,mi]=nanmin(abs(timeswrtcuestart_nextsignal-timeswrtcuestart(j)));
        sig(f(j))=temp(mi);
    end
end

end

function [signal_cue,signal_cue_times,signal_cue_times_wrt_thisTrialCueStart,signal_which_trial]=lineUpCueChunks(tbt, cuefield, cuetimes)

cue=tbt.(cuefield);
cue_times=tbt.(cuetimes);
% Each chunk begins at the start of each cue onset, chunk continues until
% the start of the next chunk onset
signal_cue=[];
signal_cue_times=[];
signal_cue_times_wrt_thisTrialCueStart=[];
signal_which_trial=[];
for i=1:size(cue,1)
    % find first cue onset 
    temp=cue(i,:);
    f=find(temp>0.5,1,'first');
    % find index before next cue onset (or last index of trial)
    f2=find(temp<0.5,1,'last');
    if i==1
        f=1; % include ITI before first trial
    end
    chunk=temp(f:f2);
    signal_cue=[signal_cue chunk];
    signal_cue_times=[signal_cue_times cue_times(i,f:f2)];
    signal_cue_times_wrt_thisTrialCueStart=[signal_cue_times_wrt_thisTrialCueStart cue_times(i,f:f2)-cue_times(i,f)];
    signal_which_trial=[signal_which_trial i*ones(size(chunk))];
end

end

function tbt=addTimesWrtTrialStarts(tbt,whichTimesField,cuename)

if ~isempty(whichTimesField)
    timesfield=tbt.(whichTimesField);
else
    % look for a times field
    % find it automatically by searching for token 'times'
    f=fieldnames(tbt);
    for i=1:length(f)
        if ~isempty(regexp(f{i},'times'))
            timesfield=tbt.(f{i});
            break
        end
    end
end

[~,f]=max(mean(tbt.(cuename),1,'omitnan'),[],2);
tbt.times_wrt_trial_start=timesfield-repmat(min(timesfield(:,1:f),[],2,'omitnan'),1,size(timesfield,2));

end

function [behphys_tbt, behphoto_tbt, tookForPhys, tookForPhoto]=trimTbtsToMatch(behphys_tbt, behphoto_tbt)

foundReference=false;
if isfield(behphys_tbt, 'this_is_which_beh')
    if behphys_tbt.this_is_which_beh==1
        if isfield(behphys_tbt, 'reference_into_beh2trialinds')
            if any(any(~isnan(behphys_tbt.reference_into_beh2trialinds),2))
                foundReference=true;
                % use this set of trials and reference
                ref_from_phys_to_photo=behphys_tbt.reference_into_beh2trialinds(:,1);
                % throw out trials that are missing
                tookForPhys=~isnan(ref_from_phys_to_photo);
                behphys_tbt=takeOnlyTheseTrials(behphys_tbt,~isnan(ref_from_phys_to_photo));
                % throw out trials that are missing 
                tookForPhoto=ref_from_phys_to_photo;
                behphoto_tbt=takeOnlyTheseTrials(behphoto_tbt,ref_from_phys_to_photo);
            else
                % this field is missing reference
            end
        end
    elseif behphys_tbt.this_is_which_beh==2
        if isfield(behphys_tbt, 'reference_into_beh1trialinds')
            if any(any(~isnan(behphys_tbt.reference_into_beh1trialinds),2))
                foundReference=true;
                % use this set of trials and reference
                ref_from_phys_to_photo=behphys_tbt.reference_into_beh1trialinds(:,1);
                % throw out trials that are missing 
                tookForPhys=~isnan(ref_from_phys_to_photo);
                behphys_tbt=takeOnlyTheseTrials(behphys_tbt,~isnan(ref_from_phys_to_photo));
                % throw out trials that are missing 
                tookForPhoto=ref_from_phys_to_photo;
                behphoto_tbt=takeOnlyTheseTrials(behphoto_tbt,ref_from_phys_to_photo);
            else
                % this field is missing reference
            end
        end
    end
elseif isfield(behphoto_tbt, 'this_is_which_beh')
    if behphoto_tbt.this_is_which_beh==1
        if isfield(behphoto_tbt, 'reference_into_beh2trialinds')
            if any(any(~isnan(behphoto_tbt.reference_into_beh2trialinds),2))
                foundReference=true;
                % use this set of trials and reference
                ref_from_photo_to_phys=behphoto_tbt.reference_into_beh2trialinds(:,1);
                % throw out trials that are missing
                tookForPhoto=~isnan(ref_from_photo_to_phys);
                behphoto_tbt=takeOnlyTheseTrials(behphoto_tbt,~isnan(ref_from_photo_to_phys));
                % throw out trials that are missing 
                tookForPhys=ref_from_photo_to_phys;
                behphys_tbt=takeOnlyTheseTrials(behphys_tbt,ref_from_photo_to_phys);
            else
                % this field is missing reference
            end
        end
    elseif behphoto_tbt.this_is_which_beh==2
        if isfield(behphoto_tbt, 'reference_into_beh1trialinds')
            if any(any(~isnan(behphoto_tbt.reference_into_beh1trialinds),2))
                foundReference=true;
                % use this set of trials and reference
                ref_from_photo_to_phys=behphoto_tbt.reference_into_beh1trialinds(:,1);
                % throw out trials that are missing 
                tookForPhoto=~isnan(ref_from_photo_to_phys);
                behphoto_tbt=takeOnlyTheseTrials(behphoto_tbt,~isnan(ref_from_photo_to_phys));
                % throw out trials that are missing 
                tookForPhys=ref_from_photo_to_phys;
                behphys_tbt=takeOnlyTheseTrials(behphys_tbt,ref_from_photo_to_phys);
            else
                % this field is missing reference
            end
        end
    end
end
if foundReference==false
    error('Could not find reference to align behavior data for photometry and physiology');
end
if size(behphys_tbt.cue,1)~=size(behphoto_tbt.cue,1)
    error('Tbts must be same size');
end

end

function tbt=takeOnlyTheseTrials(tbt,takeThese)

f=fieldnames(tbt);
for i=1:length(f)
    temp=tbt.(f{i});
    temp=temp(takeThese(~isnan(takeThese)),:);
    tbt.(f{i})=temp;
end

end