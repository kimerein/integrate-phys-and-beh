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
    behphys_tbt=takeOnlyTheseTrials(behphys_tbt,ref_from_photo_to_phys);
else
    % no neural data, only behavior data
end
[tbt, tbt_photo, signalSpikeFields, signalPhotoFields]=addTimesWrtFirstTrial(tbt, tbt_photo, signalSpikeFields, signalPhotoFields);

% Add signals and names of signals


% Trials
Trials=(1:size(tbt.cue,1))';
% Cue onset
Cue_onset=
% Pellet presentation wheel begins to move
TrialStart=
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

function [tbt, tbt_photo, signalSpikeFields, signalPhotoFields]=addTimesWrtFirstTrial(tbt, tbt_photo, signalSpikeFields, signalPhotoFields)

% have already trimmed 



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