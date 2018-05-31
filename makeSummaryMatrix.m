function summary=makeSummaryMatrix(tbt_events,tbt_classes,tbt_metadata)

% tbt means trial-by-trial organization of data

% For each trial ...
% Trial ID, start datetime, end datetime, mouse ID, trial type, and
% events in form [[event_id, event_ms],...]

% Get trial start and end times
tbt_metadata=getTrialStartAndEnd(tbt_metadata,tbt_events);

% Save all metadata to summary
f=fieldnames(tbt_metadata);
for i=1:length(f)
    summary.(f{i})=tbt_metadata.(f{i});
end

% Convert all logical fields to single
f=fieldnames(tbt_classes);
for i=1:length(f)
    tbt_classes.(f{i})=single(tbt_classes.(f{i})==1);
    summary.(f{i})=tbt_classes.(f{i});
end

% All event fields as double
f=fieldnames(tbt_events);
for i=1:length(f)
    if ~isempty(strfind(f{i},'backup'))
        continue
    end
    summary.(f{i})=tbt_events.(f{i});
end

% Also get indices of event fields

end

function metadata=getTrialStartAndEnd(metadata,tbt_events)

% Get unwrapped trial start and end times
settings=plotCueTriggered_settings();
pointsBefore=settings.pointsFromPreviousTrial;
tookSecBeforeCue=pointsBefore*mode(diff(nanmean(tbt_events.times,1)));
% For each session, calculate unwrapped trial times
[sessid,sessStartInd]=unique(metadata.sessid);
sessStartInd=[sessStartInd; length(metadata.sessid)+1];
for i=1:length(sessid)
    curr_sessid=sessid(i);
    curr_sessStartInd=sessStartInd(i);
    trialStarts=nan(1,length(curr_sessStartInd:sessStartInd(i+1)-1));
    trialEnds=nan(1,length(curr_sessStartInd:sessStartInd(i+1)-1));
    runningSum=0;
    for j=curr_sessStartInd:sessStartInd(i+1)-1
        trialStarts(j)=runningSum;
        if ~any(~isnan(tbt_events.times(j,:)))
            trialEnds(j)=runningSum;
            continue
        end
        runningSum=runningSum+(nanmax(tbt_events.times(j,:))-nanmin(tbt_events.times(j,:)));
        trialEnds(j)=runningSum;
        runningSum=runningSum-tookSecBeforeCue; % this time window is also included as baseline in next trial
    end
    % Save trial starts and ends as datetimes
    for j=curr_sessStartInd:sessStartInd(i+1)-1
        % Subtract off total duration of session, then add back unwrapped
        % start or end time of each trial
        starts_tbtDatetime{j}=metadata.sess_datetime{j}-seconds(runningSum)+seconds(trialStarts(j));
        ends_tbtDatetime{j}=metadata.sess_datetime{j}-seconds(runningSum)+seconds(trialEnds(j));
    end
end
metadata.trial_starts=starts_tbtDatetime';
metadata.trial_ends=ends_tbtDatetime';

end