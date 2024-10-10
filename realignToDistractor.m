function [alltbt,trialTypes,metadata]=realignToDistractor(alltbt,trialTypes,metadata,realignToAllDistractors)

if ~isempty(trialTypes)
    fields1=fieldnames(trialTypes);
    for i = 1:length(fields1)
        alltbt.(fields1{i}) = trialTypes.(fields1{i});
    end
end
if ~isempty(metadata)
    fields2=fieldnames(metadata);
    for i = 1:length(fields2)
        alltbt.(fields2{i}) = metadata.(fields2{i});
    end
end

% Realign alltbt to distractor instead of cue
alltbt=realignToADistractor(alltbt,'movie_distractor',realignToAllDistractors);
% Change cueZone_onVoff field to be movie_distractor
alltbt.cueZone_onVoff=alltbt.movie_distractor;

if ~isempty(trialTypes)
   for i = 1:length(fields1)
       trialTypes.(fields1{i})=alltbt.(fields1{i});
   end
end

if ~isempty(metadata)
   for i = 1:length(fields2)
       metadata.(fields2{i})=alltbt.(fields2{i});
   end
end

end

function alltbt=realignToADistractor(alltbt,distractName,realignToAllDistractors)

% realignToAllDistractors=true;
if realignToAllDistractors==true
    % Find the index of the cue in the current alignment
    [~, ma] = nanmax(nanmean(alltbt.cueZone_onVoff, 1), [], 2);
    temp = alltbt.(distractName);
    temp(temp >= 0.5) = 1; temp(temp < 0.5) = 0;
    alltbt.(distractName) = temp;
    figure(); plot(nanmean(alltbt.(distractName),1),'Color','k');
    % Initialize new alltbt structure to store expanded data
    new_alltbt = struct();
    fieldNames = fieldnames(alltbt);
    for i = 1:numel(fieldNames)
        new_alltbt.(fieldNames{i}) = [];
    end

    % List of fields to shift
    shiftFields = {'movie_distractor', 'all_reachBatch', 'cue', 'isChewing', 'isHold', 'optoOn', 'optoZone', 'pawOnWheel', 'pelletPresent', 'pelletmissingreach_reachStarts',...
       'reachBatch_all_pawOnWheel', 'reachBatch_drop_reachStarts', 'reachBatch_miss_reachStarts', 'reachBatch_success_reachStarts', 'reachStarts',...
       'reachStarts_pelletPresent', 'cueZone_onVoff'}; % Added 'cueZone_onVoff' to shift list

    % Initialize a counter for the new trials
    newTrialCounter = 0;
    
    % Loop over each trial
    for i = 1:size(temp, 1)
        if mod(i,100)==0
            disp(i);
        end
        % Find distractor onsets
        tempie = temp(i, :);
        f = find(diff(tempie) == 1);
        if isempty(f)
            continue
        end

        % For each distractor onset in the trial
        for idx = 1:length(f)
            distractor_onset = f(idx);
            shift = distractor_onset - ma; % Shift amount

            % Increment the new trial counter
            newTrialCounter = newTrialCounter + 1;

            % Loop over all fields in alltbt
            for j = 1:numel(fieldNames)
                fieldName = fieldNames{j};
                fieldData = alltbt.(fieldName)(i, :);

                % Shift fields that are in the shiftFields list
                if ismember(fieldName, shiftFields)
                    if ~isempty(fieldData)
                        if shift >= 0
                            shiftedData = [fieldData(shift + 1:end), nan(1, shift)];
                        else
                            shiftedData = [nan(1, -shift), fieldData(1:end + shift)];
                        end
                    else
                        shiftedData = nan(size(fieldData));
                    end
                else
                    % For fields not in shiftFields, copy the data as is
                    shiftedData = fieldData;
                end

                try
                    % Append the shifted data to the new alltbt structure
                    if isdatetime(shiftedData)
                        if isempty(new_alltbt.(fieldName))
                            tempda=datetime.empty;
                        else
                            tempda=new_alltbt.(fieldName);
                        end
                        tempda(newTrialCounter)=shiftedData;
                        new_alltbt.(fieldName)=tempda;
                    elseif iscell(shiftedData)
                        tempda=new_alltbt.(fieldName);
                        tempda{newTrialCounter}=shiftedData{1};
                        new_alltbt.(fieldName)=tempda;
                    else
                        new_alltbt.(fieldName)(newTrialCounter, :) = shiftedData;
                    end
                catch
                    error("problem aligning this field");
                end
            end
        end
    end

    % Replace the original alltbt with the expanded new_alltbt
    alltbt = new_alltbt;

    % Plotting for verification (you can modify or remove this section)
    hold on;
    plot(nanmean(alltbt.(distractName), 1), 'Color', 'b');
    title('Average of movie_distractor after realignment to all distractors');

else
    randomDistractor=false;

    % for each trial, randomly choose one of the distractors and align trial to
    % this instead of cue
    % where is cue currently
    [~,ma]=nanmax(nanmean(alltbt.cueZone_onVoff,1),[],2);
    temp=alltbt.(distractName);
    temp(temp>=0.5)=1; temp(temp<0.5)=0;
    alltbt.(distractName)=temp;
    figure(); plot(nanmean(alltbt.(distractName),1),'Color','k');
    if randomDistractor==true
        wasDistract=any(alltbt.(distractName)>0.5,2);
    else
        wasDistract=any(temp(:,ma:end)>0.5,2);
    end
    hold on; plot(nanmean(temp(wasDistract==1,:),1),'Color','r');
    shiftBy=nan(size(temp,1),1);
    if randomDistractor==true
        for i=1:size(temp,1)
            % find distractor onsets
            f=find(diff(temp(i,:))==1);
            if isempty(f)
                continue
            end
            f=f(randperm(length(f)));
            f=f(1);
            % realign to this
            shiftBy(i)=f-ma; % if positive, will shift backwards, else will shift forward in time
        end
    else
        % The distractor after the cue
        for i=1:size(temp,1)
            % find distractor onsets
            tempie=temp(i,:);
            tempie(1:ma)=0;
            f=find(diff(tempie)==1);
            if isempty(f)
                continue
            end
            f=f(randperm(length(f)));
            f=f(1);
            % realign to this
            shiftBy(i)=f-ma; % if positive, will shift backwards, else will shift forward in time
        end
    end
    % shift all relevant fields
    f={'movie_distractor','all_reachBatch','cue','isChewing','isHold','optoOn','optoZone','pawOnWheel','pelletPresent','pelletmissingreach_reachStarts',...
        'reachBatch_all_pawOnWheel','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','reachBatch_success_reachStarts','reachStarts',...
        'reachStarts_pelletPresent'};
    for i=1:length(f)
        if ~isfield(alltbt,f{i})
            warning([f{i} ' is not a field of alltbt']);
            continue
        end
        temp=alltbt.(f{i});
        for j=1:size(temp,1)
            if ~isnan(shiftBy(j))
                if shiftBy(j)>0
                    tempie=temp(j,:);
                    temp(j,:)=[tempie(shiftBy(j)+2:end) tempie(1:shiftBy(j)+1)]; % KR fixed 7/26/2024
                end
            else
                temp(j,:)=nan; % KR added 7/26/2024
            end
        end
        alltbt.(f{i})=temp;
    end

    hold on; plot(nanmean(alltbt.(distractName),1),'Color','b');
    temp=alltbt.(distractName);
    hold on; plot(nanmean(temp(wasDistract==1,:),1),'Color','c');
    title('black before, red only distract present trials before, blue after shift');
end

end