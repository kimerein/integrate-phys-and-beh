function [alltbt,trialTypes,metadata]=realignToDistractor(alltbt,trialTypes,metadata)

% Realign alltbt to distractor instead of cue
alltbt=realignToADistractor(alltbt,'movie_distractor');
% Change cueZone_onVoff field to be movie_distractor
alltbt.cueZone_onVoff=alltbt.movie_distractor;

end

function alltbt=realignToADistractor(alltbt,distractName)

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