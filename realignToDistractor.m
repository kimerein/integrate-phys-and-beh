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
f={'all_reachBatch','cue','isChewing','isHold','movie_distractor','optoOn','optoZone','pawOnWheel','pelletPresent','pelletmissingreach_reachStarts',...
   'reachBatch_all_pawOnWheel','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','reachBatch_success_reachStarts','reachStarts',...
   'reachStarts_pelletPresent'};
for i=1:length(f)
    temp=alltbt.(f{i});
    for j=1:size(temp,1)
        if ~isnan(shiftBy(j))
            if shiftBy(j)>0
                temp(j,:)=circshift(temp(j,:),[0 shiftBy(j)]);
            end
        end
    end
    alltbt.(f{i})=temp;
end

hold on; plot(nanmean(alltbt.(distractName),1),'Color','b');
title('black before, blue after shift');

end