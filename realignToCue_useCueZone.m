function realign_tbt=realignToCue_useCueZone(tbt,useAsCue,cueDuration,settings)

beginningOrEnd=false; % if want to just shift according to guess of whether beginning or end
% otherwise will actually try to align to cue onset using the derivative of
% the cue zone

if isempty(settings)
    settings.lowThresh=0.05;
    settings.isOrchestra=false;
end

lowThresh=settings.lowThresh;

% was each cue detection at the beginning or end of cue?
cueDurationInds=floor(cueDuration/mode(diff(nanmean(tbt.times,1))));

% line up cues

% find first cue ind on
cue=tbt.(useAsCue);
[~,cueind]=max(nanmean(cue,1));

cueZone=tbt.cueZone;

f=fieldnames(tbt);
for i=1:length(f)
    realign_tbt.(f{i})=nan(size(tbt.(f{i})));
end

realign_check=nan(size(cue,1),2*cueDurationInds+1);
fi=nan(1,size(cue,1));
for i=1:size(cue,1)
    temp=find(cue(i,:)>lowThresh,1,'first');
    % if cue is missing from this trial, drop this trial
    if isempty(temp)
        fi(i)=nan;
    else
        fi(i)=temp;
        % what does cue zone look like surrounding this point?
        if temp-2*cueDurationInds<1 || temp+cueDurationInds>size(cueZone,2)
            continue
        end

        if beginningOrEnd==false
            % Actually realign to cue onset
            % Use positive peak of derivative
            d=diff(cueZone(i,temp-2*cueDurationInds:temp+cueDurationInds));
            [~,maxie]=nanmax(d,[],2);
            fi(i)=temp-2*cueDurationInds+maxie;
        else
            % Just check if beginning or end of cue
            if nanmean(cueZone(i,temp-cueDurationInds:temp))<nanmean(cueZone(i,temp:temp+cueDurationInds)) % this is beginning of cue
                % leave alone
            else % this is end of cue
                fi(i)=temp-(cueDurationInds-2);
            end
        end

        if fi(i)-cueDurationInds<1 || fi(i)+cueDurationInds>size(cueZone,2)
            continue
        end
        realign_check(i,:)=cueZone(i,fi(i)-cueDurationInds:fi(i)+cueDurationInds); % this is with fixed temp
    end
end

% realign cue and rest of fields
for i=1:length(f)
    currfield=tbt.(f{i});
    newfield=nan(size(currfield));
    if size(currfield,1)~=size(tbt.(useAsCue),1)
        % skip this
        continue
    end
    for j=1:size(cue,1)
        % realign each trial
        if isnan(fi(j))
            % exclude this trial
            % nan out
            temp=nan(size(currfield(j,:)));
        elseif fi(j)==cueind
            % already aligned
            temp=currfield(j,:);
        elseif fi(j)<cueind
            % shift back in time
            temp=[nan(1,cueind-fi(j)) currfield(j,1:end-(cueind-fi(j)))];
        elseif fi(j)>cueind
            % shift forward in time
            temp=[currfield(j,1+(fi(j)-cueind):end) nan(1,fi(j)-cueind)];
        end
        newfield(j,:)=temp;
    end
    % save into tbt
    realign_tbt.(f{i})=newfield;
end
            
% check alignment
if settings.isOrchestra==0
    figure(); 
    plot(realign_tbt.(useAsCue)');
    title('Re-aligned cues');
end

% Make cue uniform across trials, now that aligned
av=nanmean(realign_tbt.(useAsCue),1);
ma=max(av);
av(av<ma)=0;
realign_tbt.(useAsCue)=repmat(av,size(realign_tbt.(useAsCue),1),1);

if settings.isOrchestra==0
    figure();
    plot(realign_check');
    title('Re-aligned cue zone');
end
end