function plotHallmarksOfSatiety(rr,dataset,alltbt,metadata,trialTypes)

maxTrialLength=20; % in sec
binByNTrials=20;

% rr is reachrates from plotChangeInReachProbability_fromRTdataset.m
matchesEventCond_trial_n=dataset.realDistributions.event_isSeq{1}==1;
trialLengths=nanmax(alltbt.times,[],2);
trialLengths(trialLengths>maxTrialLength)=maxTrialLength;

% check whether sessids are unique, else make them unique
u=unique(metadata.mouseid);
alreadyUsedSessIDs=[];
sessIDsAreUnique=true;
for i=1:length(u)
    currmouse=u(i);
    if any(ismember(alreadyUsedSessIDs,unique(metadata.sessid(metadata.mouseid==currmouse))))
        sessIDsAreUnique=false;
        break
    else
        alreadyUsedSessIDs=[alreadyUsedSessIDs; unique(metadata.sessid(metadata.mouseid==currmouse))];
    end
end
if sessIDsAreUnique==true
    % use sessids
else
    % make unique sessids
    u=unique(metadata.mouseid);
    j=0;
    for i=1:length(u)
        metadata.sessid(metadata.mouseid==u(i))=metadata.nth_session(metadata.mouseid==u(i))+j;
        j=j+nanmax(metadata.sessid(metadata.mouseid==u(i)));
    end
end

% get total length of each session
u=unique(metadata.sessid);
sessionLengths=nan(1,length(u));
for i=1:length(u)
    sessionLengths(i)=nansum(trialLengths(metadata.sessid==u(i)));
end

cmap=colormap('cool');
cmap_steps=floor(linspace(1,size(cmap,1),100));
figure();
bincounter=1;
rr_count=0;
avrr=[];
pelletsEaten=[];
timesintosess=[];
for i=1:size(rr.alltrials_cued,1)
    for j=1:size(rr.alltrials_cued,2)
        % for each session, for each trial
        % get reach rate on this trial
        % get fraction through session
        % get time into session
        % get number of pellets consumed so far
        curr_rr=nanmean([rr.alltrials_cued(i,j) rr.alltrials_uncued(i,j)]);
        if j+1>size(rr.fracsThroughSess,2)
            nextFrac=rr.fracsThroughSess(i,j);
        else
            nextFrac=rr.fracsThroughSess(i,j+1);
        end
        frac=rr.fracsThroughSess(i,j);
        if isnan(nextFrac) || bincounter==binByNTrials
            rr_count=rr_count+curr_rr;
            rr_count=rr_count/bincounter;
            bincounter=1;
        else
            rr_count=rr_count+curr_rr;
            bincounter=bincounter+1;
            continue
        end
        if isnan(frac)
            rr_count=0;
            bincounter=1;
            continue
        end
        timeIntoSess=frac*sessionLengths(i);
        % get approx which trial in this session
        currSess=u(i);
        nTrialsInSess=sum(metadata.sessid==currSess);
        approxCurrTrial=floor(frac*nTrialsInSess);
        if approxCurrTrial<1
            approxCurrTrial=1;
        end
        f=find(metadata.sessid==currSess);
        pelletsSoFar=nansum(trialTypes.consumed_pellet(f(1:approxCurrTrial)));
        co=ceil(pelletsSoFar)+1;
        if co>length(cmap_steps)
            co=length(cmap_steps);
        end
        scatter(timeIntoSess,rr_count,[],cmap(cmap_steps(co),:));
        avrr=[avrr rr_count];
        pelletsEaten=[pelletsEaten pelletsSoFar];
        timesintosess=[timesintosess timeIntoSess];
        hold on;
        rr_count=0;
    end
end
xlabel('Time into session (seconds)');
ylabel('Average reach rate (1/sec)');

figure();
scatter(pelletsEaten,avrr);

figure();
scatter(timesintosess,avrr);

end
        
        
        