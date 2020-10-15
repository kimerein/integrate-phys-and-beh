function alltbt=getTimesWrtSessionStart(alltbt,metadata)

alltbt.timesFromSessionStart=nan(size(alltbt.timesfromarduino));

timeStep=mode(mode(nanmean(alltbt.times(:,2:end),1)-nanmean(alltbt.times(:,1:end-1),1)));

currTime=0;
currSess=0;
for i=1:size(alltbt.timesfromarduino,1)
    temp=alltbt.timesfromarduino(i,:);
    tempSess=metadata.sessid(i);
    if i==1
        currSess=tempSess;
    elseif tempSess~=currSess
        currTime=0; % reset to 0 for new session
        currSess=tempSess;
    end
    mi=nanmin(temp);
    ma=nanmax(temp);
    if isnan(mi)
        mi=0;
    end
    if isnan(ma)
        ma=timeStep*nansum(~isnan(alltbt.reachStarts(i,:)));
    end
    alltbt.timesFromSessionStart(i,:)=temp-mi+currTime;
    currTime=currTime+(ma-mi);
    if isnan(currTime)
        currTime=0;
    end
end

% Plot first 5 sessions only
u=unique(metadata.sessid);
for i=1:5
    figure();
    plot(alltbt.timesFromSessionStart(metadata.sessid==u(i),:)');
end

