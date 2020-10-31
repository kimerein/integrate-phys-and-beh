function plotReactionTimeDotScatter(alltbt,cueName,reachName,successName,dropName,missName,pelletMissingName,trialTypes,metadata,useOptoMode,trialsToUse,omitEmptyTrials,sortIntoPair1vs2)

thresh=0.05;
% trialsToUse should be a vector of 0's and 1's -- will only plot trial
% pairs for which trialsToUse==1
plotReachesWithinSec=10; % plot only first reach if isempty, else plot all reaches within this many secs of first reach

if (isempty(trialsToUse) && size(alltbt.(reachName),1)>200) || (~isempty(trialsToUse) && nansum(trialsToUse==1)>100)
    answer=questdlg('More than 200 trials to plot. Subsample or plot all?','Trials to plot','Subsample','Plot all','Cancel','Cancel');
    switch answer
        case 'Subsample'
            if isempty(trialsToUse)
                 plotThese=randsample(size(alltbt.(reachName),1),200);
            else
                plotThese=randsample(nansum(trialsToUse==1),100);
                f=find(trialsToUse==1);
                plotThese=f(plotThese);
            end
        case 'Plot all'
            plotThese=1:size(alltbt.(reachName),1);
        case 'Cancel'
            return
    end
else
    if isempty(trialsToUse)
        plotThese=1:size(alltbt.(reachName),1);
    else
        plotThese=find(trialsToUse==1);
    end
end
if ~isempty(trialsToUse)
    % take both trials in pair
    backup_plotThese=plotThese;
    plotThese=[plotThese plotThese+1];
    plotThese=unique(plotThese);
    firstInPair=backup_plotThese;
    secondInPair=plotThese(~ismember(plotThese,backup_plotThese));
end

if useOptoMode==true % will plot the average timing of the opto rather than the exact timing per trial
    if isfield(trialTypes,'optoGroup')
        u=unique(trialTypes.optoGroup);
        u=u(~isnan(u));
        for j=1:length(u)
            optoStarts=nan(size(alltbt.(reachName),1),1);
            optoEnds=nan(size(alltbt.(reachName),1),1);
            for i=1:size(alltbt.(reachName),1)
                if trialTypes.optoGroup(i)==u(j)
                    temp=find(alltbt.optoOn(i,:)>thresh,1,'first');
                    if isempty(temp)
                        continue
                    end
                    optoStarts(i)=temp;
                    temp=find(alltbt.optoOn(i,optoStarts(i):end)<thresh,1,'first');
                    if isempty(temp)
                        continue
                    end
                    optoEnds(i)=temp;
                end
            end
            allStart(j)=mode(optoStarts);
            allEnd(j)=mode(optoEnds);
        end
    else
        optoStarts=nan(size(alltbt.(reachName),1),1);
        optoEnds=nan(size(alltbt.(reachName),1),1);
        for i=1:size(alltbt.(reachName),1)
            if trialTypes.led(i)==1
               temp=find(alltbt.optoOn(i,:)>thresh,1,'first');
               if isempty(temp)
                   continue
               end
               optoStarts(i)=temp;
               temp=find(alltbt.optoOn(i,optoStarts(i):end)<thresh,1,'first');
               if isempty(temp)
                   continue
               end
               optoEnds(i)=temp;
            end
        end
        allStart=mode(optoStarts);
        allEnd=mode(optoEnds);
    end
end

cueInd=find(nanmean(alltbt.(cueName),1)>thresh,1,'first');
timeStep=mode(diff(nanmean(alltbt.times,1)));
figure();
reaches=alltbt.(reachName);
successes=alltbt.(successName);
drops=alltbt.(dropName);
misses=alltbt.(missName);
pelletMissing=alltbt.(pelletMissingName);
hold on;
useTinc=1;
for incr=1:length(plotThese)
    disp(incr);
    i=plotThese(incr);
    if omitEmptyTrials==true
        yi=useTinc;
        useTinc=useTinc+1;
    else
        yi=i;
    end    
    if ~isempty(trialsToUse)
        if ismember(i,firstInPair)
            line([0 size(alltbt.optoOn,2)*timeStep],[yi yi],'Color',[0.8 0.8 0.8],'LineWidth',1);
        elseif ismember(i,secondInPair)
            line([0 size(alltbt.optoOn,2)*timeStep],[yi yi],'Color',[0.8 1 0.8],'LineWidth',1);
        end
    end
    if trialTypes.led(i)==1
        if useOptoMode==false
            startOptoOn=find(alltbt.optoOn(i,:)>thresh,1,'first');
            endOptoOn=find(alltbt.optoOn(i,startOptoOn:end)<thresh,1,'first');
        else
            if isfield(trialTypes,'optoGroup')
                currgrp=trialTypes.optoGroup(i);
                f=find(u==currgrp);
                startOptoOn=allStart(f);
                endOptoOn=allEnd(f);
            else
                startOptoOn=allStart;
                endOptoOn=allEnd;
            end
        end
        line([startOptoOn*timeStep (startOptoOn+endOptoOn)*timeStep],[yi yi],'Color',[1 0.5 0.5],'LineWidth',1);
    end
    scatter(cueInd*timeStep,yi,[],'b','fill');
    temp=reaches(i,:);
    firstReachInd=find(temp(cueInd:end)>thresh,1,'first');
    firstReachInd=cueInd+firstReachInd-1;
    if isempty(plotReachesWithinSec)
        nextReachesInds=[];
    else
        if firstReachInd+floor(plotReachesWithinSec/timeStep)>size(alltbt.optoOn,2)
            upTo=size(alltbt.optoOn,2);
        else
            upTo=firstReachInd+floor(plotReachesWithinSec/timeStep);
        end
        nextReachesInds=find(temp(firstReachInd+1:upTo)>thresh);
        nextReachesInds=firstReachInd+1+nextReachesInds-1;
    end
    allReachesInds=[firstReachInd nextReachesInds];
    for k=1:length(allReachesInds)
        firstReachInd=allReachesInds(k);
        type=nan;
        if successes(i,firstReachInd)>thresh
            type=1;
        elseif drops(i,firstReachInd)>thresh
            type=2;
        elseif misses(i,firstReachInd)>thresh
            type=3;
        elseif pelletMissing(i,firstReachInd)>thresh
            type=4;
        end
        if isnan(type)
            continue
        end
        switch type
            case 1
                scatter(firstReachInd*timeStep,yi,[],'g','fill');
            case 2
                scatter(firstReachInd*timeStep,yi,[],'r','fill');
            case 3
                scatter(firstReachInd*timeStep,yi,[],'c','fill');
            case 4
                scatter(firstReachInd*timeStep,yi,[],[0.8 0.8 0.8]);
        end
    end
end