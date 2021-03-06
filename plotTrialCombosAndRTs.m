function [reactionTimes,sequenceMatchStarts,forPairs_sequenceMatchStarts,RT_pairs]=plotTrialCombosAndRTs(templateSequence,metadata,trialTypes,tbt,whichReach,useAsCue,nNext,doPlot)

% have templateSequence take on a form like
% templateSequence{1}=trialsMatchingCondA;
% templateSequence{2}=nan;      nan is a wildcard indicating that 2nd trial
% in sequence can have any value
% templateSequence{3}=trialsMatchingCondB;
% templateSequence{4}=trialsMatchingCondC;

% doPlot 1 if want to plot these selected trials matching templateSequence

f=fieldnames(trialTypes);
% get matches to templateSequence
sequencesInSet=zeros(1,length(trialTypes.(f{1}))-(length(templateSequence)-1));
for i=1:length(templateSequence)
    if ~isnan(templateSequence{i})
        temp=templateSequence{i}==1;
        % make row vector
        if size(temp,1)>1
            temp=temp';
        end
        if size(temp,1)>1
            error('templateSequence{i} must be 1D vector');
        end
        sequencesInSet(1:length(trialTypes.(f{1}))-(length(templateSequence)-1))=sequencesInSet(1:length(trialTypes.(f{1}))-(length(templateSequence)-1))+temp(i:end-((length(templateSequence)-1)-(i-1)));
    else
        sequencesInSet=sequencesInSet+1; % nan in templateSequence indicates wildcard, all trial types acceptable
    end
end
% look for spans of trials where all conditions in templateSequence are
% true sequentially
sequenceMatchStarts=sequencesInSet>=length(templateSequence);

if doPlot==1
    maxTrialsToPlot=350;
    useThese=find(sequenceMatchStarts==1);
    for i=1:length(templateSequence)-1
        useThese=[useThese useThese+i];
    end
    useThese=sort(useThese);
    useThese=unique(useThese);
    if length(useThese)>maxTrialsToPlot
        startAt=randsample(length(useThese)-maxTrialsToPlot,1);
        useThese=useThese(startAt:startAt+maxTrialsToPlot);
    end
    plotEventScatter(tbt,useThese,1);
    title('templateSequence TRUE');
%     useThese=find(out==0);
%     if length(useThese)>maxTrialsToPlot
%         useThese=useThese(randsample(length(useThese),maxTrialsToPlot));
%     end
%     plotEventScatter(tbt,useThese,1);
%     title([curr_type.name ' FALSE']);
%     pause;
end

forPairs_sequenceMatchStarts=sequenceMatchStarts;
sequenceMatchStarts=[sequenceMatchStarts zeros(1,length(templateSequence)-1)];
[reactionTimes,RT_pairs]=getRTpairs_contingentOnTrialType(tbt,whichReach,useAsCue,metadata,sequenceMatchStarts,[],1,[],'contingency1',[],nNext,doPlot);

end


function plotEventScatter(tbt,showRows,doReachBatch)

timespertrial=nanmean(tbt.times,1);
showRows=showRows(~isnan(showRows));
plotSettings=plotCueTriggered_settings();
figure();
k=1;
plotfields=plotSettings.plotevents;
plotcolors=plotSettings.eventColors;
plotoutlines=plotSettings.eventOutlines;
plotthresh=plotSettings.eventThresh;
plotfirstN=plotSettings.firstN;
% Only plot reaches, cue and pellet presented events
takeThese_plotfields=zeros(1,length(plotfields));
for i=1:length(plotfields)
    if ~isempty(strfind(plotfields{i},'cue')) || ~isempty(strfind(plotfields{i},'pellet')) || ~isempty(strfind(plotfields{i},'reach'))
        takeThese_plotfields(i)=1;
        if doReachBatch==1
            currname=plotfields{i};
            if isfield(tbt,['reachBatch_' currname])
                plotfields{i}=['reachBatch_' currname];
            end
        end
    end
end
plotfields=plotfields(takeThese_plotfields==1);
plotcolors=plotcolors(takeThese_plotfields==1);
plotoutlines=plotoutlines(takeThese_plotfields==1);
plotthresh=plotthresh(takeThese_plotfields==1);
plotfirstN=plotfirstN(takeThese_plotfields==1);
for i=1:length(showRows)
    for j=1:length(plotfields)
        if ~isfield(tbt,plotfields{j})
            error([plotfields{j} ' field absent from tbt. See plotCueTriggered_settings.m to specify fields to plot.']);
        end
        currEvents=tbt.(plotfields{j});
        event_thresh=plotthresh{j};
        event_ind=find(currEvents(showRows(i),:)>event_thresh);
        n=length(event_ind);
        if ischar(plotfirstN{j})
            if strcmp('all',plotfirstN{j})
                % plot all events
                if n>500
                    event_ind=event_ind(1:10:end);
                    n=length(event_ind);
                end
            end
        else
            % plot first n events
            n=plotfirstN{j};
        end
        if isempty(event_ind)
            continue
        end
        for l=1:n
            scatter([timespertrial(event_ind(l)) timespertrial(event_ind(l))],[k k],[],'MarkerEdgeColor',plotoutlines{j},...
                'MarkerFaceColor',plotcolors{j},...
                'LineWidth',plotSettings.eventLineWidth);
            hold on;
        end
    end
    k=k+1;
end

end