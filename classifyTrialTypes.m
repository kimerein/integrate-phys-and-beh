function [out,tbt]=classifyTrialTypes(tbt,settings)

% User-defined settings in trialTypeSettings.m

% Evaluate boolean tests
for i=1:length(settings.bool_test)
    curr_test=settings.bool_test(i);
    out.bool_test(i).testResults=evaluateTest(tbt,curr_test,settings);
end

% Combine outcomes of boolean tests to get trial type classifications
for i=1:length(settings.trialtype)
    curr_type=settings.trialtype(i);
    out.trialtype(i).name=curr_type.name;
    out.trialtype(i).isThisType=classifyUsingBoolTests(out,curr_type,settings,0,tbt); % second-to-last argument is whether to plot
end

end

function out=classifyUsingBoolTests(bool_out,curr_type,settings,doPlot,tbt)

parts=cell(1,length(curr_type.outcomes));
useThese=zeros(1,length(curr_type.outcomes));
for i=1:length(curr_type.outcomes)
    if isnan(curr_type.outcomes(i))
        % skip this criterion
        parts{i}='';
    elseif curr_type.outcomes(i)==0
        parts{i}=[' (bool_out.bool_test(' num2str(i) ').testResults==0) '];
        useThese(i)=1;
    elseif curr_type.outcomes(i)==1
        parts{i}=[' (bool_out.bool_test(' num2str(i) ').testResults==1) ']; 
        useThese(i)=1;
    end
end
parts=parts(useThese==1);

switch curr_type.bool
    case '&'
        statementToEval=parts{1};
        for i=2:length(parts)
            statementToEval=[statementToEval '&' parts{i}];
        end
    case '|'
        statementToEval=parts{1};
        for i=2:length(parts)
            statementToEval=[statementToEval '|' parts{i}];
        end
    otherwise
        error('unrecognized value in trialtype.bool from trialTypeSettings.m');
end
out=eval(statementToEval);
disp('Evaluating the following statement');
disp(statementToEval);

if doPlot==1
    maxTrialsToPlot=350;
    useThese=find(out==1);
    if length(useThese)>maxTrialsToPlot
        useThese=useThese(randsample(length(useThese),maxTrialsToPlot));
    end
    plotEventScatter(tbt,useThese,1);
    title([curr_type.name ' TRUE']);
    useThese=find(out==0);
    if length(useThese)>maxTrialsToPlot
        useThese=useThese(randsample(length(useThese),maxTrialsToPlot));
    end
    plotEventScatter(tbt,useThese,1);
    title([curr_type.name ' FALSE']);
    pause;
end

end

function out=evaluateTest(tbt,curr_test,settings)

% i=1;
% bool_test(i).testwhat='reach batch';
% bool_test(i).fieldname='success_reachStarts';
% bool_test(i).test='first';
% bool_test(i).inEventSet='reachStarts_noPawOnWheel';
% bool_test(i).thresh=settings.lowThresh;
% bool_test(i).comparator='>';
% bool_test(i).window='after_cue';

getTimeWindowInd=nan;
for i=1:length(settings.timeWindow)
    currWindow=settings.timeWindow(i);
    if strcmp(currWindow.name,curr_test.window)
        getTimeWindowInd=i;
        break
    end
end
if isnan(getTimeWindowInd)
    error(['did not specify time window ' curr_test.window ' in trialTypeSettings.m']);
end
timeWindow=settings.timeWindow(getTimeWindowInd);

% Convert time window wrt cue onset into indices into data
cueInd=find(nanmean(tbt.(settings.nameOfCue),1)>settings.lowThresh,1,'first');
startInds=floor(abs(timeWindow.start)/mode(diff(nanmean(tbt.times,1))));
if timeWindow.start<0
    startInds=-startInds;
end
endInds=floor(abs(timeWindow.end)/mode(diff(nanmean(tbt.times,1))));
if timeWindow.end<0
    endInds=-endInds;
end
startInds=cueInd+startInds;
endInds=cueInd+endInds;

switch curr_test.testwhat
    case 'single reach'
        data=tbt.(curr_test.fieldname);
        testingthisfield=curr_test.fieldname;
    case 'reach batch'
        data=tbt.(['reachBatch_' curr_test.fieldname]);  
        testingthisfield=['reachBatch_' curr_test.fieldname];
    otherwise
        error(['do not recognize ' curr_test.fieldname ' as setting for bool_test from trialTypeSettings.m']);
end

if startInds<1
    startInds=1;
elseif startInds>size(data,2)
    startInds=size(data,2);
end
if endInds<1
    endInds=1;
elseif endInds>size(data,2)
    endInds=size(data,2);
end

switch curr_test.test
    case 'any'
        out=evaluateAnyTest(data(:,startInds:endInds),curr_test,settings,2);
        disp(['Determining whether any data in field tbt.' testingthisfield ' occurring between ' num2str(timeWindow.start) ' and ' num2str(timeWindow.end) ' from cue onset is ' curr_test.comparator ' ' num2str(curr_test.thresh)]);
    case 'first'
        % Check whether first in event set is of indicated type
        largerSet=tbt.(curr_test.inEventSet);
        out=isFirstEventThis(data(:,startInds:endInds),largerSet(:,startInds:endInds),curr_test,settings);
        out=out';
        disp(['Determining whether first event ' curr_test.comparator ' ' num2str(curr_test.thresh) ' in each trial in data in field tbt.' testingthisfield ' occurring between ' num2str(timeWindow.start) ' and ' num2str(timeWindow.end) ' from cue onset is' ...
                      ' at the same index location as first event in tbt.' curr_test.inEventSet ' ' curr_test.comparator ' ' num2str(curr_test.thresh) ' also occurring between ' num2str(timeWindow.start) ' and ' num2str(timeWindow.end) ' from cue onset']);
    otherwise
        % let 'any' be default mode
        out=evaluateAnyTest(data(:,startInds:endInds),curr_test,settings,2);
end

end

function isThisEventType=isFirstEventThis(theseEvents,largerEventSet,curr_test,settings)

isThisEventType=zeros(1,size(theseEvents,1));
for i=1:size(theseEvents,1)
   f=find(eval(['theseEvents(i,:) ' curr_test.comparator num2str(curr_test.thresh)]),1,'first');
   f2=find(eval(['largerEventSet(i,:) ' curr_test.comparator num2str(curr_test.thresh)]),1,'first');
   if f==f2
       isThisEventType(i)=1;
   end
end

end

function out=evaluateAnyTest(data,curr_test,settings,dim)

out=eval(['any(data' curr_test.comparator num2str(curr_test.thresh) ',' num2str(dim) ')']);

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