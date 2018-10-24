function [out,tbt]=classifyTrialTypes(tbt)

% User-defined settings in trialTypeSettings.m
settings=trialTypeSettings();

% Find and convert reach batches
tbt=findReachBatches(tbt,settings.lowThresh,0); % last argument is 1 if plot, else 0

% Evaluate boolean tests
for i=1:length(settings.bool_test)
    curr_test=settings.bool_test(i);
    out.bool_test(i).testResults=evaluateTest(tbt,curr_test,settings);
end

% Combine outcomes of boolean tests to get trial type classifications
for i=1:length(settings.trialtype)
    curr_type=settings.trialtype(i);
    out.trialtype(i).name=curr_type.name;
    out.trialtype(i).isThisType=classifyUsingBoolTests(out,curr_type,settings);
end

end

function out=classifyUsingBoolTests(bool_out,curr_type,settings)

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

end

function out=evaluateTest(tbt,curr_test,settings)

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
    case 'reach batch'
        data=tbt.(['reachBatch_' curr_test.fieldname]);  
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
    case 'first'
        % Check whether first in event set is of indicated type
        largerSet=tbt.(curr_test.inEventSet);
        out=isFirstEventThis(data(:,startInds:endInds),largerSet(:,startInds:endInds),curr_test,settings);
        out=out';
    otherwise
        % let 'any' be default mode
        out=evaluateBoolTest(data(:,startInds:endInds),curr_test,settings,2);
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

function tbt=findReachBatches(tbt,lowThresh,showExampleChanges)

% Need to transpose matrix for operations in this function
f=fieldnames(tbt);
for i=1:length(f)
    temp=tbt.(f{i});
    temp=temp';
    tbt.(f{i})=temp;
end

% Finds reach batches according to definition in trialTypeSettings.m
settings=trialTypeSettings();
reach_batch=settings.reach_batch;

% Make new fields for reach batches
for i=1:length(reach_batch.firstreach_type)
    currtype=reach_batch.firstreach_type{i};
    newstr=['reachBatch_' currtype];
    tbt.(newstr)=tbt.(currtype);
end
for i=1:length(reach_batch.secondreach_type)
    currtype=reach_batch.secondreach_type{i};
    newstr=['reachBatch_' currtype];
    tbt.(newstr)=tbt.(currtype);
end
if reach_batch.removePawOnWheel==1
    if reach_batch.take_first_or_second_type==1
        for i=1:length(reach_batch.firstreach_type)
            currtype=reach_batch.firstreach_type{i};
            if ~isempty(strfind(currtype,'pawOnWheel'))
                r=regexp(currtype,'pawOnWheel');
                newstr=['reachBatch_' currtype(1:r-2)];
                tbt.(newstr)=tbt.(currtype);
            end
        end
    elseif reach_batch.take_first_or_second_type==2
        for i=1:length(reach_batch.secondreach_type)
            currtype=reach_batch.secondreach_type{i};
            if ~isempty(strfind(currtype,'pawOnWheel'))
                r=regexp(currtype,'pawOnWheel');
                newstr=['reachBatch_' currtype(1:r-2)];
                tbt.(newstr)=tbt.(currtype(1:r-2));
            end
        end
    end
end

% Get times of various reach types
secondtype_f=cell(1,length(reach_batch.secondreach_type));
for i=1:length(reach_batch.secondreach_type)
    temp=tbt.(reach_batch.secondreach_type{i});
    tbtsecondtype=temp(1:end);
    f=find(tbtsecondtype>lowThresh);
    secondtype_f{i}=f;
end

showRows=nan(1,size(tbt.times,1));
showRows_inc=1;
for i=1:length(reach_batch.firstreach_type)
    currtype=reach_batch.firstreach_type{i};
    % Find all reaches of this type
    % For each reach of this type, check whether a second reach of an appropriate type occurs
    % within window seconds
    temp=tbt.(currtype);
    tbtcurrtype=temp(1:end);
    f=find(tbtcurrtype>lowThresh);
    for j=1:length(f)
        [row,col]=ind2sub(size(tbt.(currtype)),f(j));
        candidate_secondreaches.inds=[];
        candidate_secondreaches.types=[]; 
        for k=1:length(reach_batch.secondreach_type)
            f2=secondtype_f{k};
            ne_ind=find(f2>f(j),1,'first');
            ne=f2(ne_ind);
            % is ne in the same trial and within window?
            [row2,col2]=ind2sub(size(tbt.(currtype)),ne);
            if col==col2 % reaches occur in same trial
                % do reaches occur within time window defining reach batch?
                timediff=tbt.times(row2,col2)-tbt.times(row,col);
                if timediff<reach_batch.window
                    % candidate second reach
                    candidate_secondreaches.inds=[candidate_secondreaches.inds ne];
                    candidate_secondreaches.types=[candidate_secondreaches.types k];
                end
            end
        end
        % find reach immediately following first reach
        if ~isempty(candidate_secondreaches.inds)
            % found a reach qualifying for reach batch
            [secondreach.ind,srind]=min(candidate_secondreaches.inds);
            secondreach.type=candidate_secondreaches.types(srind);
            % save for plotting
            if showExampleChanges==1
                showRows(showRows_inc)=col; % save trial that I am changing
                showRows_inc=showRows_inc+1;
            end
            % fix reach batch 
            if reach_batch.take_first_or_second_type==1
                % take the reach type of first reach
                if reach_batch.take_first_or_second_timing==1
                    % take timing of first reach in batch
                    % thus zero out second reach in batch
                elseif reach_batch.take_first_or_second_timing==2
                    % take timing of second reach in batch
                    % thus switch time of first reach in batch
                    % and zero out second reach in batch
                    newstr=['reachBatch_' currtype];
                    if reach_batch.removePawOnWheel==1
                        if ~isempty(strfind(newstr,'pawOnWheel'))
                            r=regexp(newstr,'pawOnWheel');
                            newstr=newstr(1:r-2);
                        end
                    end
                    tempreaches=tbt.(newstr);
                    [row2,col2]=ind2sub(size(tbt.(currtype)),secondreach.ind);
                    tempreaches(row2,col2)=1;
                    tempreaches(row,col)=0;
                    tbt.(newstr)=tempreaches;
                    if reach_batch.removePawOnWheel==1
                        newstr=['reachBatch_' currtype];
                        if ~isempty(strfind(newstr,'pawOnWheel'))
                            tempreaches=tbt.(newstr);
                            tempreaches(row,col)=0;
                            tbt.(newstr)=tempreaches;
                        end
                    end
                end
                % zero out second reach in batch
                if isfield(tbt,['reachBatch_' reach_batch.secondreach_type{secondreach.type}])
                    if reach_batch.removePawOnWheel==1 && ~isempty(strfind(reach_batch.secondreach_type{secondreach.type},currtype))
                    else
                        newstr=['reachBatch_' reach_batch.secondreach_type{secondreach.type}];
                        tempreaches=tbt.(newstr);
                        [row2,col2]=ind2sub(size(tbt.(currtype)),secondreach.ind);
                        tempreaches(row2,col2)=0;
                        tbt.(newstr)=tempreaches;
                    end
                end
            elseif reach_batch.take_first_or_second_type==2
                % take the reach type of second reach
                if reach_batch.take_first_or_second_timing==1
                    % take timing of first reach in batch
                    % thus switch time of second reach in batch
                    % and zero out first reach in batch
                    newstr=['reachBatch_' reach_batch.secondreach_type{secondreach.type}];
                    if reach_batch.removePawOnWheel==1
                        if ~isempty(strfind(newstr,'pawOnWheel'))
                            r=regexp(newstr,'pawOnWheel');
                            newstr=newstr(1:r-2);
                        end
                    end
                    tempreaches=tbt.(newstr);
                    tempreaches(row,col)=1;
                    [row2,col2]=ind2sub(size(tbt.(currtype)),secondreach.ind);
                    tempreaches(row2,col2)=0;
                    tbt.(newstr)=tempreaches;
                    if reach_batch.removePawOnWheel==1
                        newstr=['reachBatch_' reach_batch.secondreach_type{secondreach.type}];
                        if ~isempty(strfind(newstr,'pawOnWheel'))
                            tempreaches=tbt.(newstr);
                            tempreaches(row2,col2)=0;
                            tbt.(newstr)=tempreaches;
                        end
                    end
                elseif reach_batch.take_first_or_second_timing==2
                    % take timing of second reach in batch
                    % thus zero out first reach in batch
                end
                % zero out first reach in batch
                if isfield(tbt,['reachBatch_' currtype])
                    if reach_batch.removePawOnWheel==1 && ~isempty(strfind(reach_batch.secondreach_type{secondreach.type},currtype))
                    else
                        newstr=['reachBatch_' currtype];
                        tempreaches=tbt.(newstr);
                        tempreaches(row,col)=0;
                        tbt.(newstr)=tempreaches;
                    end
                end                
            end
        else
            % no reaches in reach batch
        end
    end
end

% Transpose matrix again to return to original format
f=fieldnames(tbt);
for i=1:length(f)
    temp=tbt.(f{i});
    temp=temp';
    tbt.(f{i})=temp;
end

% Plot changes to ensure proper functioning
if showExampleChanges==1
    plotEventScatter(tbt,showRows,0);
    title('BEFORE fixing reach batches');
    plotEventScatter(tbt,showRows,1);
    title('AFTER fixing reach batches');
end

end

function plotEventScatter(tbt,showRows,doReachBatch)

timespertrial=nanmean(tbt.times,1);
showRows=showRows(~isnan(showRows));
plotSettings=plotCueTriggered_settings();
% Plot old version before correcting reach batches
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