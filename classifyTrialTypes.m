function classifyTrialTypes(tbt)

% User-defined settings in trialTypeSettings.m
settings=trialTypeSettings();

% Find and convert reach batches
tbt=findReachBatches(tbt,settings.lowThresh,1);

for i=1:length(settings.bool_test)
    curr_test=settings.bool_test(i);
    out.bool_test(i).testResults=nan(1,size(tbt.times,1)); % test each trial
    switch curr_test.test
        case 'any'
            useFunc='any';
        case 'first'
            useFunc='find';
            n=1;
        otherwise
            useFunc='any'; % default mode
    end
    switch curr_test.testwhat
        case 'single reach'
            temp=tbt.(fieldname);
            
            
        case 'reach batch'
            
            
            
        otherwise
    end
end

bool_test(i).testwhat='single reach';
bool_test(i).fieldname='reach_ongoing';
bool_test(i).test='any';
bool_test(i).thresh=settings.lowThresh;
bool_test(i).comparator='>';
bool_test(i).window='wheel_turning';



trialtype(1).outcomes=[0 1 0 0 0 nan nan nan nan nan];
trialtype(1).name='after_cue_success';
trialtype(1).bool='&';



% Exclude trials where paw was on wheel while wheel turning
if excludePawOnWheelTrials==1
    % Find trials where paw was on wheel while wheel turning
    plot_cues=[];
    for i=1:size(tbt.(nameOfCue),1)
        presentInd=find(tbt.pelletPresented(i,:)>0.5,1,'first');
        temp=tbt.(nameOfCue);
        cueInd=find(temp(i,:)>0.5,1,'first');
        pawWasOnWheel=0;
        if any(tbt.pawOnWheel(i,presentInd:cueInd)>0.5)
            pawWasOnWheel=1;
        else
            plot_cues=[plot_cues i];
        end
    end
else
    plot_cues=1:size(tbt.(nameOfCue),1);
end
if settings.excludeFirstTrial==1
    plot_cues=plot_cues(~ismember(plot_cues,1));
end




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
    case 'firstReach'
        
    otherwise
        % let 'any' be default mode
        out=evaluateBoolTest(data(:,startInds:endInds),curr_test,settings,2);
end

end

function out=evaluateAnyTest(data,curr_test,settings,dim)

out=eval([useFunc '(data' curr_test.comparator num2str(curr_test.thresh) ',' num2str(dim) ')']);

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