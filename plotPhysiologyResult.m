function [fout,dataout,n_events_in_av,alignmentCompanion,f_heatmap,plotBehFieldOut,su]=plotPhysiologyResult(photometry_tbt,behavior_tbt,outdata,alignTo,plotPhotoField,plotBehField,firstInTrialOrLast,withinRange,hax)

% alignTo can be 'cue','success_reachStarts', etc.
fout=[];
dataout=[];
n_events_in_av=[];
alignmentCompanion=[];
f_heatmap=[];
plotBehFieldOut=[];
su=[];

if strcmp(plotPhotoField,'unit_by_unit')
    doingMultipleUnits=true;
    nonUnitFields={'unitTimes','unitsum','multiunit','sum_over_singleunit','av_over_singleunit','unitTimes_wrt_trial_start','grp1_unitav','grp2_unitav','grp3_unitav','grpDA1_unitav','grpDA2_unitav','grpDA3_unitav','grpnegDA1_unitav','grpnegDA2_unitav','grpnegDA3_unitav','grpconsensusDA1_unitav','grpconsensusDA2_unitav','grpconsensusDA3_unitav'}; % every other field that contains 'unit' is assumed to be one single unit
    plotPhotoField='unitsum';
else
    doingMultipleUnits=false;
end

if strcmp(firstInTrialOrLast,'first')
    firstInTrial=true;
    lastInTrial=false;
elseif strcmp(firstInTrialOrLast,'last')
    firstInTrial=false;
    lastInTrial=true;
else
    firstInTrial=false;
    lastInTrial=false;
end

if ~isfield(behavior_tbt,'all_reachBatch')
    disp('behavior_tbt is missing all_reachBatch ... WILL SUBSTITUTE WITH REACHSTARTS!');
    behavior_tbt.all_reachBatch=behavior_tbt.reachStarts;
    behavior_tbt.reachBatch_success_reachStarts=behavior_tbt.success_reachStarts;
    behavior_tbt.reachBatch_drop_reachStarts=behavior_tbt.drop_reachStarts;
    behavior_tbt.reachBatch_miss_reachStarts=behavior_tbt.miss_reachStarts;
end
    
cutBeforeNextCue=true; % if is true, will only plot inter-trial interval after first cue
alwaysRealignToCue=false;
minITI=9; % in seconds
behaviorCue='cueZone_onVoff';
dosmooth=true;
% lowPassCutoff=5; % in Hz
lowPassCutoff=100000; % in Hz
alignPeaks=false;
indsFromPeak=5;
% withinRange=[-0.1 1.5];  % take reaches in this range, time is in seconds from onset of cue
% withinRange=[1.5 16];  % take reaches in this range, time is in seconds from onset of cue
% withinRange=[-0.1 1.5];  % take reaches in this range, time is in seconds from onset of cue
% withinRange=[-0.1 16];  % take reaches in this range, time is in seconds from onset of cue
downSamp=true;
ds=100;
maxTimeForHeatmap=9;

if downSamp==true
    downSampWhichFields={'cue','opto','cue_times','distractor','cuetimes_wrt_trial_start'};
    discardWhichFields={'phys_timepoints','phys_inds','from_first_video','from_second_video','from_third_video'};
    f=fieldnames(photometry_tbt);
    for i=1:length(f)
        if ismember(f{i},discardWhichFields)
            rmfield(photometry_tbt,f{i});
            continue
        end
        if ismember(f{i},downSampWhichFields)
            photometry_tbt.(f{i})=downSampMatrix(photometry_tbt.(f{i}),ds);
        end
    end
end

thesePhotoFieldsUseTimeField1={'cue'};
timeField1='cue_times';
f=fieldnames(photometry_tbt);
thesePhotoFieldsUseTimeField2={};
for i=1:length(f)
    if ~isempty(regexp(f{i},'unit'))
        if strcmp(f{i},'unitTimes')
            continue
        end
        thesePhotoFieldsUseTimeField2{length(thesePhotoFieldsUseTimeField2)+1}=f{i};
    end
end
timeField2='unitTimes';

temp=photometry_tbt.(timeField1);
photometry_tbt.([timeField1 '_wrt_trial_start'])=temp-repmat(temp(:,1),1,size(temp,2));
temp=photometry_tbt.(timeField2);
photometry_tbt.([timeField2 '_wrt_trial_start'])=temp-repmat(temp(:,1),1,size(temp,2));

if strcmp(alignTo,'cue') || alwaysRealignToCue==true
    [photometry_tbt,alignedCueTo]=realignToCue(photometry_tbt,'cue',thesePhotoFieldsUseTimeField1,timeField1,timeField2);
end
[~,fpe]=nanmax(nanmean(behavior_tbt.cueZone_onVoff,1));
temp=nanmean(behavior_tbt.times_wrt_trial_start,1);
beh_cue=temp(fpe);
ma=nanmax(nanmean(photometry_tbt.cue,1));
fpe=find(nanmean(photometry_tbt.cue,1)>ma*0.75,1,'first');
temp=nanmean(photometry_tbt.cue_times_wrt_trial_start,1);
photo_cue=temp(fpe);
f=fieldnames(photometry_tbt);
for i=1:length(f)
    if ~isempty(regexp(f{i},'time'))
        photometry_tbt.(f{i})=photometry_tbt.(f{i})+(beh_cue-photo_cue);
    end
end
figure();
plot(nanmean(photometry_tbt.cue_times_wrt_trial_start,1),nanmean(photometry_tbt.cue,1),'Color','k');
hold on;
plot(nanmean(behavior_tbt.times_wrt_trial_start,1),nanmean(behavior_tbt.cueZone_onVoff,1),'Color','r');
title('realigned to cue');

if ismember(plotPhotoField,thesePhotoFieldsUseTimeField1)
    getCorrectTime=timeField1;
    getCorrectTime_wrtTrialStart=[timeField1 '_wrt_trial_start'];
elseif ismember(plotPhotoField,thesePhotoFieldsUseTimeField2)
    getCorrectTime=timeField2;
    getCorrectTime_wrtTrialStart=[timeField2 '_wrt_trial_start'];
else
    error('Do not recognize this photometry_tbt field name.');
end
temp=photometry_tbt.(getCorrectTime);
% find first non-nan column
f=find(~isnan(temp(1,:)),1,'first');
phototimes=nanmean(photometry_tbt.(getCorrectTime)-repmat(temp(:,f),1,size(temp,2)),1);

if cutBeforeNextCue==true
    % find second cue time
    temp=mode(photometry_tbt.cue_times(2,2:end)-photometry_tbt.cue_times(2,1:end-1));
    if temp==0
        temp=mode(photometry_tbt.cue_times(3,2:end)-photometry_tbt.cue_times(3,1:end-1));
    end
    firstSearchInd=floor(minITI./temp);
    f=find(nanmean(photometry_tbt.cue(:,firstSearchInd:end),1)>0.01,1,'first');
    f=f+firstSearchInd;
    
    % find first non-nan column
    ftostart=find(~isnan(photometry_tbt.cue_times(1,:)),1,'first');
    temp=photometry_tbt.(getCorrectTime);
    % find first non-nan column
    ftostart_photo=find(~isnan(temp(1,:)),1,'first');
    
    f_times=nanmean(photometry_tbt.cue_times(:,f)-photometry_tbt.cue_times(:,ftostart),1);
    thisTimeField=photometry_tbt.(timeField2);
    [~,f_photometry_ind]=nanmin(abs(nanmean(thisTimeField-repmat(thisTimeField(:,ftostart_photo),1,size(thisTimeField,2)),1)-f_times));
    if ismember(plotPhotoField,thesePhotoFieldsUseTimeField1)
        plotUntilInd=f;
    elseif ismember(plotPhotoField,thesePhotoFieldsUseTimeField2)
        plotUntilInd=f_photometry_ind;
    else
        error('Do not recognize this photometry_tbt field name.');
    end
    % fix phototimes
    if ismember(plotPhotoField,thesePhotoFieldsUseTimeField1)
        temp=phototimes(2:f)-phototimes(1:f-1);
        temp=temp(~isnan(temp));
        ts=mode(temp);
    elseif ismember(plotPhotoField,thesePhotoFieldsUseTimeField2)
        temp=phototimes(2:f_photometry_ind)-phototimes(1:f_photometry_ind-1);
        temp=temp(~isnan(temp));
        ts=mode(temp);
    else
        error('Do not recognize this photometry_tbt field name.');
    end
    phototimes=0:ts:(length(phototimes)-1)*ts;
else
    plotUntilInd=[];
end

if dosmooth==true
    photometry_tbt=smoothPhotometry(photometry_tbt,1/mode(diff(phototimes)),lowPassCutoff,plotPhotoField);
    disp(['using photometry Fs ' num2str(1/mode(diff(phototimes)))]);
end

typeOfReach=false;
useCombo=false;
switch alignTo
    case 'cue'
        ftostart=find(~isnan(photometry_tbt.cue_times(1,:)),1,'first');
        f_times=nanmean(photometry_tbt.cue_times-repmat(photometry_tbt.cue_times(:,ftostart),1,size(photometry_tbt.cue_times,2)),1);
        if ~isempty(hax)
            if ~isnumeric(hax{1})
                axes(hax{1});
                f_heatmap=nan;
            else
                f_heatmap=figure();
            end
        else
            f_heatmap=figure();
        end
        plotTrialsAsHeatmap(photometry_tbt.(plotPhotoField),phototimes,photometry_tbt.cue,f_times,10,maxTimeForHeatmap,false);
        hold on;
        plotEventsScatter(behavior_tbt,1:size(behavior_tbt.cue,1),'cueZone_onVoff','all_reachBatch','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts',[]);
        if ~isempty(hax)
            if ~isnumeric(hax{2})
                axes(hax{2});
                fout=nan;
            else
                fout=figure();
            end
        else
            fout=figure();
        end
        plotWStderr(photometry_tbt.(plotPhotoField),phototimes,'k',plotUntilInd,size(photometry_tbt.cue,1));
        dataout.x=phototimes;
        dataout.y=photometry_tbt.(plotPhotoField);
        hold on;
        temp=photometry_tbt.(plotPhotoField);
        plotWStderr(nanmean(photometry_tbt.cue,1)*(nanmax(nanmean(temp(:,1:plotUntilInd),1))/nanmax(nanmean(photometry_tbt.cue,1))),f_times,'b',f,size(photometry_tbt.cue,1));
        n_events_in_av=size(photometry_tbt.(plotPhotoField),1);
        alignmentCompanion.x=f_times;
        alignmentCompanion.y=photometry_tbt.cue.*nanmax(nanmean(photometry_tbt.(plotPhotoField),1));
        if ~isempty(plotBehField)
            plotBehFieldOut.x=nanmean(behavior_tbt.times_wrt_trial_start,1);
            plotBehFieldOut.y=behavior_tbt.(plotBehField);
        end
        if doingMultipleUnits==true
            tryforSUfields=fieldnames(photometry_tbt);
            j=1;
            for i=1:length(tryforSUfields)
                if ~isempty(regexp(tryforSUfields{i},'unit','ONCE')) && ~ismember(tryforSUfields{i},nonUnitFields)
                    % this is a single unit
                    SUfields{j}=tryforSUfields{i};
                    j=j+1;
                end
            end
            disp('will plot units in this order:');
            for i=1:length(SUfields)
                su(i).alignedData=photometry_tbt.(SUfields{i});
                disp(SUfields{i});
            end
        end
    case 'success_fromPerchOrWheel'
        typeOfReach=true;
        useReach='combo';
        useCombo=behavior_tbt.reachBatch_success_reachStarts+behavior_tbt.reachBatch_success_reachStarts_pawOnWheel;
    case 'drop_fromPerchOrWheel'
        typeOfReach=true;
        useReach='combo';
        useCombo=behavior_tbt.reachBatch_drop_reachStarts+behavior_tbt.reachBatch_drop_reachStarts_pawOnWheel;
    case 'miss_fromPerchOrWheel'
        typeOfReach=true;
        useReach='combo';
        useCombo=behavior_tbt.reachBatch_miss_reachStarts+behavior_tbt.reachBatch_miss_reachStarts_pawOnWheel;
    case 'misses_and_pelletMissing'
        typeOfReach=true;
        useReach='combo';
        excludeAfterSuccess=true;
        excludeWithinTimeWindow=6; % in sec
        excludeWithinInds=floor(excludeWithinTimeWindow/mode(diff(nanmean(behavior_tbt.times,1))));
        useCombo=behavior_tbt.reachBatch_miss_reachStarts+behavior_tbt.reachBatch_miss_reachStarts_pawOnWheel+behavior_tbt.pelletmissingreach_reachStarts;
        if excludeAfterSuccess==true
            for i=1:size(useCombo,1)
                f=find(useCombo(i,:)>0.5,1,'first');
                if isempty(f)
                    continue
                end
                for j=1:length(f)
                    currf=f(j);
                    starter=currf-excludeWithinInds;
                    if starter<1
                        starter=1;
                    end
                    if any(behavior_tbt.success_reachStarts(i,starter:currf)>0.5 | behavior_tbt.reachBatch_success_reachStarts(i,starter:currf)>0.5 | behavior_tbt.reachBatch_success_reachStarts_pawOnWheel(i,starter:currf)>0.5)
                        % success before this, exclude
                        useCombo(i,f)=0;
                    end
                end
            end
        end
    case 'misses_and_pelletMissing_and_drop'
        typeOfReach=true;
        useReach='combo';
        excludeAfterSuccess=true;
        excludeWithinTimeWindow=6; % in sec
        excludeWithinInds=floor(excludeWithinTimeWindow/mode(diff(nanmean(behavior_tbt.times,1))));
        useCombo=behavior_tbt.reachBatch_miss_reachStarts+behavior_tbt.reachBatch_miss_reachStarts_pawOnWheel+behavior_tbt.pelletmissingreach_reachStarts+behavior_tbt.reachBatch_drop_reachStarts+behavior_tbt.reachBatch_drop_reachStarts_pawOnWheel;
        if excludeAfterSuccess==true
            for i=1:size(useCombo,1)
                f=find(useCombo(i,:)>0.05,1,'first');
                if isempty(f)
                    continue
                end
                for j=1:length(f)
                    currf=f(j);
                    starter=currf-excludeWithinInds;
                    if starter<1
                        starter=1;
                    end
                    if any(behavior_tbt.success_reachStarts(i,starter:currf)>0.5 | behavior_tbt.reachBatch_success_reachStarts(i,starter:currf)>0.5 | behavior_tbt.reachBatch_success_reachStarts_pawOnWheel(i,starter:currf)>0.5)
                        % success before this, exclude
                        upto=currf+10; 
                        if upto>length(useCombo(i,:))
                            upto=length(useCombo(i,:));
                        end
                        useCombo(i,starter:upto)=0;
                    end
                end
            end
        end
    case 'reachBatch_success_reachStarts'
        typeOfReach=true;
        useReach=alignTo;
        behavior_tbt=findPelletDislodgedAfterSuccess(behavior_tbt,'reachBatch_success_reachStarts','pelletPresent');
    case 'reachBatch_drop_reachStarts'
        typeOfReach=true;
        useReach=alignTo;
        behavior_tbt=findPelletDislodgedAfterSuccess(behavior_tbt,'reachBatch_drop_reachStarts','pelletPresent');
    case 'success batch when pellet dislodged'
        typeOfReach=true;
        useReach='combo';
        behavior_tbt=findPelletDislodgedAfterSuccess(behavior_tbt,'reachBatch_success_reachStarts','pelletPresent');
        useCombo=behavior_tbt.pelletDislodgedAfterSuccess;
    case 'success batch av reach and pellet dislodged'
        typeOfReach=true;
        useReach='combo';
        behavior_tbt=findPelletDislodgedAvWithSuccess(behavior_tbt,'reachBatch_success_reachStarts','pelletPresent');
        useCombo=behavior_tbt.pelletDislodgedAfterSuccess;
    case 'drop batch when pellet dislodged'
        typeOfReach=true;
        useReach='combo';
        behavior_tbt=findPelletDislodgedAfterSuccess(behavior_tbt,'reachBatch_drop_reachStarts','pelletPresent');
        useCombo=behavior_tbt.pelletDislodgedAfterSuccess;
    case 'drop batch av reach and pellet dislodged'
        typeOfReach=true;
        useReach='combo';
        behavior_tbt=findPelletDislodgedAvWithSuccess(behavior_tbt,'reachBatch_drop_reachStarts','pelletPresent');
        useCombo=behavior_tbt.pelletDislodgedAfterSuccess;
    otherwise 
        typeOfReach=true;
        useReach=alignTo;
end
if typeOfReach==true
    % in behavior_tbt
    if strcmp(useReach,'combo')
        temp=useCombo;
    else
        temp=behavior_tbt.(useReach);
    end
    % note that pellet is only present after cue
    nIndsToTake=200;
    timeStep=mode(diff(nanmean(behavior_tbt.times,1)));
    cueInd=find(nanmean(behavior_tbt.(behaviorCue),1)>0,1,'first');
    withinRange_inds=[cueInd+ceil(withinRange(1)/timeStep) cueInd+ceil(withinRange(2)/timeStep)];
    if firstInTrial==true
        [behAligned,alignInds,fromInputRow]=alignToEventFirstInTrial(temp,withinRange_inds,temp,nIndsToTake,90);
    elseif lastInTrial==true
        [behAligned,alignInds,fromInputRow]=alignToEventLastInTrial(temp,withinRange_inds,temp,nIndsToTake,90);
    else
        [behAligned,alignInds,fromInputRow]=alignToEvent(temp,withinRange_inds,temp,nIndsToTake,90);
    end
    alignTimes=getValsFromRows(behavior_tbt.times_wrt_trial_start,alignInds,fromInputRow);
    indsIntoPhoto=getIndsFromRows(photometry_tbt.(getCorrectTime_wrtTrialStart),alignTimes,fromInputRow);
    [alignedData,alignedAt]=alignRowsToInds(photometry_tbt.(plotPhotoField),indsIntoPhoto,nanmin(indsIntoPhoto),fromInputRow,alignPeaks,indsFromPeak);
    if doingMultipleUnits==true
        tryforSUfields=fieldnames(photometry_tbt);
        j=1;
        for i=1:length(tryforSUfields)
            if ~isempty(regexp(tryforSUfields{i},'unit','ONCE')) && ~ismember(tryforSUfields{i},nonUnitFields)
                % this is a single unit
                SUfields{j}=tryforSUfields{i};
                j=j+1;
            end
        end
        disp('will plot units in this order:');
        for i=1:length(SUfields)
            [su(i).alignedData,~]=alignRowsToInds(photometry_tbt.(SUfields{i}),indsIntoPhoto,nanmin(indsIntoPhoto),fromInputRow,alignPeaks,indsFromPeak);
            disp(SUfields{i});
        end
    end
    
    % check other reaches with respect to this reach type
    [allReachesAlignedData,allReachesAlignedAt]=alignRowsToInds(behavior_tbt.reachStarts,alignInds,300,fromInputRow,false,indsFromPeak);
    behTimesAllReaches=nanmean(behavior_tbt.times_wrt_trial_start(:,1:size(allReachesAlignedData,2)),1);
    behTimes=nanmean(behavior_tbt.times_wrt_trial_start(:,1:size(behAligned,2)),1);
    [~,bmax]=nanmax(nanmean(behAligned,1));
    figure();
    plotWStderr(behAligned,behTimes,'g',[],size(behAligned,1));
    hold on;
    plotWStderr(allReachesAlignedData,behTimesAllReaches-(behTimesAllReaches(allReachesAlignedAt)-behTimes(bmax)),'k',[],size(allReachesAlignedData,1));
    
    if ~isempty(plotBehField)
        % plot another behavior event aligned in same way 
        [allReachesAlignedData,allReachesAlignedAt]=alignRowsToInds(behavior_tbt.(plotBehField),alignInds,300,fromInputRow,false,indsFromPeak);
        figure();
        plotWStderr(behAligned,behTimes,'g',[],size(behAligned,1));
        hold on;
        plotWStderr(allReachesAlignedData,behTimesAllReaches-(behTimesAllReaches(allReachesAlignedAt)-behTimes(bmax)),'k',[],size(allReachesAlignedData,1));
        title([plotBehField ' (black) aligned in same way']);
        plotBehFieldOut.x=behTimesAllReaches-(behTimesAllReaches(allReachesAlignedAt)-behTimes(bmax));
        plotBehFieldOut.y=allReachesAlignedData;
    end
    
    figure();
    plotTrialsAsHeatmap(alignedData,phototimes(1:size(alignedData,2))-(phototimes(alignedAt)-behTimes(bmax)),behAligned,behTimes,10,maxTimeForHeatmap,true);
    
    if ~isempty(hax)
        if ~isnumeric(hax{1})
            axes(hax{1});
            f_heatmap=nan;
        else
            f_heatmap=figure();
        end
    else
        f_heatmap=figure();
    end
    plotTrialsAsHeatmap(alignedData,phototimes(1:size(alignedData,2))-(phototimes(alignedAt)-behTimes(bmax)),behAligned,behTimes,10,maxTimeForHeatmap,true);
    hold on;
    whichFields={'cueZone_onVoff','all_reachBatch','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts','success_reachStarts_pawOnWheel'};
    if ~isempty(plotBehField)
        whichFields{length(whichFields)+1}=plotBehField;
    end
    if isfield(behavior_tbt,'pelletDislodgedAfterSuccess')
        whichFields{length(whichFields)+1}='pelletDislodgedAfterSuccess';
    end
    [newBehavior_tbt,allbalignedAt]=makeShiftedTbt(behavior_tbt,whichFields,alignInds,300,fromInputRow);
    newBehavior_tbt.times=behTimesAllReaches-(behTimesAllReaches(allbalignedAt)-behTimes(bmax));
    plotEventsScatter(newBehavior_tbt,fromInputRow(~isnan(fromInputRow)),'cueZone_onVoff','all_reachBatch','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts','success_reachStarts_pawOnWheel');
    
    if ~isempty(plotBehField)
        figure();
        plotTrialsAsHeatmap(alignedData,phototimes(1:size(alignedData,2))-(phototimes(alignedAt)-behTimes(bmax)),behAligned,behTimes,10,maxTimeForHeatmap,true);
        hold on;
        plotEventsScatter(newBehavior_tbt,fromInputRow(~isnan(fromInputRow)),'cueZone_onVoff','all_reachBatch','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts',plotBehField);
        title(['Black dots are ' plotBehField]);
    end
    
    figure();
    hold on;
    plotEventsScatter(behavior_tbt,fromInputRow(~isnan(fromInputRow)),'cueZone_onVoff','all_reachBatch','reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts','success_reachStarts_pawOnWheel');
    if ~isempty(hax)
        if ~isnumeric(hax{2})
            axes(hax{2});
            fout=nan;
        else
            fout=figure();
        end
    else
        fout=figure();
    end
    plotWStderr(alignedData,phototimes(1:size(alignedData,2))-(phototimes(alignedAt)-behTimes(bmax)),'k',[],size(behAligned,1));
    dataout.x=phototimes(1:size(alignedData,2))-(phototimes(alignedAt)-behTimes(bmax));
    dataout.y=alignedData;
    if doingMultipleUnits==true
        dataout.y=su;
    end
    hold on;
    te=nanmean(alignedData,1);
    te=te(~isnan(te));
    if nansum(te<=0)<length(te)*0.9
        % scale behavior
        plotWStderr((behAligned./nanmax(nanmean(behAligned,1))).*range(te(te>0))+nanmin(te(te>0)),behTimes,'g',[],size(behAligned,1));
        alignmentCompanion.x=behTimes;
        alignmentCompanion.y=(behAligned./nanmax(nanmean(behAligned,1))).*range(te(te>0))+nanmin(te(te>0));
    else
        plotWStderr(behAligned,behTimes,'g',[],size(behAligned,1));
        alignmentCompanion.x=behTimes;
        alignmentCompanion.y=behAligned;
    end
    disp([num2str(size(behAligned,1)) ' events averaged']);
    n_events_in_av=size(behAligned,1);
end

end

function behavior_tbt=findPelletDislodgedAvWithSuccess(behavior_tbt,successField,pelletPresent)

successes=behavior_tbt.(successField);
pellet=behavior_tbt.(pelletPresent);
behavior_tbt.pelletDislodgedAfterSuccess=nan(size(successes));
for i=1:size(successes,1)
    f=find(successes(i,:)>0.5);
    pelletGoneInds=nan(1,length(f));
    currpellet=pellet(i,:);
    for j=1:length(f)
        % find first index when pellet dislodged after success
        temp=f(j)-1+find(currpellet(f(j):end)<0.5,1,'first');
        if isempty(temp)
            continue
        end
        pelletGoneInds(j)=round(nanmean([temp f(j)],2));
    end
    behavior_tbt.pelletDislodgedAfterSuccess(i,:)=zeros(size(behavior_tbt.pelletDislodgedAfterSuccess(i,:)));
    behavior_tbt.pelletDislodgedAfterSuccess(i,pelletGoneInds)=1;
end

end

function behavior_tbt=findPelletDislodgedAfterSuccess(behavior_tbt,successField,pelletPresent)

successes=behavior_tbt.(successField);
pellet=behavior_tbt.(pelletPresent);
behavior_tbt.pelletDislodgedAfterSuccess=nan(size(successes));
for i=1:size(successes,1)
    f=find(successes(i,:)>0.5);
    pelletGoneInds=nan(1,length(f));
    currpellet=pellet(i,:);
    for j=1:length(f)
        % find first index when pellet dislodged after success
        temp=f(j)-1+find(currpellet(f(j):end)<0.5,1,'first');
        if isempty(temp)
            continue
        end
        pelletGoneInds(j)=temp;
    end
    behavior_tbt.pelletDislodgedAfterSuccess(i,:)=zeros(size(behavior_tbt.pelletDislodgedAfterSuccess(i,:)));
    if isnan(pelletGoneInds)
        behavior_tbt.pelletDislodgedAfterSuccess(i,end)=1;
    else
        behavior_tbt.pelletDislodgedAfterSuccess(i,pelletGoneInds)=1;
    end
end

end

function data=smoothPhotometry(data,Fs,lowPassCutoff,currField)

disp('Smoothing fields');
temp=data.(currField);
% truncate after nans begin AFTER real signal ends
for j=1:size(temp,1)
    frealsig=find(~isnan(temp(j,:)) & temp(j,:)~=0,1,'last');
    fna=find(isnan(temp(j,frealsig+1:end)),1,'first');
    if ~isempty(fna)
        fna=frealsig+fna;
        temp(j,fna:end)=nan;
    end
end
data.(currField)=temp;

disp(['Smoothing field named ' currField]);
% can't put in nans, pad with a low value
padval=prctile(temp(1:end),10);
temp(isnan(temp))=padval;
newtemp=fftFilter(temp',Fs,lowPassCutoff,1);
beforenans=real(newtemp');
% put nans back in
beforenans(isnan(data.(currField)))=nan;
data.(currField)=beforenans;
    
end

function [newBehavior_tbt,alignedAt]=makeShiftedTbt(behavior_tbt,whichFields,alignInds,nBefore,fromInputRow)

for i=1:length(whichFields)
    [newBehavior_tbt.(whichFields{i}),alignedAt]=alignRowsToInds(behavior_tbt.(whichFields{i}),alignInds,nBefore,fromInputRow,false,[]);
end

end

function plotEventsScatter(alltbt,plotThese,cueName,reachName,successName,dropName,missName,pelletMissingName,otherName)

plotReachesWithinSec=16; % plot only first reach if isempty, else plot all reaches within this many secs of first reach
plotAllEvents=true;
thresh=0.05;
cueInd=find(nanmean(alltbt.(cueName),1)>thresh,1,'first');
if isvector(alltbt.times)
    useTimes=alltbt.times;
else
    if length(size(alltbt.times))>1
        timeStep=mode(diff(nanmean(alltbt.times,1)));
    else
        timeStep=mode(diff(alltbt.times));
    end
    useTimes=0:timeStep:(size(alltbt.times,2)-1)*timeStep;
end
reaches=alltbt.(reachName);
successes=alltbt.(successName);
drops=alltbt.(dropName);
misses=alltbt.(missName);
pelletMissing=alltbt.(pelletMissingName);
if ~isempty(otherName)
    others=alltbt.(otherName);
else
    others=[];
end

xlabel('Time (sec)');

if length(plotThese)==size(alltbt.(cueName),1)
    plotThese=1:length(plotThese);
end

yi=0;
for incr=1:length(plotThese)
    i=plotThese(incr);
    yi=yi+1;
%     if isfield(alltbt,'optoOn')
%         startOptoOn=find(alltbt.optoOn(i,:)>thresh,1,'first');
%         endOptoOn=find(alltbt.optoOn(i,startOptoOn:end)<thresh,1,'first');
%         if ~isempty(startOptoOn)
%             line([useTimes(startOptoOn) useTimes(startOptoOn+endOptoOn)],[yi yi],'Color',[1 0.5 0.5],'LineWidth',1);
%         end
%     end
    if isfield(alltbt,'pelletDislodgedAfterSuccess')
        fdis=find(alltbt.pelletDislodgedAfterSuccess(i,:)>thresh,1,'first');
        scatter(useTimes(fdis),yi,[],[0.8 0.8 0.8],'fill');
    end
    if length(size(alltbt.times))>2
        scatter(useTimes(cueInd),yi,[],'b','fill');
    else
        temp=alltbt.(cueName);
        cueInd=find(temp(i,:)>thresh,1,'first');
        if isempty(cueInd)
        else
            scatter(useTimes(cueInd),yi,[],'b','fill');
        end
    end
    temp=reaches(i,:);
    firstReachInd=find(temp(cueInd:end)>thresh,1,'first');
    firstReachInd=cueInd+firstReachInd-1;
    if plotAllEvents==true
        allReachesInds=find(temp>thresh);
    else
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
    end
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
                scatter(useTimes(firstReachInd),yi,[],'g','fill');
            case 2
                scatter(useTimes(firstReachInd),yi,[],'r','fill');
            case 3
                scatter(useTimes(firstReachInd),yi,[],'c','fill');
            case 4
                scatter(useTimes(firstReachInd),yi,[],[0.8 0.8 0.8]);
        end
    end
    if ~isempty(others)
        f=find(others(i,:)>thresh);
        scatter(useTimes(f),yi*ones(size(f)),2,'k');
    end
end

end

function plotTrialsAsHeatmap(alignedData,times,behAligned,behTimes,ceilAbove,maxTimeInSec,plotBehLines)

dormoutlier=true;

if ~isempty(ceilAbove)
    alignedData(alignedData>ceilAbove)=ceilAbove;
end

useTheseTrials=~(nansum(isnan(alignedData),2)>0.9*size(alignedData,2));
ft=find(useTheseTrials);
for i=1:length(ft)
    temp=alignedData(ft(i),1:end-floor(size(alignedData,2)/2));
    if any(isnan(temp))
        temp=temp(~isnan(temp));
    end
%     alignedData(ft(i),1:length(temp))=smooth(temp,60);
    alignedData(ft(i),1:length(temp))=smooth(temp,50);
    alignedData(ft(i),length(temp)+1:end)=nan;
%     alignedData(ft(i),1:end-floor(size(alignedData,2)/2))=smooth(temp,60);
%     alignedData(ft(i),end-floor(size(alignedData,2)/2)-60:end)=nan;
end
[~,mi]=nanmin(abs(times-maxTimeInSec));

if dormoutlier==true
    temp=alignedData(useTheseTrials,1:mi);
    bottomcut=prctile(temp(1:end),1);
    uppercut=prctile(temp(1:end),99);
    temp(temp<bottomcut)=bottomcut;
    temp(temp>uppercut)=uppercut;
    alignedData(useTheseTrials,1:mi)=temp;
end    
imAlpha=ones(size(alignedData(useTheseTrials,1:mi)));
imAlpha(isnan(alignedData(useTheseTrials,1:mi)))=0;
% colormap('jet');
imagesc(times(1:mi),1:nansum(useTheseTrials),alignedData(useTheseTrials,1:mi),'AlphaData',imAlpha);
f=find(useTheseTrials);
if plotBehLines==true
    hold on;
    for i=1:length(f)
        [~,fp]=nanmax(behAligned(i,:));
        line([behTimes(fp) behTimes(fp)],[i-0.5 i+0.5],'Color','w','LineWidth',1);
        fp=find(behAligned(i,:)>0.5);
        for j=1:length(fp)
            line([behTimes(fp(j)) behTimes(fp(j))],[i-0.5 i+0.5],'Color','w');
        end
    end
end

end

function [data,alignedCueTo]=realignToCue(data,cueField,fieldsLikeCue,cueTimesName,otherTimesName)

[ma,cueInd]=nanmax(nanmean(data.(cueField),1)); % align all to this mode
cueHalfMax=ma/2;
temp=data.(cueTimesName);
tim=nanmean(temp-repmat(temp(:,1),1,size(temp,2)),1);
if 100>length(tim)
    ts_cue=mode(tim(2:end)-tim(1:end-1)); % cue times
else
    ts_cue=mode(tim(2:100)-tim(1:100-1)); % cue times
end
alignedCueTo=ts_cue*cueInd;
temp=data.(otherTimesName);
tim=nanmean(temp-repmat(temp(:,1),1,size(temp,2)),1);
% throw out after stops being monotonic
f=find(diff(tim)<=0,1,'first');
if isempty(f)
    f=length(tim);
end
tim=tim(1:f);
ts_other=mode(tim(2:end)-tim(1:end-1)); % other fields times

tempdata=data.(cueField);
cueLikeShiftInds=nan(1,size(tempdata,1));
otherLikeShiftInds=nan(1,size(tempdata,1));
for i=1:size(tempdata,1)
    f=find(tempdata(i,:)>cueHalfMax/4,1,'first');
    if isempty(f)
        continue
    end
    % shift to align f with cueInd
    cueLikeShiftInds(i)=cueInd-f; % if positive, pad with nan at front
    otherLikeShiftInds(i)=round((cueLikeShiftInds(i)*ts_cue)./ts_other);
end

% shift all fields accordingly
f=fieldnames(data);
for i=1:length(f)
    if size(data.(f{i}),2)==1
        continue
    end
    temp=data.(f{i});
    if ismember(f{i},fieldsLikeCue)
        for j=1:length(cueLikeShiftInds)
            if isnan(cueLikeShiftInds(j))
                continue
            end
            if cueLikeShiftInds(j)<0
                temprow=temp(j,:);
                temprow=[temprow(abs(cueLikeShiftInds(j))+1:end) nan(1,abs(cueLikeShiftInds(j)))];
                temp(j,:)=temprow;
            elseif cueLikeShiftInds(j)>0
                temprow=temp(j,:);
                temprow=[nan(1,cueLikeShiftInds(j)) temprow(1:end-cueLikeShiftInds(j))];
                temp(j,:)=temprow;
            end
        end
        data.(f{i})=temp;
    else
        for j=1:length(otherLikeShiftInds)
            if isnan(otherLikeShiftInds(j))
                continue
            end
            if otherLikeShiftInds(j)<0
                temprow=temp(j,:);
                temprow=[temprow(abs(otherLikeShiftInds(j))+1:end) nan(1,abs(otherLikeShiftInds(j)))];
                temp(j,:)=temprow;
            elseif otherLikeShiftInds(j)>0
                temprow=temp(j,:);
                temprow=[nan(1,otherLikeShiftInds(j)) temprow(1:end-otherLikeShiftInds(j))];
                temp(j,:)=temprow;
            end
        end
        data.(f{i})=temp;
    end
end

end

function [alignedData,alignedAt]=alignRowsToInds(data,alignTo,indsBeforeAlign,fromInputRow,alignPeaks,withinInds)

mi=nanmin(alignTo);
if indsBeforeAlign>=mi
    indsBeforeAlign=mi-1;
end
alignedData=nan(length(alignTo),indsBeforeAlign+size(data,2)-mi);
for i=1:length(alignTo)
    if isnan(alignTo(i))
        continue
    end
    if alignPeaks==true
        peaksAreIn=data(fromInputRow(i),:);
        [~,newf]=nanmax(peaksAreIn(alignTo(i)-withinInds:alignTo(i)+withinInds));
        % find when first gets to half-max of peak
        %newf=find(peaksAreIn(alignTo(i)-withinInds:alignTo(i)+withinInds)>=pe/2,1,'first');
        alignTo(i)=alignTo(i)-1+newf;
    end
    alignedData(i,1:length(data(fromInputRow(i),alignTo(i)-indsBeforeAlign:end)))=data(fromInputRow(i),alignTo(i)-indsBeforeAlign:end);
end
alignedAt=indsBeforeAlign+1;

end

function out=getIndsFromRows(dataVals,valsInData,fromInputRow)

out=nan(1,length(valsInData));
for i=1:length(valsInData)
    if ~isnan(valsInData(i))
        [~,out(i)]=nanmin(abs(dataVals(fromInputRow(i),:)-valsInData(i))); 
    end
end 

end

function out=getValsFromRows(dataVals,indsIntoRows,fromInputRow)

out=nan(1,length(indsIntoRows));
for i=1:length(indsIntoRows)
    if ~isnan(indsIntoRows(i))
        out(i)=dataVals(fromInputRow(i),indsIntoRows(i));
    end
end

end

function [alignedOut,eventInds,fromInputRow]=alignToEvent(event,withinRange_inds,whatToAlign,howMuchToTake,indsBeforeEvent)

% for each trial, get first event in this range
% align to this event

if withinRange_inds(1)<1
    withinRange_inds(1)=1;
end
if withinRange_inds(2)>size(event,2)
    withinRange_inds(2)=size(event,2);
end
event=event(:,withinRange_inds(1):withinRange_inds(2));

eventInds=[];
fromInputRow=[];
for i=1:size(event,1)
    if any(event(i,:)>0.5,2)
        f=find(event(i,:)>0.5);
        eventInds=[eventInds f];
        fromInputRow=[fromInputRow ones(1,length(f))*i];
    end
end
eventInds=eventInds+withinRange_inds(1)-1;

j=1;
mi=nanmin(eventInds);
if indsBeforeEvent>mi
    indsBeforeEvent=mi-1;
end
alignedOut=nan(sum(~isnan(eventInds)),length(mi-indsBeforeEvent:mi+howMuchToTake-1));
for i=1:length(eventInds)
    if ~isnan(eventInds(i))
        if eventInds(i)+howMuchToTake-1>size(whatToAlign,2)
            alignedOut(j,1:length(whatToAlign(fromInputRow(i),eventInds(i)-indsBeforeEvent:end)))=whatToAlign(fromInputRow(i),eventInds(i)-indsBeforeEvent:end);
        else
            alignedOut(j,:)=whatToAlign(fromInputRow(i),eventInds(i)-indsBeforeEvent:eventInds(i)+howMuchToTake-1);
        end
        j=j+1;
    end
end

end

function [alignedOut,eventInds,fromInputRow]=alignToEventFirstInTrial(event,withinRange_inds,whatToAlign,howMuchToTake,indsBeforeEvent)

% for each trial, get first event in this range
% align to this event

if withinRange_inds(1)<1
    withinRange_inds(1)=1;
end
if withinRange_inds(2)>size(event,2)
    withinRange_inds(2)=size(event,2);
end
event=event(:,withinRange_inds(1):withinRange_inds(2));

eventInds=nan(1,size(event,1));
fromInputRow=nan(1,size(event,1));
for i=1:size(event,1)
    if any(event(i,:)>0.5,2)
        f=find(event(i,:)>0.5,1,'first');
        eventInds(i)=f;
        fromInputRow(i)=i;
    end
end
eventInds=eventInds+withinRange_inds(1)-1;

j=1;
mi=nanmin(eventInds);
if indsBeforeEvent>mi
    indsBeforeEvent=mi-1;
end
alignedOut=nan(sum(~isnan(eventInds)),length(mi-indsBeforeEvent:mi+howMuchToTake-1));
for i=1:size(event,1)
    if ~isnan(eventInds(i))
        if eventInds(i)+howMuchToTake-1>size(whatToAlign,2)
            alignedOut(j,1:length(whatToAlign(i,eventInds(i)-indsBeforeEvent:end)))=whatToAlign(i,eventInds(i)-indsBeforeEvent:end);
        else
            alignedOut(j,:)=whatToAlign(i,eventInds(i)-indsBeforeEvent:eventInds(i)+howMuchToTake-1);
        end
        j=j+1;
    end
end

end

function [alignedOut,eventInds,fromInputRow]=alignToEventLastInTrial(event,withinRange_inds,whatToAlign,howMuchToTake,indsBeforeEvent)

% for each trial, get first event in this range
% align to this event

if withinRange_inds(1)<1
    withinRange_inds(1)=1;
end
if withinRange_inds(2)>size(event,2)
    withinRange_inds(2)=size(event,2);
end
event=event(:,withinRange_inds(1):withinRange_inds(2));

eventInds=nan(1,size(event,1));
fromInputRow=nan(1,size(event,1));
for i=1:size(event,1)
    if any(event(i,:)>0.5,2)
        f=find(event(i,:)>0.5,1,'last');
        eventInds(i)=f;
        fromInputRow(i)=i;
    end
end
eventInds=eventInds+withinRange_inds(1)-1;

j=1;
mi=nanmin(eventInds);
if indsBeforeEvent>mi
    indsBeforeEvent=mi-1;
end
alignedOut=nan(sum(~isnan(eventInds)),length(mi-indsBeforeEvent:mi+howMuchToTake-1));
for i=1:size(event,1)
    if ~isnan(eventInds(i))
        if eventInds(i)+howMuchToTake-1>size(whatToAlign,2)
            alignedOut(j,1:length(whatToAlign(i,eventInds(i)-indsBeforeEvent:end)))=whatToAlign(i,eventInds(i)-indsBeforeEvent:end);
        else
            alignedOut(j,:)=whatToAlign(i,eventInds(i)-indsBeforeEvent:eventInds(i)+howMuchToTake-1);
        end
        j=j+1;
    end
end

end

function plotWStderr(dataMatrix,times,c,plotUntilInd,nEvents)

showStdevInstead=false;

if isempty(plotUntilInd)
    plotUntilInd=size(dataMatrix,2);
end

plot(times(1:plotUntilInd),nanmean(dataMatrix(:,1:plotUntilInd),1),'Color',c,'LineWidth',1);
hold on;
if showStdevInstead==true
    plot(times(1:plotUntilInd),nanmean(dataMatrix(:,1:plotUntilInd),1)-nanstd(dataMatrix(:,1:plotUntilInd),[],1),'Color',c,'LineWidth',0.5);
    plot(times(1:plotUntilInd),nanmean(dataMatrix(:,1:plotUntilInd),1)+nanstd(dataMatrix(:,1:plotUntilInd),[],1),'Color',c,'LineWidth',0.5);
else
    plot(times(1:plotUntilInd),nanmean(dataMatrix(:,1:plotUntilInd),1)-nanstd(dataMatrix(:,1:plotUntilInd),[],1)./sqrt(nEvents),'Color',c,'LineWidth',0.5);
    plot(times(1:plotUntilInd),nanmean(dataMatrix(:,1:plotUntilInd),1)+nanstd(dataMatrix(:,1:plotUntilInd),[],1)./sqrt(nEvents),'Color',c,'LineWidth',0.5);
end
    
end