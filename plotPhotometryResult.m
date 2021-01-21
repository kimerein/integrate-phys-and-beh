function plotPhotometryResult(photometry_tbt,behavior_tbt,outdata,alignTo,plotPhotoField,plotBehField)

% alignTo can be 'cue','success_reachStarts'

cutBeforeNextCue=true; % if is true, will only plot inter-trial interval after first cue
minITI=9; % in seconds
behaviorCue='cueZone_onVoff';

thesePhotoFieldsUseTimeField1={'green_mod','red_mod','opto','cue','cue_times','distractor','from_first_video','from_second_video'};
timeField1='cue_times';
thesePhotoFieldsUseTimeField2={'green_ch','red_ch','green_time','raw_green_ch','red_time','raw_red_ch'};
timeField2='red_time';

temp=photometry_tbt.(timeField1);
photometry_tbt.([timeField1 '_wrt_trial_start'])=temp-repmat(temp(:,1),1,size(temp,2));
temp=photometry_tbt.(timeField2);
photometry_tbt.([timeField2 '_wrt_trial_start'])=temp-repmat(temp(:,1),1,size(temp,2));

photometry_tbt=realignToCue(photometry_tbt,'cue',thesePhotoFieldsUseTimeField1,timeField1,timeField2);
% figure();
% plot(nanmean(photometry_tbt.cue,1),'Color','k');
% title('realigned to cue');

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
    firstSearchInd=floor(minITI./mode(photometry_tbt.cue_times(2,2:end)-photometry_tbt.cue_times(2,1:end-1)));
    f=find(nanmean(photometry_tbt.cue(:,firstSearchInd:end),1)>0.01,1,'first');
    f=f+firstSearchInd;
    
    % find first non-nan column
    ftostart=find(~isnan(photometry_tbt.cue_times(1,:)),1,'first');
    temp=photometry_tbt.(getCorrectTime);
    % find first non-nan column
    ftostart_photo=find(~isnan(temp(1,:)),1,'first');
    
    f_times=nanmean(photometry_tbt.cue_times(:,f)-photometry_tbt.cue_times(:,ftostart),1);
    [~,f_photometry_ind]=nanmin(abs(nanmean(photometry_tbt.green_time-repmat(photometry_tbt.green_time(:,ftostart_photo),1,size(photometry_tbt.green_time,2)),1)-f_times));
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

typeOfReach=false;
switch alignTo
    case 'cue'
        figure();
        plotWStderr(photometry_tbt.(plotPhotoField),phototimes,'k',plotUntilInd,size(photometry_tbt.cue,1));
        hold on;
        plotWStderr(nanmean(photometry_tbt.cue,1),nanmean(photometry_tbt.cue_times-repmat(photometry_tbt.cue_times(:,1),1,size(photometry_tbt.cue_times,2)),1),'b',f,size(photometry_tbt.cue,1));
    case 'success_fromPerchOrWheel'
        typeOfReach=true;
        useReach='combo';
        useCombo=behavior_tbt.success_reachStarts+behavior_tbt.success_reachStarts_pawOnWheel;
    case 'drop_fromPerchOrWheel'
        typeOfReach=true;
        useReach='combo';
        useCombo=behavior_tbt.drop_reachStarts+behavior_tbt.drop_reachStarts_pawOnWheel;
    case 'miss_fromPerchOrWheel'
        typeOfReach=true;
        useReach='combo';
        useCombo=behavior_tbt.miss_reachStarts+behavior_tbt.miss_reachStarts_pawOnWheel;
    otherwise 
        typeOfReach=true;
        useReach=alignTo;
end
if typeOfReach==true
    % in behavior_tbt
    if ischar(useReach)
        temp=useCombo;
    else
        temp=behavior_tbt.(useReach);
    end
    % note that pellet is only present after cue
    withinRange=[-1.5 16];  % take reaches in this range, time is in seconds from onset of cue
%     withinRange=[5 16];  % take reaches in this range, time is in seconds from onset of cue
%     withinRange=[0 5];  % take reaches in this range, time is in seconds from onset of cue
    nIndsToTake=200;
    timeStep=mode(diff(nanmean(behavior_tbt.times,1)));
    cueInd=find(nanmean(behavior_tbt.(behaviorCue),1)>0,1,'first');
    withinRange_inds=[cueInd+ceil(withinRange(1)/timeStep) cueInd+ceil(withinRange(2)/timeStep)];
    [behAligned,alignInds]=alignToEvent(temp,withinRange_inds,temp,nIndsToTake,90);
    alignTimes=getValsFromRows(behavior_tbt.times_wrt_trial_start,alignInds);
    indsIntoPhoto=getIndsFromRows(photometry_tbt.(getCorrectTime_wrtTrialStart),alignTimes);
    alignedData=alignRowsToInds(photometry_tbt.(plotPhotoField),indsIntoPhoto,nanmin(indsIntoPhoto));
    figure();
    plotWStderr(alignedData,phototimes(1:size(alignedData,2)),'k',[],size(behAligned,1));
    hold on;
    plotWStderr(behAligned,nanmean(behavior_tbt.times_wrt_trial_start(:,1:size(behAligned,2)),1),'g',[],size(behAligned,1));
    disp([num2str(size(behAligned,1)) ' events averaged']);
end

end

function data=realignToCue(data,cueField,fieldsLikeCue,cueTimesName,otherTimesName)

[ma,cueInd]=nanmax(nanmean(data.(cueField),1)); % align all to this mode
cueHalfMax=ma/2;
temp=data.(cueTimesName);
tim=nanmean(temp-repmat(temp(:,1),1,size(temp,2)),1);
if 100>length(tim)
    ts_cue=mode(tim(2:end)-tim(1:end-1)); % cue times
else
    ts_cue=mode(tim(2:100)-tim(1:100-1)); % cue times
end
temp=data.(otherTimesName);
tim=nanmean(temp-repmat(temp(:,1),1,size(temp,2)),1);
ts_other=mode(tim(2:end)-tim(1:end-1)); % other fields times

tempdata=data.(cueField);
cueLikeShiftInds=nan(1,size(tempdata,1));
otherLikeShiftInds=nan(1,size(tempdata,1));
for i=1:size(tempdata,1)
    f=find(tempdata(i,:)>cueHalfMax,1,'first');
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

function alignedData=alignRowsToInds(data,alignTo,indsBeforeAlign)

mi=nanmin(alignTo);
if indsBeforeAlign>=mi
    indsBeforeAlign=mi-1;
end
alignedData=nan(size(data,1),indsBeforeAlign+size(data,2)-mi);
for i=1:length(alignTo)
    if isnan(alignTo(i))
        continue
    end
    alignedData(i,1:length(data(i,alignTo(i)-indsBeforeAlign:end)))=data(i,alignTo(i)-indsBeforeAlign:end);
end

end

function out=getIndsFromRows(dataVals,valsInData)

out=nan(1,length(valsInData));
for i=1:length(valsInData)
    if ~isnan(valsInData(i))
        [~,out(i)]=nanmin(abs(dataVals(i,:)-valsInData(i))); 
    end
end 

end

function out=getValsFromRows(dataVals,indsIntoRows)

out=nan(1,length(indsIntoRows));
for i=1:length(indsIntoRows)
    if ~isnan(indsIntoRows(i))
        out(i)=dataVals(i,indsIntoRows(i));
    end
end

end

function [alignedOut,eventInds]=alignToEvent(event,withinRange_inds,whatToAlign,howMuchToTake,indsBeforeEvent)

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
for i=1:size(event,1)
    if any(event(i,:)>0.5,2)
        f=find(event(i,:)>0.5,1,'first');
        eventInds(i)=f;
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