function phys_tbt=alignToOpto(phys_tbt)

minITI=9; % in seconds

thesePhotoFieldsUseTimeField1={'cue','cuetimes_wrt_trial_start'};
timeField1='cue_times';
f=fieldnames(phys_tbt);
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

if strcmp(timeField1,'cue_times') 
    if isfield(phys_tbt,'cuetimes_wrt_trial_start')
        phys_tbt.([timeField1 '_wrt_trial_start'])=phys_tbt.('cuetimes_wrt_trial_start');
    else
        temp=phys_tbt.(timeField1);
        phys_tbt.([timeField1 '_wrt_trial_start'])=temp-repmat(temp(:,1),1,size(temp,2));
    end
else
    temp=phys_tbt.(timeField1);
    phys_tbt.([timeField1 '_wrt_trial_start'])=temp-repmat(temp(:,1),1,size(temp,2));
end
temp=phys_tbt.(timeField2);
phys_tbt.([timeField2 '_wrt_trial_start'])=temp-repmat(temp(:,1),1,size(temp,2));

f=fieldnames(phys_tbt);
for i=1:length(f)
    if size(phys_tbt.(f{i}),2)==size(phys_tbt.(thesePhotoFieldsUseTimeField1{1}),2)
        if ~ismember(f{i},thesePhotoFieldsUseTimeField1)
            thesePhotoFieldsUseTimeField1{length(thesePhotoFieldsUseTimeField1)+1}=f{i};
        end
    end
end

[phys_tbt,alignedOptoTo,hasOpto]=realignToCue(phys_tbt,'opto',thesePhotoFieldsUseTimeField1,timeField1,timeField2,minITI);
phys_tbt.hasOpto=hasOpto;

% Compare firing rate at opto onset versus firing rate at baseline
phys_tbt=isOptoTagged(phys_tbt,hasOpto,thesePhotoFieldsUseTimeField2);

end

function phys_tbt=isOptoTagged(phys_tbt,hasOpto,thesePhotoFieldsUseTimeField2)

% binomial in a bin at baseline 
% vs binomial in same bin at opto onset

binsize=20; % in ms
binsize=binsize/1000; % in sec
% for baseline, two approaches:
% 1. take all time windows before opto 
% 2. same time point as opto but on trials where no opto

% get index and time of opto onset after have aligned all trials
maind=find(mean(phys_tbt.opto,1,'omitnan')>0.01,1,'first');
optoOff=find(mean(phys_tbt.opto(hasOpto==1,maind:end),1,'omitnan')<0.01,1,'first')+maind-1;
timestep=mode(diff(mean(phys_tbt.cue_times,1,'omitnan')));
times_wrt_trial_start=0:timestep:(size(phys_tbt.cue_times,2)-1)*timestep;
optoTime=times_wrt_trial_start(maind);
optoOffTime=times_wrt_trial_start(optoOff);
timestep=mode(diff(mean(phys_tbt.unitTimes_wrt_trial_start,1,'omitnan')));
unittimes_wrt_trial_start=0:timestep:(size(phys_tbt.unitTimes_wrt_trial_start,2)-1)*timestep;
[~,optounitind]=min(abs(mean(unittimes_wrt_trial_start,1,'omitnan')-optoTime));
[~,optooffunitind]=min(abs(mean(unittimes_wrt_trial_start,1,'omitnan')-optoOffTime));

nindsperbin=ceil(binsize/(unittimes_wrt_trial_start(2)-unittimes_wrt_trial_start(1)));
binsForBaseline=1:nindsperbin:optounitind-1;
binsDuringOpto=optounitind:nindsperbin:optooffunitind-1;
for i=1:length(thesePhotoFieldsUseTimeField2)
    if strcmp(thesePhotoFieldsUseTimeField2{i},'unitsum')
        % skip
        continue
    end
    nCountsPerBin_duringOptoOnset=nan(sum(hasOpto,2,'omitnan'),length(binsForBaseline)-1);
    nCountsPerBin_duringOptoNoOnset=nan(sum(hasOpto,2,'omitnan'),length(binsForBaseline)-1);
    nCountsPerBin_beforeOptoOnset=nan(sum(hasOpto,2,'omitnan'),length(binsForBaseline)-1);
    temp=phys_tbt.(thesePhotoFieldsUseTimeField2{i});
    countsForApproach1Baseline=nan(size(temp,1),length(binsForBaseline)-1);
    for j=1:length(binsForBaseline)-1
        % approach 1 baseline
        currbin=binsForBaseline(j):binsForBaseline(j+1);
        countsForApproach1Baseline(:,j)=countEventsInBin(temp,currbin);
    end
    countsForApproach2Baseline=nan(sum(hasOpto==0),length(binsDuringOpto)-1);
    for j=1:length(binsDuringOpto)-1
        % approach 2 baseline
        currbin=binsDuringOpto(j):binsDuringOpto(j+1);
        countsForApproach2Baseline(:,j)=countEventsInBin(temp(hasOpto==0,:),currbin);
    end
    countsDuringOpto=nan(sum(hasOpto==1),length(binsDuringOpto)-1);
    pval_approach2=nan(1,length(binsDuringOpto)-1);
    pval_approach1=nan(1,length(binsDuringOpto)-1);
    for j=1:length(binsDuringOpto)-1
        currbin=binsDuringOpto(j):binsDuringOpto(j+1);
        countsDuringOpto(:,j)=countEventsInBin(temp(hasOpto==1,:),currbin);
        pval_approach2(j)=compareCounts(countsForApproach2Baseline(:,j),countsDuringOpto(:,j));
        pval_approach1(j)=compareCounts(countsForApproach1Baseline(1:end),countsDuringOpto(:,j));
    end
    phys_tbt.([thesePhotoFieldsUseTimeField2{i} '_pval1'])=pval_approach1;
    phys_tbt.([thesePhotoFieldsUseTimeField2{i} '_pval2'])=pval_approach2;
    phys_tbt.([thesePhotoFieldsUseTimeField2{i} '_avAlignedToOpto'])=mean(temp(hasOpto==1,:),1,'omitnan');
    if ~isfield(phys_tbt,'optoOnInUnitTimes')
        phys_tbt.optoOnInUnitTimes=zeros(size(phys_tbt.([thesePhotoFieldsUseTimeField2{i} '_avAlignedToOpto'])));
        phys_tbt.optoOnInUnitTimes(optounitind:optooffunitind-1)=1;
    end
end

end

function pval=compareCounts(base,test)

% binomial stats
% what is the probability of observing the number of spikes during opto 
% if the real probability of a spike were unchanged from baseline
% one-sided test that opto spiking is greater than baseline

p=sum(base>0)/length(base);
k=sum(test>0);
n=length(test);
pval=0;
for i=k:n
    pval=pval+binopdf(i,n,p);
end

end

function out=countEventsInBin(dataMatrix,binInInds)

dataMatrix(dataMatrix>0)=1;
out=sum(dataMatrix(:,binInInds),2,'omitnan');

end

function [data,alignedCueTo,hasCue]=realignToCue(data,cueField,fieldsLikeCue,cueTimesName,otherTimesName,minITI)

temp=nanmean(data.(cueField),1);
[ma,cueInd]=nanmax(temp); % align all to this mode
[~,fbeyond]=nanmin(abs(nanmean(data.('cue_times_wrt_trial_start'),1)-minITI));
if cueInd>fbeyond
    % cue max cannot occur after minITI
    tempie=data.(cueField);
    tempie(:,fbeyond:end)=0;
    [ma,cueInd]=nanmax(nanmean(tempie,1)); % align all to this mode
end
% cueHalfMax=ma/2;
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
hasCue=zeros(1,size(tempdata,1));
for i=1:size(tempdata,1)
    tempdata(i,fbeyond:end)=0;
    f=find(tempdata(i,:)>0.01,1,'first');
    if isempty(f)
        continue
    end
    hasCue(i)=1;
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