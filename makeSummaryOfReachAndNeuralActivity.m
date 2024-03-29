function physiology_tbt=makeSummaryOfReachAndNeuralActivity(physiology_tbt,phys_beh_tbt,photometry_tbt,photo_beh_tbt,spikesDir,spikes,skipSpikes)

% skipSpikes is true if physiology_tbt already has unit fields

spikeName='spikes';
binsize=10; % in ms, for PSTHs
goodUnitLabel=2; 
bsmooth=1;
maxTrialLength=9; % in seconds
normalizeSU=true;

% load spikes from directory
spike_d=dir(spikesDir);
all_wvfms=[];
all_halfWidths=[];
all_depths=[];
all_assigns=[];
unitnames={};
whichTrode=[];
firstloadedspikes=true;
if ~isempty(spikes) || skipSpikes==true
    if skipSpikes==true && isfield(physiology_tbt,'sum_over_singleunit')
        % physiology_tbt already has unit fields
        unitnames=physiology_tbt.details.unitnames;
        all_wvfms=physiology_tbt.details.wvfms;
        all_halfWidths=physiology_tbt.details.halfWidths;
        all_depths=physiology_tbt.details.depths;
        whichTrode=physiology_tbt.details.trodes;
    else
        if ~isfield(spikes,'sweeps')
            [~,spikes]=organizeSpikesToMatch_physiology_tbt(spikes,physiology_tbt);
        end
        % get PSTH for each unit
        if firstloadedspikes==true
            [physiology_tbt,~,unit_fieldnames]=addGoodUnitsAsFields(physiology_tbt,spikes,goodUnitLabel,binsize,bsmooth,true);
            firstloadedspikes=false;
        else
            [physiology_tbt,~,unit_fieldnames]=addGoodUnitsAsFields(physiology_tbt,spikes,goodUnitLabel,binsize,bsmooth,false);
        end
        if ~isfield(physiology_tbt,'sum_over_singleunit')
            physiology_tbt.sum_over_singleunit=zeros(size(physiology_tbt.unitsum));
        end
        tmp=cat(3,physiology_tbt.sum_over_singleunit,physiology_tbt.unitsum); 
        physiology_tbt.sum_over_singleunit=nansum(tmp,3);
        unitnames(length(unitnames)+1:length(unitnames)+length(unit_fieldnames))=unit_fieldnames;
        % get info to classify units
        [unit_wvfms,unit_halfWidths,unit_depths,useAssigns]=getUnitWaveformAndDepth(spikes);
        all_wvfms=[all_wvfms; unit_wvfms(spikes.labels(:,2)==goodUnitLabel,:)];
        all_assigns=[all_assigns useAssigns(spikes.labels(:,2)==goodUnitLabel)];
        all_halfWidths=[all_halfWidths unit_halfWidths(spikes.labels(:,2)==goodUnitLabel)];
        all_depths=[all_depths unit_depths(spikes.labels(:,2)==goodUnitLabel)];
        whichTrode=[whichTrode ones(size(unit_depths(spikes.labels(:,2)==goodUnitLabel)))*1];
        physiology_tbt.details.wvfms=all_wvfms;
        physiology_tbt.details.halfWidths=all_halfWidths;
        physiology_tbt.details.depths=all_depths;
        physiology_tbt.details.trodes=whichTrode;
        physiology_tbt.details.unitnames=unitnames;
    end
else
    for i=1:length(spike_d)
        spind=regexp(spike_d(i).name,spikeName);
        if ~isempty(spind)
            % load spikes
            [sta,en]=regexp(spike_d(i).name,'(?<name>\d+)');
            if isempty(sta)
                currtrode=1;
            else
                currtrode=str2double(spike_d(i).name(sta:en));
            end
            disp(['loading ' spike_d(i).name]);
            a=load([spikesDir '\' spike_d(i).name]);
            spikes=a.spikes;
        else
            continue
        end
        [~,spikes]=organizeSpikesToMatch_physiology_tbt(spikes,physiology_tbt);
        % get PSTH for each unit
        if firstloadedspikes==true
            [physiology_tbt,~,unit_fieldnames]=addGoodUnitsAsFields(physiology_tbt,spikes,goodUnitLabel,binsize,bsmooth,true);
            firstloadedspikes=false;
        else
            [physiology_tbt,~,unit_fieldnames]=addGoodUnitsAsFields(physiology_tbt,spikes,goodUnitLabel,binsize,bsmooth,false);
        end
        if ~isfield(physiology_tbt,'sum_over_singleunit')
            physiology_tbt.sum_over_singleunit=zeros(size(physiology_tbt.unitsum));
        end
        tmp=cat(3,physiology_tbt.sum_over_singleunit,physiology_tbt.unitsum); 
        physiology_tbt.sum_over_singleunit=nansum(tmp,3);
        unitnames(length(unitnames)+1:length(unitnames)+length(unit_fieldnames))=unit_fieldnames;
        % get info to classify units
        % if spikes has fewer than 4 ev channels, fix this
        if size(spikes.waveforms,3)<4
            for j=size(spikes.waveforms,3)+1:4
                spikes.waveforms(:,:,j)=zeros(size(spikes.waveforms(:,:,1)));
            end
        end
        [unit_wvfms,unit_halfWidths,unit_depths,useAssigns]=getUnitWaveformAndDepth(spikes);
        all_wvfms=[all_wvfms; unit_wvfms(spikes.labels(:,2)==goodUnitLabel,:)];
        all_assigns=[all_assigns useAssigns(spikes.labels(:,2)==goodUnitLabel)];
        all_halfWidths=[all_halfWidths; unit_halfWidths(spikes.labels(:,2)==goodUnitLabel)];
        all_depths=[all_depths; unit_depths(spikes.labels(:,2)==goodUnitLabel)];
        whichTrode=[whichTrode; ones(size(unit_depths(spikes.labels(:,2)==goodUnitLabel)))*currtrode];
    end
    physiology_tbt.details.wvfms=all_wvfms;
    physiology_tbt.details.halfWidths=all_halfWidths;
    physiology_tbt.details.depths=all_depths;
    physiology_tbt.details.trodes=whichTrode;
    physiology_tbt.details.unitnames=unitnames;
end
physiology_tbt.av_over_singleunit=physiology_tbt.sum_over_singleunit./length(physiology_tbt.details.unitnames);

% downSamp=20; % downSamp=6; % if don't want to downSamp, set to 1, else set to index binsize for downSamp
downSamp=6; % if don't want to downSamp, set to 1, else set to index binsize for downSamp
% downSamp=1; % if don't want to downSamp, set to 1, else set to index binsize for downSamp
physiology_tbt=makeUnitSubgroups(physiology_tbt,downSamp);

% physiology_tbt=addUnitSubgroups(physiology_tbt,downSamp,physiology_tbt.details.whichDAgrp,physiology_tbt.details.whichgrp==1,'DA');
% physiology_tbt=addUnitSubgroups(physiology_tbt,downSamp,physiology_tbt.details.whichDAgrp,physiology_tbt.details.whichgrp==1,'negDA');
% temp=physiology_tbt.details.whichDAgrp;
% temp(physiology_tbt.details.whichDAgrp~=physiology_tbt.details.whichnegDAgrp)=3;
% physiology_tbt.details.whichConsensusDAgrp=temp;
% physiology_tbt=addUnitSubgroups(physiology_tbt,downSamp,physiology_tbt.details.whichConsensusDAgrp,physiology_tbt.details.whichgrp==1,'consensusDA');

% figure 5
% photothresh=0.75;
% lowPassCutoff=5; % in Hz
% dosmooth=true;
% smoothFields={'green_ch'};
% triggerOnField=smoothFields{1};
% if dosmooth==true
%     phototimes=nanmean(photometry_tbt.times_wrt_trial_start,1);
%     photometry_tbt=smoothPhotometry(photometry_tbt,1/mode(diff(phototimes)),lowPassCutoff,smoothFields);
%     disp(['using photometry Fs ' num2str(1/mode(diff(phototimes)))]);
% end
% photometry_tbt.isPhotoEvent=double(photometry_tbt.(triggerOnField)>photothresh);
% photo_beh_tbt=makeSameFieldInBeh(photometry_tbt,photo_beh_tbt,'isPhotoEvent',photometry_tbt.times_wrt_trial_start,photo_beh_tbt.times_wrt_trial_start);
% phys_beh_tbt=putPhotoBehFieldIntoPhysBeh(photo_beh_tbt,phys_beh_tbt,'isPhotoEvent');
% alignments={'isPhotoEvent'};
% xranges={[0 maxTrialLength]};
% timewindows={[-0.25 16]};
% withintimewindow={'last'};
% % beh_fields={'all_reachBatch','fidgetData'};
% % photo_fields={'green_ch'};
% phys_fields={'unit_by_unit','av_over_singleunit','grp1_unitav','grp2_unitav','grp3_unitav'};
% beh_fields={'all_reachBatch'};
% photo_fields={'green_ch'};
% % phys_fields={'unit_by_unit'};
% physthenphoto_fields(1:length(phys_fields))=phys_fields;
% physthenphoto_fields(length(phys_fields)+1:length(phys_fields)+length(photo_fields))=photo_fields;
% isPhysField=[ones(size(1:length(phys_fields))) zeros(size(1:length(photo_fields)))];
% outSU=makeSummaryFig(beh_fields,photo_fields,phys_fields,alignments,physthenphoto_fields,withintimewindow,timewindows,isPhysField,xranges,photometry_tbt,photo_beh_tbt,physiology_tbt,phys_beh_tbt,normalizeSU,maxTrialLength);
% firstDAchange=classifyUnitsAsDAincreaseOrDecrease(outSU);
% timewindows={[-0.25 16]};
% withintimewindow={'first'};
% beh_fields={'all_reachBatch','fidgetData'};
% photo_fields={'green_ch'};
% phys_fields={'unit_by_unit','av_over_singleunit','grp1_unitav','grp2_unitav','grp3_unitav'};
% % beh_fields={'all_reachBatch'};
% % photo_fields={'green_ch'};
% % phys_fields={'unit_by_unit'};
% physthenphoto_fields(1:length(phys_fields))=phys_fields;
% physthenphoto_fields(length(phys_fields)+1:length(phys_fields)+length(photo_fields))=photo_fields;
% isPhysField=[ones(size(1:length(phys_fields))) zeros(size(1:length(photo_fields)))];
% outSU=makeSummaryFig(beh_fields,photo_fields,phys_fields,alignments,physthenphoto_fields,withintimewindow,timewindows,isPhysField,xranges,photometry_tbt,photo_beh_tbt,physiology_tbt,phys_beh_tbt,normalizeSU,maxTrialLength);
% % lastDAchange=classifyUnitsAsDAincreaseOrDecrease(outSU);
% % DAfieldname='whichDAgrp';
% % physiology_tbt=compareFirstAndLastDAchange(firstDAchange,lastDAchange,physiology_tbt.details.whichgrp,physiology_tbt,DAfieldname);
% % physiology_tbt=addUnitSubgroups(physiology_tbt,1,physiology_tbt.details.whichDAgrp,physiology_tbt.details.whichgrp==1,'DA');

% figure 6
% photothresh=-0.75;
% lowPassCutoff=5; % in Hz
% dosmooth=true;
% smoothFields={'green_ch'};
% triggerOnField=smoothFields{1};
% if dosmooth==true
%     phototimes=nanmean(photometry_tbt.times_wrt_trial_start,1);
%     photometry_tbt=smoothPhotometry(photometry_tbt,1/mode(diff(phototimes)),lowPassCutoff,smoothFields);
%     disp(['using photometry Fs ' num2str(1/mode(diff(phototimes)))]);
% end
% photometry_tbt.isPhotoEvent=double(photometry_tbt.(triggerOnField)<photothresh);
% photo_beh_tbt=makeSameFieldInBeh(photometry_tbt,photo_beh_tbt,'isPhotoEvent',photometry_tbt.times_wrt_trial_start,photo_beh_tbt.times_wrt_trial_start);
% phys_beh_tbt=putPhotoBehFieldIntoPhysBeh(photo_beh_tbt,phys_beh_tbt,'isPhotoEvent');
% alignments={'isPhotoEvent'};
% xranges={[0 maxTrialLength]};
% timewindows={[-0.25 16]};
% withintimewindow={'first'};
% % beh_fields={'all_reachBatch','fidgetData'};
% % photo_fields={'green_ch'};
% % phys_fields={'unit_by_unit','av_over_singleunit','grp1_unitav','grp2_unitav','grp3_unitav'};
% beh_fields={'all_reachBatch'};
% photo_fields={'green_ch'};
% phys_fields={'unit_by_unit'};
% physthenphoto_fields(1:length(phys_fields))=phys_fields;
% physthenphoto_fields(length(phys_fields)+1:length(phys_fields)+length(photo_fields))=photo_fields;
% isPhysField=[ones(size(1:length(phys_fields))) zeros(size(1:length(photo_fields)))];
% outSU=makeSummaryFig(beh_fields,photo_fields,phys_fields,alignments,physthenphoto_fields,withintimewindow,timewindows,isPhysField,xranges,photometry_tbt,photo_beh_tbt,physiology_tbt,phys_beh_tbt,normalizeSU,maxTrialLength);
% firstDAchange=classifyUnitsAsDAincreaseOrDecrease(outSU);
% timewindows={[-0.25 16]};
% withintimewindow={'last'};
% % beh_fields={'all_reachBatch','fidgetData'};
% % photo_fields={'green_ch'};
% % phys_fields={'unit_by_unit','av_over_singleunit','grp1_unitav','grp2_unitav','grp3_unitav'};
% beh_fields={'all_reachBatch'};
% photo_fields={'green_ch'};
% phys_fields={'unit_by_unit'};
% physthenphoto_fields(1:length(phys_fields))=phys_fields;
% physthenphoto_fields(length(phys_fields)+1:length(phys_fields)+length(photo_fields))=photo_fields;
% isPhysField=[ones(size(1:length(phys_fields))) zeros(size(1:length(photo_fields)))];
% outSU=makeSummaryFig(beh_fields,photo_fields,phys_fields,alignments,physthenphoto_fields,withintimewindow,timewindows,isPhysField,xranges,photometry_tbt,photo_beh_tbt,physiology_tbt,phys_beh_tbt,normalizeSU,maxTrialLength);
% lastDAchange=classifyUnitsAsDAincreaseOrDecrease(outSU);
% DAfieldname='whichnegDAgrp';
% physiology_tbt=compareFirstAndLastDAchange(firstDAchange,lastDAchange,physiology_tbt.details.whichgrp,physiology_tbt,DAfieldname);
% physiology_tbt=addUnitSubgroups(physiology_tbt,1,physiology_tbt.details.whichnegDAgrp,physiology_tbt.details.whichgrp==1,'negDA');

% figure 1
% alignments={'cue','all_reachBatch'};
% alignments={'cue'};
% timewindows={[-1 16]};
% withintimewindow={[]};
% beh_fields={'all_reachBatch','fidgetData'};
% photo_fields={'green_ch'};
% phys_fields={'unit_by_unit','av_over_singleunit','grp1_unitav','grp2_unitav','grp3_unitav'}; % each field must contain "unit"
% xranges={[0 maxTrialLength],[0 maxTrialLength]};
% physthenphoto_fields(1:length(phys_fields))=phys_fields;
% physthenphoto_fields(length(phys_fields)+1:length(phys_fields)+length(photo_fields))=photo_fields;
% isPhysField=[ones(size(1:length(phys_fields))) zeros(size(1:length(photo_fields)))];
% makeSummaryFig(beh_fields,photo_fields,phys_fields,alignments,physthenphoto_fields,withintimewindow,timewindows,isPhysField,xranges,photometry_tbt,photo_beh_tbt,physiology_tbt,phys_beh_tbt,normalizeSU,maxTrialLength);

% figure 2
alignments={'success_fromPerchOrWheel','success batch when pellet dislodged','drop_fromPerchOrWheel','misses_and_pelletMissing'};
alignments={'success_fromPerchOrWheel','success batch when pellet dislodged','drop batch when pellet dislodged','pelletmissingreach_reachStarts'};
xranges={[0 maxTrialLength],[0 maxTrialLength],[0 maxTrialLength],[0 maxTrialLength]};
timewindows={[-1 16],[-1 16],[-1 16],[-1 16]};
withintimewindow={'first','first','first','first'};
beh_fields={'all_reachBatch','fidgetData'};
photo_fields={'green_ch'};
% phys_fields={'unit_by_unit','av_over_singleunit','grp1_unitav','grp2_unitav','grp3_unitav','grpDA1_unitav','grpDA2_unitav','grpDA3_unitav'};
phys_fields={'unit_by_unit','av_over_singleunit','grp1_unitav','grp2_unitav','grp3_unitav'}; %,'grpconsensusDA1_unitav','grpconsensusDA2_unitav','grpconsensusDA3_unitav'};
physthenphoto_fields(1:length(phys_fields))=phys_fields;
physthenphoto_fields(length(phys_fields)+1:length(phys_fields)+length(photo_fields))=photo_fields;
isPhysField=[ones(size(1:length(phys_fields))) zeros(size(1:length(photo_fields)))];
makeSummaryFig(beh_fields,photo_fields,phys_fields,alignments,physthenphoto_fields,withintimewindow,timewindows,isPhysField,xranges,photometry_tbt,photo_beh_tbt,physiology_tbt,phys_beh_tbt,normalizeSU,maxTrialLength);

% figure 3
% alignments={'success_fromPerchOrWheel','drop_fromPerchOrWheel','misses_and_pelletMissing'};
% xranges={[0 maxTrialLength],[0 maxTrialLength],[0 maxTrialLength]};
% timewindows={[0 3],[0 3],[0 3]};
% withintimewindow={'first','first','first'};
% beh_fields={'all_reachBatch','fidgetData'};
% photo_fields={'green_ch'};
% phys_fields={'unit_by_unit','av_over_singleunit','grp1_unitav','grp2_unitav','grp3_unitav'};
% physthenphoto_fields(1:length(phys_fields))=phys_fields;
% physthenphoto_fields(length(phys_fields)+1:length(phys_fields)+length(photo_fields))=photo_fields;
% isPhysField=[ones(size(1:length(phys_fields))) zeros(size(1:length(photo_fields)))];
% makeSummaryFig(beh_fields,photo_fields,phys_fields,alignments,physthenphoto_fields,withintimewindow,timewindows,isPhysField,xranges,photometry_tbt,photo_beh_tbt,physiology_tbt,phys_beh_tbt,normalizeSU,maxTrialLength);

% figure 3.5
% alignments={'success batch when pellet dislodged','success batch when pellet dislodged','success batch when pellet dislodged','misses_and_pelletMissing','misses_and_pelletMissing','misses_and_pelletMissing'};
% xranges={[0 maxTrialLength],[0 maxTrialLength],[0 maxTrialLength],[0 maxTrialLength],[0 maxTrialLength],[0 maxTrialLength]};
% timewindows={[-1.5 0],[0 1],[1 16],[-1.5 0],[0 1],[1 16]};
% withintimewindow={'first','first','first','first','first','first'};
% beh_fields={'all_reachBatch','fidgetData'};
% photo_fields={'green_ch'};
% phys_fields={'unit_by_unit','av_over_singleunit','grp1_unitav','grp2_unitav','grp3_unitav'};
% physthenphoto_fields(1:length(phys_fields))=phys_fields;
% physthenphoto_fields(length(phys_fields)+1:length(phys_fields)+length(photo_fields))=photo_fields;
% isPhysField=[ones(size(1:length(phys_fields))) zeros(size(1:length(photo_fields)))];
% makeSummaryFig(beh_fields,photo_fields,phys_fields,alignments,physthenphoto_fields,withintimewindow,timewindows,isPhysField,xranges,photometry_tbt,photo_beh_tbt,physiology_tbt,phys_beh_tbt,normalizeSU,maxTrialLength);

% figure 4
% alignments={'success_fromPerchOrWheel','drop_fromPerchOrWheel','misses_and_pelletMissing'};
% xranges={[0 maxTrialLength],[0 maxTrialLength],[0 maxTrialLength]};
% timewindows={[5 16],[5 16],[5 16]};
% withintimewindow={'first','first','first'};
% beh_fields={'all_reachBatch','fidgetData'};
% photo_fields={'green_ch'};
% phys_fields={'unit_by_unit','av_over_singleunit','grp1_unitav','grp2_unitav','grp3_unitav'};
% physthenphoto_fields(1:length(phys_fields))=phys_fields;
% physthenphoto_fields(length(phys_fields)+1:length(phys_fields)+length(photo_fields))=photo_fields;
% isPhysField=[ones(size(1:length(phys_fields))) zeros(size(1:length(photo_fields)))];
% makeSummaryFig(beh_fields,photo_fields,phys_fields,alignments,physthenphoto_fields,withintimewindow,timewindows,isPhysField,xranges,photometry_tbt,photo_beh_tbt,physiology_tbt,phys_beh_tbt,normalizeSU,maxTrialLength);

plotSpikeWvfmsByFR(all_wvfms,all_halfWidths,unitnames,physiology_tbt);

end

function physiology_tbt=compareFirstAndLastDAchange(firstDAchange,lastDAchange,whichUnitGrp,physiology_tbt,DAfieldname)

% out.changeIn1=changeIn1;
% out.changeIn2=changeIn2;
% out.after1=after1;
% out.after2=after2;
% out.before1=before1;
% out.before2=before2;

cmap=colormap('jet');
grp1c=cmap(1,:); grp2c=cmap(floor(size(cmap,1)/2),:); grp3c=cmap(end,:);
c=nan(length(whichUnitGrp),size(cmap,2));
for i=1:length(whichUnitGrp)
    switch whichUnitGrp(i)
        case 1
            c(i,:)=grp1c;
        case 2
            c(i,:)=grp2c;
        case 3
            c(i,:)=grp3c;
    end
end

% figure();
% scatter(firstDAchange.changeIn2./firstDAchange.before2,lastDAchange.changeIn2./lastDAchange.before2,[],c);
% xlabel('First post-DA change');
% ylabel('Last post-DA change');
% 
cs=cmap(1:floor(size(cmap,1)/length(whichUnitGrp)):end,:);
cs=cs(1:length(whichUnitGrp),:);
cs(whichUnitGrp~=1,:)=repmat([1 1 1],nansum(whichUnitGrp~=1),1);
% figure();
% scatter(firstDAchange.changeIn1./firstDAchange.before1,firstDAchange.changeIn2./firstDAchange.before2,100,cs);
% hold on;
% scatter(lastDAchange.changeIn1./lastDAchange.before1,lastDAchange.changeIn2./lastDAchange.before2,100,cs);
% 
% figure();
% scatter((firstDAchange.around1-firstDAchange.before1),(lastDAchange.around1-lastDAchange.before1),100,cs);
% 
% figure();
% scatter((firstDAchange.around1-firstDAchange.before1)./firstDAchange.before1,(lastDAchange.around1-lastDAchange.before1)./lastDAchange.before1,100,cs);

% figure(); 
% scatter((lastDAchange.around1-lastDAchange.before1)./lastDAchange.before1,firstDAchange.changeIn2./firstDAchange.before2,50,cs);
% grpA=((lastDAchange.around1-lastDAchange.before1)./lastDAchange.before1)>0.05 & (firstDAchange.changeIn2./firstDAchange.before2)>0.5;
% grpB=((lastDAchange.around1-lastDAchange.before1)./lastDAchange.before1)<0.05 & (firstDAchange.changeIn2./firstDAchange.before2)<=0.5;

% if strcmp(DAfieldname,'whichDAgrp')
%     figure();
%     scatter(lastDAchange.changeIn2./lastDAchange.before2,firstDAchange.changeIn2./firstDAchange.before2,50,cs);
%     grpA=(lastDAchange.changeIn2./lastDAchange.before2)>0;
%     grpB=(lastDAchange.changeIn2./lastDAchange.before2)<0;
%     grpC=~(grpA | grpB);
% elseif strcmp(DAfieldname,'whichnegDAgrp')
%     figure();
%     scatter(firstDAchange.changeIn2./firstDAchange.before2,lastDAchange.changeIn2./lastDAchange.before2,50,cs);
%     grpB=(firstDAchange.changeIn2./firstDAchange.before2)>0;
%     grpA=~grpB;
%     grpC=~(grpA | grpB);
% end

temp=nan(size(physiology_tbt.details.whichgrp));
temp(grpA==true)=1;
temp(grpB==true)=2;
temp(grpC==true)=3;
physiology_tbt.details.(DAfieldname)=temp;

end

function physiology_tbt=addUnitSubgroups(physiology_tbt,downSamp,whichgrp,whichToUse,subgroupName)

unitnames=physiology_tbt.details.unitnames;

for i=1:3 % 3 groups
    physiology_tbt.(['grp' subgroupName num2str(i) '_unitav'])=zeros(size(physiology_tbt.unitsum));
    for j=1:length(whichgrp)
        if whichToUse(j)==false
            continue
        end
        if whichgrp(j)==i
            tmp=cat(3,physiology_tbt.(['grp' subgroupName num2str(i) '_unitav']),physiology_tbt.(unitnames{j})); 
            physiology_tbt.(['grp' subgroupName num2str(i) '_unitav'])=nansum(tmp,3);
        end
    end
    physiology_tbt.(['grp' subgroupName num2str(i) '_unitav'])=physiology_tbt.(['grp' subgroupName num2str(i) '_unitav'])./nansum(whichgrp==i);
end

if downSamp~=1
    for i=1:3
        temp=downSampMatrix(physiology_tbt.(['grp' subgroupName num2str(i) '_unitav']),downSamp);
        reinterp_temp=interp1(downSampAv(nanmean(physiology_tbt.unitTimes,1),downSamp)',temp',nanmean(physiology_tbt.unitTimes,1)');
        physiology_tbt.(['grp' subgroupName num2str(i) '_unitav'])=reinterp_temp';
    end
end

end

function out=classifyUnitsAsDAincreaseOrDecrease(outSU)

% outSU.su
% outSU.times
% outSU.alignmentCompanion.y
% outSU.alignmentCompanion.x

% find time of DA event
outSU.alignmentCompanion.x=nanmean(outSU.alignmentCompanion.x,1);
[~,fDA]=nanmax(nanmean(outSU.alignmentCompanion.y,1));
fDA_time=outSU.alignmentCompanion.x(fDA);
[~,fInSuTimes]=nanmin(abs(fDA_time-outSU.times));

% get change after DA event
changeWindow1=[0 0.25];
changeWindow2=[0.25 1];
temp=diff(outSU.times);
temp=temp(~isnan(temp));
timestep=mode(temp);
changeWindow1_inds=[ceil(changeWindow1(1)/timestep) ceil(changeWindow1(2)/timestep)];
changeWindow2_inds=[ceil(changeWindow2(1)/timestep) ceil(changeWindow2(2)/timestep)];
changeIn1=nan(1,length(outSU.su));
changeIn2=nan(1,length(outSU.su));
after1=nan(1,length(outSU.su));
after2=nan(1,length(outSU.su));
before1=nan(1,length(outSU.su));
before2=nan(1,length(outSU.su));
for i=1:length(outSU.su)
    temp=outSU.su(i).alignedData;
    after1(i)=nanmean(nanmean(temp(:,fInSuTimes+changeWindow1_inds(1)-1:fInSuTimes+changeWindow1_inds(2)-1),1),2);
    after2(i)=nanmean(nanmean(temp(:,fInSuTimes+changeWindow2_inds(1)-1:fInSuTimes+changeWindow2_inds(2)-1),1),2);
    before1(i)=nanmean(nanmean(temp(:,1:fInSuTimes+changeWindow1_inds(1)-changeWindow1_inds(2)),1),2); % same
    before2(i)=nanmean(nanmean(temp(:,1:fInSuTimes+changeWindow1_inds(1)-changeWindow1_inds(2)),1),2);
%     around1(i)=nanmean(nanmean(temp(:,fInSuTimes-changeWindow1_inds(2):fInSuTimes-1),1),2); % same
%     around2(i)=nanmean(nanmean(temp(:,fInSuTimes-changeWindow1_inds(2):fInSuTimes-1),1),2);
    around1(i)=nanmean(nanmean(temp(:,fInSuTimes-changeWindow1_inds(2):fInSuTimes-1+changeWindow1_inds(2)),1),2); % same
    around2(i)=nanmean(nanmean(temp(:,fInSuTimes-changeWindow1_inds(2):fInSuTimes-1+changeWindow1_inds(2)),1),2);
    changeIn1(i)=after1(i)-before1(i);
    changeIn2(i)=after2(i)-before2(i);
end

% plot scatter
% figure(); 
% scatter(changeIn1,changeIn2);
% 
% figure();
% scatter(changeIn1./before1,changeIn2./before2);

out.changeIn1=changeIn1;
out.changeIn2=changeIn2;
out.after1=after1;
out.after2=after2;
out.before1=before1;
out.before2=before2;
out.around1=around1;
out.around2=around2;

end

function phys_beh_tbt=putPhotoBehFieldIntoPhysBeh(photo_beh_tbt,phys_beh_tbt,putInField)

if photo_beh_tbt.this_is_which_beh==1
    ref=photo_beh_tbt.reference_into_beh2trialinds;
else
    ref=photo_beh_tbt.reference_into_beh1trialinds;
end
phys_beh_tbt_temp=zeros(size(phys_beh_tbt.times));
temp=photo_beh_tbt.(putInField);
for i=1:size(temp,1)
    whichrow=ref(i,1);
    if isnan(whichrow)
        continue
    end
    if whichrow>size(phys_beh_tbt_temp,1)
        continue
    end
    phys_beh_tbt_temp(whichrow,:)=temp(i,:);
end
phys_beh_tbt.(putInField)=phys_beh_tbt_temp;

end

function beh_tbt=makeSameFieldInBeh(photo_tbt,beh_tbt,usePhotoField,phototimes,behtimes)

temp=photo_tbt.(usePhotoField);
tempbehfield=zeros(size(behtimes));
for i=1:size(temp,1)
    currevents=temp(i,:);
    currtimes=phototimes(i,:);
    currbehtimes=behtimes(i,:);
    % find corresponding times in beh_tbt
    f=find(currevents>0.5);
    for j=1:length(f)
        [~,mi]=nanmin(abs(currtimes(f(j))-currbehtimes));
        if ~isnan(mi)
            tempbehfield(i,mi)=1;
        end
    end
end
beh_tbt.(usePhotoField)=tempbehfield;

end

function data=smoothPhotometry(data,Fs,lowPassCutoff,smoothFields)

disp('Smoothing photometry fields');
for i=1:length(smoothFields)
    currField=smoothFields{i};
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
    for j=1:size(beforenans,1)
        beforenans(j,:)=smooth(beforenans(j,:),60);
    end
    % put nans back in
    beforenans(isnan(data.(currField)))=nan;
    data.(currField)=beforenans;
end
    
end

function physiology_tbt=makeUnitSubgroups(physiology_tbt,downSamp)

grp1_fr_thresh=2; % in Hz
grp1_fr_contingency='<=';
grp2_fr_thresh=5; % in Hz
grp2_fr_contingency='>=';
grp1_halfwidth_thresh=2*10^-4; % in sec
grp1_halfwidth_contingency='>';
grp2_halfwidth_thresh=2*10^-4; % in sec
grp2_halfwidth_contingency='<=';

% get unit FRs (firing rates)
unitnames=physiology_tbt.details.unitnames;
fr=nan(1,length(unitnames));
for i=1:length(unitnames)
    fr(i)=nanmean(nanmean(physiology_tbt.(unitnames{i})));
end

% classify units as grp1, grp2 or other (grp3)
whichgrp=nan(1,length(unitnames));
for i=1:length(unitnames)
    if eval(['fr(i)' grp1_fr_contingency num2str(grp1_fr_thresh)]) && eval(['physiology_tbt.details.halfWidths(i)' grp1_halfwidth_contingency num2str(grp1_halfwidth_thresh)])
        whichgrp(i)=1;
    elseif eval(['fr(i)' grp2_fr_contingency num2str(grp2_fr_thresh)]) && eval(['physiology_tbt.details.halfWidths(i)' grp2_halfwidth_contingency num2str(grp2_halfwidth_thresh)])
        whichgrp(i)=2;
    else
        whichgrp(i)=3;
    end
end
physiology_tbt.details.whichgrp=whichgrp;

for i=1:3 % 3 groups
    physiology_tbt.(['grp' num2str(i) '_unitav'])=zeros(size(physiology_tbt.unitsum));
    for j=1:length(unitnames)
        if whichgrp(j)==i
            tmp=cat(3,physiology_tbt.(['grp' num2str(i) '_unitav']),physiology_tbt.(unitnames{j})); 
            physiology_tbt.(['grp' num2str(i) '_unitav'])=nansum(tmp,3);
        end
    end
    physiology_tbt.(['grp' num2str(i) '_unitav'])=physiology_tbt.(['grp' num2str(i) '_unitav'])./nansum(whichgrp==i);
end

if downSamp~=1
    for i=1:3
        temp=downSampMatrix(physiology_tbt.(['grp' num2str(i) '_unitav']),downSamp);
        reinterp_temp=interp1(downSampAv(nanmean(physiology_tbt.unitTimes,1),downSamp)',temp',nanmean(physiology_tbt.unitTimes,1)');
        physiology_tbt.(['grp' num2str(i) '_unitav'])=reinterp_temp';
    end
end

end

function plotSpikeWvfmsByFR(all_wvfms,all_halfWidths,unitnames,physiology_tbt)

whichcolormap='jet';

fr=nan(1,length(unitnames));
for i=1:length(unitnames)
    fr(i)=nanmean(nanmean(physiology_tbt.(unitnames{i})));
end

cmap=colormap(whichcolormap);
steps=floor(size(cmap,1)/length(unitnames));
indIntoCmap=1:steps:size(cmap,1);
indIntoCmap=indIntoCmap(1:length(unitnames));
figure();
for i=1:length(unitnames)
    plot(all_wvfms(i,:),'Color',cmap(indIntoCmap(i),:));
    hold on;
end

figure();
xoffset=0;
for i=1:length(unitnames)
    plot([1:size(all_wvfms,2)]++xoffset,all_wvfms(i,:),'Color',cmap(indIntoCmap(i),:));
    strtoplot=unitnames{i};
    f_=regexp(strtoplot,'_');
    strtoplot(f_)=' ';
    text(xoffset,0,strtoplot); 
    xoffset=xoffset+size(all_wvfms,2);
    hold on;
end

figure();
scatter(all_halfWidths,fr,[],cmap(indIntoCmap,:));
set(gca, 'YScale', 'log');
for i=1:length(all_halfWidths)
    strtoplot=unitnames{i};
    f_=regexp(strtoplot,'_');
    strtoplot(f_)=' ';
    text(all_halfWidths(i),fr(i),strtoplot); 
end

end

function outSU=makeSummaryFig(beh_fields,photo_fields,phys_fields,alignments,physthenphoto_fields,withintimewindow,timewindows,isPhysField,xranges,photometry_tbt,photo_beh_tbt,physiology_tbt,phys_beh_tbt,normalizeSU,maxTrialLength)

outSU=[];

% make figures
% Set up figure layout 1
mainfig1=figure();
Nh=length(beh_fields)+length(photo_fields)*2+length(phys_fields); % number of rows
Nw=length(alignments); % number of columns
gap=[.01 .03]; % between plots
marg_h=[.1 .01]; % margin
marg_w=[.1 .01];% marg_w=[.01 .01]; % margin
[ha,pos]=tight_subplot(Nh,Nw,gap,marg_h,marg_w);
% reorder ha so populates down rows within column first
reorder=[];
for i=0:Nw-1
    reorder=[reorder i+[1:Nw:Nw*Nh]];
end
ha=ha(reorder);
indintoha=1;
for i=1:Nw
    % each alignment
    for j=1:length(physthenphoto_fields)
        % each type of data to plot
        for k=1:length(beh_fields) % get behavior fields aligned
            if j==1
                if isPhysField(j)==1
                    if k==1
                        [fout,dataout,n_events_in_av,alignmentCompanion,f_heatmap,plotBehFieldOut,su]=plotPhysiologyResult(physiology_tbt,phys_beh_tbt,[],alignments{i},physthenphoto_fields{j},beh_fields{k},withintimewindow{i},timewindows{i},{nan ha(indintoha)});
                        plottedIntoWhichAxes=indintoha;
                        if ~strcmp(physthenphoto_fields{j},'unit_by_unit')
                            strtoplot=[physthenphoto_fields{j} ' aligned to ' alignments{i}];
                            f_=regexp(strtoplot,'_');
                            strtoplot(f_)=' ';
                            text(nanmin(xranges{i}),nanmax(nanmean(dataout.y,1)),strtoplot);
                        end
                        indintoha=indintoha+1;
                        closeAllBut(mainfig1);
                    else
                        [fout,dataout,n_events_in_av,alignmentCompanion,f_heatmap,plotBehFieldOut,su]=plotPhysiologyResult(physiology_tbt,phys_beh_tbt,[],alignments{i},physthenphoto_fields{j},beh_fields{k},withintimewindow{i},timewindows{i},{nan nan});
                        closeAllBut(mainfig1);
                    end
                    axes(ha(indintoha));
                    indintoha=indintoha+1;
                    plotWStderr(plotBehFieldOut.y,plotBehFieldOut.x,'k',[],size(plotBehFieldOut.y,1));
                    hold on;
                    rescale=nanmax(nanmean(plotBehFieldOut.y,1))/nanmax(nanmean(alignmentCompanion.y,1));
                    plotWStderr(alignmentCompanion.y.*rescale,alignmentCompanion.x,'g',[],size(alignmentCompanion.y,1));
                    strtoplot=[beh_fields{k} ' aligned to ' alignments{i}];
                    f_=regexp(strtoplot,'_');
                    strtoplot(f_)=' ';
                    text(nanmin(xranges{i}),nanmax(nanmean(plotBehFieldOut.y,1)),strtoplot);
                else
                    if k==1
                        [fout,dataout,n_events_in_av,alignmentCompanion,f_heatmap,plotBehFieldOut]=plotPhotometryResult(photometry_tbt,photo_beh_tbt,[],alignments{i},physthenphoto_fields{j},beh_fields{k},withintimewindow{i},timewindows{i},{ha(indintoha) ha(indintoha+1)});
                        text(nanmin(xranges{i}),nanmax(dataout.y),[physthenphoto_fields{j} ' aligned to ' alignments{i}]);
                        indintoha=indintoha+2;
                        closeAllBut(mainfig1);
                    else
                        [fout,dataout,n_events_in_av,alignmentCompanion,f_heatmap,plotBehFieldOut]=plotPhotometryResult(photometry_tbt,photo_beh_tbt,[],alignments{i},physthenphoto_fields{j},beh_fields{k},withintimewindow{i},timewindows{i},{nan nan});
                        closeAllBut(mainfig1);
                    end
                    axes(ha(indintoha));
                    indintoha=indintoha+1;
                    plotWStderr(plotBehFieldOut.y,plotBehFieldOut.x,'k',[],size(plotBehFieldOut.y,1));
                    hold on;
                    rescale=nanmax(nanmean(plotBehFieldOut.y,1))/nanmax(nanmean(alignmentCompanion.y,1));
                    plotWStderr(alignmentCompanion.y.*rescale,alignmentCompanion.x,'g',[],size(alignmentCompanion.y,1));
                    strtoplot=[beh_fields{k} ' aligned to ' alignments{i}];
                    f_=regexp(strtoplot,'_');
                    strtoplot(f_)=' ';
                    text(nanmin(xranges{i}),nanmax(nanmean(plotBehFieldOut.y,1)),strtoplot);
                end
            else
                if isPhysField(j)==1
                    if k==1
                        [fout,dataout,n_events_in_av,alignmentCompanion,f_heatmap,plotBehFieldOut,su]=plotPhysiologyResult(physiology_tbt,phys_beh_tbt,[],alignments{i},physthenphoto_fields{j},[],withintimewindow{i},timewindows{i},{nan ha(indintoha)});
                        plottedIntoWhichAxes=indintoha;
                        if ~strcmp(physthenphoto_fields{j},'unit_by_unit')
                            strtoplot=[physthenphoto_fields{j} ' aligned to ' alignments{i}];
                            f_=regexp(strtoplot,'_');
                            strtoplot(f_)=' ';
                            text(nanmin(xranges{i}),nanmax(nanmean(dataout.y,1)),strtoplot);
                        end
                        indintoha=indintoha+1;
                        closeAllBut(mainfig1);
                    else
                        continue
                    end
                else
                    if k==1
                        [fout,dataout,n_events_in_av,alignmentCompanion,f_heatmap,plotBehFieldOut]=plotPhotometryResult(photometry_tbt,photo_beh_tbt,[],alignments{i},physthenphoto_fields{j},[],withintimewindow{i},timewindows{i},{ha(indintoha) ha(indintoha+1)});
                        strtoplot=[physthenphoto_fields{j} ' aligned to ' alignments{i}];
                        f_=regexp(strtoplot,'_');
                        strtoplot(f_)=' ';
                        text(nanmin(xranges{i}),nanmax(nanmean(dataout.y,1)),strtoplot);
                        indintoha=indintoha+2;
                        closeAllBut(mainfig1);
                    else
                        continue
                    end
                end
            end
            if strcmp(physthenphoto_fields{j},'unit_by_unit') && k==1
                % replace with unit by unit plot
                axes(ha(plottedIntoWhichAxes));
                cla(ha(plottedIntoWhichAxes));
                if normalizeSU==true
                    alignmentCompanion.y=alignmentCompanion.y-nanmin(nanmean(alignmentCompanion.y,1));
                    alignmentCompanion.y=length(su).*alignmentCompanion.y./nanmax(nanmean(alignmentCompanion.y,1));
                end
                plotWStderr(alignmentCompanion.y,alignmentCompanion.x,'g',[],size(alignmentCompanion.y,1));
                hold on;
                unitByUnitPlot(su,dataout.x,maxTrialLength,normalizeSU);
                strtoplot=[physthenphoto_fields{j} ' aligned to ' alignments{i}];
                f_=regexp(strtoplot,'_');
                strtoplot(f_)=' ';
                text(nanmin(xranges{i}),nanmax(nanmean(alignmentCompanion.y,1)),strtoplot);
                if strcmp(alignments{i},'isPhotoEvent')
                    outSU.su=su;
                    outSU.times=dataout.x;
                    outSU.alignmentCompanion.y=alignmentCompanion.y;
                    outSU.alignmentCompanion.x=alignmentCompanion.x;
                end
            end
        end
    end
end

k=1;
for i=1:Nw
    % each alignment
    for j=1:Nh
        set(ha(k),'XLim',xranges{i});
        % spawn individual figures
        f=figure();
        newax=copyobj(ha(k),f);
        set(newax,'Position',[0.1 0.1 0.85 0.85]);
        k=k+1;
    end
end

end

function unitByUnitPlot(su,su_times,maxTrialLength,normalizeSU)

ds=5;
spaceBetween=0;
offset=0;
endPlotAtInd=find(downSampAv(su_times,ds)>maxTrialLength,1,'first')-1;
for i=1:length(su)
    temp=downSampMatrix(su(i).alignedData,ds);
    if normalizeSU==true
        % normalize to y range 0 to 1
        sumin=nanmin(nanmean(temp(:,1:endPlotAtInd),1));
        temp=temp-sumin;
        sumax=nanmax(nanmean(temp(:,1:endPlotAtInd),1));
        temp=temp./sumax;
    end
    thisismax=plotWStderr(temp,downSampAv(su_times,ds),'k',endPlotAtInd,size(temp,1),offset);
    hold on;
    if isnan(thisismax)
        continue
    else
        offset=thisismax+spaceBetween;
    end
end
ylim([0 offset]);

end

function thisismax=plotWStderr(varargin)

dataMatrix=varargin{1};
times=varargin{2};
c=varargin{3};
plotUntilInd=varargin{4};
nEvents=varargin{5};
if length(varargin)>5
    offset=varargin{6};
else
    offset=0;
end

showStdevInstead=false;

if isempty(plotUntilInd)
    plotUntilInd=size(dataMatrix,2);
end

plot(times(1:plotUntilInd),offset+nanmean(dataMatrix(:,1:plotUntilInd),1),'Color',c,'LineWidth',1);
hold on;
if showStdevInstead==true
    plot(times(1:plotUntilInd),offset+nanmean(dataMatrix(:,1:plotUntilInd),1)-nanstd(dataMatrix(:,1:plotUntilInd),[],1),'Color',c,'LineWidth',0.5);
    plot(times(1:plotUntilInd),offset+nanmean(dataMatrix(:,1:plotUntilInd),1)+nanstd(dataMatrix(:,1:plotUntilInd),[],1),'Color',c,'LineWidth',0.5);
    thisismax=offset+nanmax(nanmean(dataMatrix(:,1:plotUntilInd),1)+nanstd(dataMatrix(:,1:plotUntilInd),[],1));
else
    plot(times(1:plotUntilInd),offset+nanmean(dataMatrix(:,1:plotUntilInd),1)-nanstd(dataMatrix(:,1:plotUntilInd),[],1)./sqrt(nEvents),'Color',c,'LineWidth',0.5);
    plot(times(1:plotUntilInd),offset+nanmean(dataMatrix(:,1:plotUntilInd),1)+nanstd(dataMatrix(:,1:plotUntilInd),[],1)./sqrt(nEvents),'Color',c,'LineWidth',0.5);
    thisismax=offset+nanmax(nanmean(dataMatrix(:,1:plotUntilInd),1)+nanstd(dataMatrix(:,1:plotUntilInd),[],1)./sqrt(nEvents));
end
    
end

function closeAllBut(f)

figs2keep = f;
all_figs = findobj(0, 'type', 'figure');
delete(setdiff(all_figs, figs2keep));

end