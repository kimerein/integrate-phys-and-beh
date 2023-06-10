function [data,timepoints]=reachAlignedPSTH(phys_tbt,behavior_tbt,spikes,assign,unitdets)

binsize=10; % in ms
bsmooth=false;
suppressFigs=true;
cueOffset=-0.16; % in sec
indsbefore=200;
indsafter=500;

% Get opto trial types
[wasOptoDuringCue,wasOptoBeforeCue,wasOptoAfterCue,noOpto,f,fBefore,fAfter,f_noOpto]=getTrialTypes_basedOnOpto(phys_tbt);

% Add unit
[~,spikes]=organizeSpikesToMatch_physiology_tbt(spikes,phys_tbt);
% Bin spikes into PSTHs for behavior alignments
phys_tbt=addGoodUnitsAsFields(phys_tbt,spikes,2,binsize,bsmooth,true,suppressFigs);

% opto during cue
behavior_tbt.reach_onTrial=behavior_tbt.all_reachBatch;
behavior_tbt.reach_onTrial(wasOptoDuringCue~=1,:)=0; % zero out non-matching trials
behavior_tbt.reach_onTrial(unitdets.dontUseTrials==1,:)=0; % zero out trials where unit missing
event='reach_onTrial';
timeWindow=[0.1 0.4]; % in seconds from cue onset
[~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
[temp1,temp2]=getSameRange(dataout,alignComp,indsbefore,indsafter);
data.optoduringcue=temp1; 
timepoints.optoduringcue=temp2;

% opto before cue
behavior_tbt.reach_onTrial=behavior_tbt.all_reachBatch;
behavior_tbt.reach_onTrial(wasOptoBeforeCue~=1,:)=0; % zero out non-matching trials
behavior_tbt.reach_onTrial(unitdets.dontUseTrials==1,:)=0; % zero out trials where unit missing
event='reach_onTrial';
timeWindow=[-0.4 -0.1]; % in seconds from cue onset
[~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
[temp1,temp2]=getSameRange(dataout,alignComp,indsbefore,indsafter);
data.optobeforecue=temp1; 
timepoints.optobeforecue=temp2;

% opto after cue
behavior_tbt.reach_onTrial=behavior_tbt.all_reachBatch;
behavior_tbt.reach_onTrial(wasOptoAfterCue~=1,:)=0; % zero out non-matching trials
behavior_tbt.reach_onTrial(unitdets.dontUseTrials==1,:)=0; % zero out trials where unit missing
event='reach_onTrial';
timeWindow=[0.4 0.7]; % in seconds from cue onset
[~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
[temp1,temp2]=getSameRange(dataout,alignComp,indsbefore,indsafter);
data.optoaftercue=temp1; 
timepoints.optoaftercue=temp2;

% all opto on together
data.allopto=mean([data.optobeforecue; data.optoduringcue; data.optoaftercue],1,'omitnan');
timepoints.allopto=mean([timepoints.optobeforecue; timepoints.optoduringcue; timepoints.optoaftercue],1,'omitnan');

% no opto
behavior_tbt.reach_onTrial=behavior_tbt.all_reachBatch;
behavior_tbt.reach_onTrial(noOpto~=1,:)=0; % zero out non-matching trials
behavior_tbt.reach_onTrial(unitdets.dontUseTrials==1,:)=0; % zero out trials where unit missing
event='reach_onTrial';
timeWindow=[-0.4 -0.1]; % in seconds from cue onset
[~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
[temp1,temp2]=getSameRange(dataout,alignComp,indsbefore,indsafter);
data.nooptobeforecue=temp1; 
timepoints.nooptobeforecue=temp2;
timeWindow=[0.1 0.4]; % in seconds from cue onset
[~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
[temp1,temp2]=getSameRange(dataout,alignComp,indsbefore,indsafter);
data.nooptoduringcue=temp1; 
timepoints.nooptoduringcue=temp2;
timeWindow=[0.4 0.7]; % in seconds from cue onset
[~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
[temp1,temp2]=getSameRange(dataout,alignComp,indsbefore,indsafter);
data.nooptoaftercue=temp1; 
timepoints.nooptoaftercue=temp2;

% all no opto together
data.allNOopto=mean([data.nooptobeforecue; data.nooptoduringcue; data.nooptoaftercue],1,'omitnan');
timepoints.allNOopto=mean([timepoints.nooptobeforecue; timepoints.nooptoduringcue; timepoints.nooptoaftercue],1,'omitnan');

figure();
ds=10;
plot(downSampAv(timepoints.allNOopto,ds),downSampAv(data.allNOopto,ds),'Color','k'); hold on;
plot(downSampAv(timepoints.allopto,ds),downSampAv(data.allopto,ds),'Color','r');

end

function [subdata,subtimes]=getSameRange(dataout,alignComp,indsbefore,indsafter)

timesalign=alignComp.x;
[~,indof]=nanmax(nanmean(alignComp.y,1));
timeof=timesalign(indof);
[~,mi]=nanmin(abs(nanmean(dataout.x,1)-timeof));
subdata=nanmean(dataout.y(:,mi-indsbefore:mi+indsafter),1);
subtimes=nanmean(dataout.x(:,mi-indsbefore:mi+indsafter),1)-timeof;

end

function [wasOptoDuringCue,wasOptoBeforeCue,wasOptoAfterCue,noOpto,f,fBefore,fAfter,f_noOpto]=getTrialTypes_basedOnOpto(physiology_tbt)

wasOptoDuringCue=nan(size(physiology_tbt.cue,1),1);
for i=1:size(physiology_tbt.cue,1)
    cueInds=find(physiology_tbt.cue(i,1:200)>0.5);
    optoOn=physiology_tbt.opto(i,cueInds);
    if nansum(optoOn)>3
        wasOptoDuringCue(i)=1;
    else
        wasOptoDuringCue(i)=0;
    end
end
f=find(wasOptoDuringCue==1);

wasOptoBeforeCue=nan(size(physiology_tbt.cue,1),1);
for i=1:size(physiology_tbt.cue,1)
    cueInds=find(physiology_tbt.cue(i,1:200)>0.5,1,'first');
    optoOn=physiology_tbt.opto(i,cueInds-10:cueInds-3);
    if nansum(optoOn)>3
        wasOptoBeforeCue(i)=1;
    else
        wasOptoBeforeCue(i)=0;
    end
end
fBefore=find(wasOptoBeforeCue==1);

wasOptoAfterCue=nan(size(physiology_tbt.cue,1),1);
for i=1:size(physiology_tbt.cue,1)
    cueInds=find(physiology_tbt.cue(i,1:200)>0.5,1,'first');
    optoOn=physiology_tbt.opto(i,cueInds+8:cueInds+15);
    if nansum(optoOn)>3
        wasOptoAfterCue(i)=1;
    else
        wasOptoAfterCue(i)=0;
    end
end
fAfter=find(wasOptoAfterCue==1);

noOpto=nan(size(physiology_tbt.cue,1),1);
for i=1:size(physiology_tbt.cue,1)
    cueInds=1:200;
    optoOn=physiology_tbt.opto(i,cueInds);
    if any(optoOn>0.5)
        noOpto(i)=0;
    else
        noOpto(i)=1;
    end
end
f_noOpto=find(noOpto==1);

end