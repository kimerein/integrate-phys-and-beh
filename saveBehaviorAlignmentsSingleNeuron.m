function unit_data=saveBehaviorAlignmentsSingleNeuron(phys_tbt,spikes,assign,behavior_tbt,trodeChs,unit_on_channel,saveDir,saveName,rawDataDir,rawDataBin)

binsize=10; % in ms
bsmooth=false;
unit_data=[];

% SU QC
if isfield(spikes, 'skipQC')
else
    spikes.skipQC=false; 
end
if spikes.skipQC==false
    if isempty(rawDataBin)
        varAdd.onVPN=true;
    else
        varAdd.onVPN=false;
    end
    varAdd.trodeChs=trodeChs;
    varAdd.furtherProcessData=@furtherProcessWHISPER;
    unit_data=SU_QC(spikes, assign, unit_on_channel, rawDataBin, rawDataDir, [], varAdd);
    saveas(gcf,[saveDir sep 'unit' num2str(assign) 'onCh' num2str(unit_on_channel) '_' saveName '_QC.fig']);
    % if varAdd.onVPN==true
    %     fid=fopen([saveDir '\VPN_too_slow_to_populate_raw_data.txt'],'wt');
    %     fclose(fid);
    % end
end

% Add unit
[~,spikes]=organizeSpikesToMatch_physiology_tbt(spikes,phys_tbt);

% Align to opto
optoAligned_phys_tbt=alignToOpto(addGoodUnitsAsFields(phys_tbt,spikes,2,1,false,true));
checkWaveformsDuringOpto(optoAligned_phys_tbt,spikes);

% Bin spikes into PSTHs for behavior alignments
phys_tbt=addGoodUnitsAsFields(phys_tbt,spikes,2,binsize,bsmooth,true);

% Save various alignments to behavior events
% Cue
event='cue';
timeWindow=[-1 16]; % in seconds from cue onset
subDir='cue';
[~,dataout,~,alignComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[]); close all;
saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'cueAligned',dataout,alignComp);

% Cue followed by success
event='cue_followedby_success';
timeWindow=[-1 16]; % in seconds from cue onset
subDir='cue_followedby_success';
[~,dataout,~,alignComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[]); close all;
saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'cueFollowedBySuccess',dataout,alignComp);

% Cued success
event='success_fromPerchOrWheel';
timeWindow=[0 3]; % in seconds from cue onset
subDir='cued_success';
[~,dataout,~,alignComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[]); close all;
saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'cuedSuccess',dataout,alignComp);

% Cued failure
event='misses_and_pelletMissing_and_drop';
timeWindow=[0 3]; % in seconds from cue onset
subDir='cued_failure';
[~,dataout,~,alignComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[]); close all;
saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'cuedFailure',dataout,alignComp);

% Uncued success
event='success_fromPerchOrWheel';
timeWindow=[5 16]; % in seconds from cue onset
subDir='uncued_success';
[~,dataout,~,alignComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[]); close all;
saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'uncuedSuccess',dataout,alignComp);

% Uncued failure
event='misses_and_pelletMissing_and_drop';
timeWindow=[5 16]; % in seconds from cue onset
subDir='uncued_failure';
[~,dataout,~,alignComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[]); close all;
saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'uncuedFailure',dataout,alignComp);

end

function saveStuff(saveDir,subDir,assign,saveName,alignmentName,dataout,alignComp)

if ~exist([saveDir sep subDir],'dir')
    mkdir([saveDir sep subDir]);
end
save([saveDir sep subDir sep 'unit' assign '_' saveName '_' alignmentName '.mat'],'dataout','alignComp');

end
