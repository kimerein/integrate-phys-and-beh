function unit_data=saveBehaviorAlignmentsSingleNeuron(phys_tbt,spikes,assign,behavior_tbt,trodeChs,unit_on_channel,saveDir,saveName,rawDataDir,rawDataBin)

binsize=10; % in ms
bsmooth=false;
unit_data=[];
suppressFigs=true;
cueOffset=-0.16; % in sec

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

% Bin spikes into PSTHs for behavior alignments
phys_tbt=addGoodUnitsAsFields(phys_tbt,spikes,2,binsize,bsmooth,true,suppressFigs);

% Save various alignments to behavior events
% % Cue
% event='cue';
% timeWindow=[-1 16]; % in seconds from cue onset
% subDir='cue';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'cueAligned',dataout,alignComp,phys_timepointsComp);
% 
% % Cued reach
% event='all_reachBatch';
% timeWindow=[0+cueOffset 3]; % in seconds from cue onset
% subDir='cued_reach';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'cuedReach',dataout,alignComp,phys_timepointsComp);
% 
% % Uncued reach
% event='all_reachBatch';
% timeWindow=[5 16]; % in seconds from cue onset
% subDir='uncued_reach';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'uncuedReach',dataout,alignComp,phys_timepointsComp);
% 
% % Cue preceded by no reach and followed by no reach (for at least 4 sec)
% event='cue_noReach';
% timeWindow=[-1 16]; % in seconds from cue onset
% subDir='cue_noReach';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'cueNoReach',dataout,alignComp,phys_timepointsComp);
% 
% % Cue followed by success
% event='cue_followedby_success';
% timeWindow=[-1 16]; % in seconds from cue onset
% subDir='cue_followedby_success';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'cueFollowedBySuccess',dataout,alignComp,phys_timepointsComp);
% 
% % Cued success
% event='success_fromPerchOrWheel';
% timeWindow=[0+cueOffset 3]; % in seconds from cue onset
% subDir='cued_success';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'cuedSuccess',dataout,alignComp,phys_timepointsComp);
% 
% % Cued failure
% event='misses_and_pelletMissing_and_drop';
% timeWindow=[0+cueOffset 3]; % in seconds from cue onset
% subDir='cued_failure';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'cuedFailure',dataout,alignComp,phys_timepointsComp);
% 
% % Cued drop
% event='drop_fromPerchOrWheel';
% timeWindow=[0+cueOffset 3]; % in seconds from cue onset
% subDir='cued_drop';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'cuedDrop',dataout,alignComp,phys_timepointsComp);
% 
% % Cued miss
% event='misses_and_pelletMissing';
% timeWindow=[0+cueOffset 3]; % in seconds from cue onset
% subDir='cued_miss';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'cuedMiss',dataout,alignComp,phys_timepointsComp);
% 
% % Uncued success
% event='success_fromPerchOrWheel';
% timeWindow=[3 16]; % in seconds from cue onset
% subDir='uncued_success';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'uncuedSuccess',dataout,alignComp,phys_timepointsComp);
% 
% % Uncued failure
% event='misses_and_pelletMissing_and_drop';
% timeWindow=[3 16]; % in seconds from cue onset
% subDir='uncued_failure';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'uncuedFailure',dataout,alignComp,phys_timepointsComp);
% 
% % Uncued drop
% event='drop_fromPerchOrWheel';
% timeWindow=[3 16]; % in seconds from cue onset
% subDir='uncued_drop';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'uncuedDrop',dataout,alignComp,phys_timepointsComp);
% 
% % Uncued miss
% event='misses_and_pelletMissing';
% timeWindow=[3 16]; % in seconds from cue onset
% subDir='uncued_miss';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'uncuedMiss',dataout,alignComp,phys_timepointsComp);
% 
% % Cued failure and no reaching afterward
% event='failure_noSuccessBeforeAndNoReachingAfter';
% timeWindow=[0+cueOffset 3]; % in seconds from cue onset
% subDir='cued_failure_then_noReach';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'cuedFailureThenNoReach',dataout,alignComp,phys_timepointsComp);
% 
% % Uncued failure and no reaching afterward
% event='failure_noSuccessBeforeAndNoReachingAfter';
% timeWindow=[3 16]; % in seconds from cue onset
% subDir='uncued_failure_then_noReach';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'uncuedFailureThenNoReach',dataout,alignComp,phys_timepointsComp);
% 
% % Cued reach when pellet missing
% event='pelletmissingreach_reachStarts';
% timeWindow=[0+cueOffset 3]; % in seconds from cue onset
% subDir='cued_pelletMissing';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'cuedReachPelletMissing',dataout,alignComp,phys_timepointsComp);
% 
% % Uncued reach when pellet missing
% event='pelletmissingreach_reachStarts';
% timeWindow=[3 16]; % in seconds from cue onset
% subDir='uncued_pelletMissing';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'uncuedReachPelletMissing',dataout,alignComp,phys_timepointsComp);
% 
% % Distractor
% event='movie_distractor';
% behavior_tbt.movie_distractor=behavior_tbt.movie_distractor>0.5;
% timeWindow=[-1 16]; % in seconds from cue onset
% subDir='distractor';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'distractorAligned',dataout,alignComp,phys_timepointsComp);
% 
% % All success
% event='success_fromPerchOrWheel';
% timeWindow=[-1 16]; % in seconds from cue onset
% subDir='all_success';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'allSuccess',dataout,alignComp,phys_timepointsComp);
% 
% % All failure
% event='misses_and_pelletMissing_and_drop';
% timeWindow=[-1 16]; % in seconds from cue onset
% subDir='all_failure';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'allFailure',dataout,alignComp,phys_timepointsComp);
% 
% % All drop
% event='drop_fromPerchOrWheel';
% timeWindow=[-1 16]; % in seconds from cue onset
% subDir='all_drop';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'allDrop',dataout,alignComp,phys_timepointsComp);
% 
% % All miss
% event='misses_and_pelletMissing';
% timeWindow=[-1 16]; % in seconds from cue onset
% subDir='all_miss';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'allMiss',dataout,alignComp,phys_timepointsComp);
% 
% % All reach when pellet missing
% event='pelletmissingreach_reachStarts';
% timeWindow=[-1 16]; % in seconds from cue onset
% subDir='all_pelletMissing';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'allReachPelletMissing',dataout,alignComp,phys_timepointsComp);

% Current cued success & cued 1 forward
event='success_to_reachinwindow';
timeWindow={[0+cueOffset 3],[0+cueOffset 3]}; % in seconds from cue onset, trial n then trial n+1
subDir='cued_success_to_cued';
[~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'cued_succ_to_cued',dataout,alignComp,phys_timepointsComp);

% Current cued success & uncued 1 forward
event='success_to_reachinwindow';
timeWindow={[0+cueOffset 3],[3 16]}; % in seconds from cue onset, trial n then trial n+1
subDir='cued_success_to_uncued';
[~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'cued_succ_to_uncued',dataout,alignComp,phys_timepointsComp);

% Current uncued success & cued 1 forward
event='success_to_reachinwindow';
timeWindow={[3 16],[0+cueOffset 3]}; % in seconds from cue onset, trial n then trial n+1
subDir='uncued_success_to_cued';
[~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'uncued_succ_to_cued',dataout,alignComp,phys_timepointsComp);

% Current uncued success & uncued 1 forward
event='success_to_reachinwindow';
timeWindow={[3 16],[3 16]}; % in seconds from cue onset, trial n then trial n+1
subDir='uncued_success_to_uncued';
[~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(phys_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'uncued_succ_to_uncued',dataout,alignComp,phys_timepointsComp);

end

function saveStuff(saveDir,subDir,assign,saveName,alignmentName,dataout,alignComp,phys_timepointsComp)

if ~exist([saveDir sep subDir],'dir')
    mkdir([saveDir sep subDir]);
end
save([saveDir sep subDir sep 'unit' assign '_' saveName '_' alignmentName '.mat'],'dataout','alignComp','phys_timepointsComp');

end
