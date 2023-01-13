function unit_data=saveBehaviorAlignmentsPhotometry(photo_tbt,behavior_tbt,saveDir,saveName)

% binsize=10; % in ms
% bsmooth=false;
% unit_data=[];
% suppressFigs=true;
% cueOffset=-0.16; % in sec
% 
% % SU QC
% if isfield(spikes, 'skipQC')
% else
%     spikes.skipQC=false; 
% end
% if spikes.skipQC==false
%     if isempty(rawDataBin)
%         varAdd.onVPN=true;
%     else
%         varAdd.onVPN=false;
%     end
%     varAdd.trodeChs=trodeChs;
%     varAdd.furtherProcessData=@furtherProcessWHISPER;
%     unit_data=SU_QC(spikes, assign, unit_on_channel, rawDataBin, rawDataDir, [], varAdd);
%     saveas(gcf,[saveDir sep 'unit' num2str(assign) 'onCh' num2str(unit_on_channel) '_' saveName '_QC.fig']);
%     % if varAdd.onVPN==true
%     %     fid=fopen([saveDir '\VPN_too_slow_to_populate_raw_data.txt'],'wt');
%     %     fclose(fid);
%     % end
% end
% 
% % Add unit
% [~,spikes]=organizeSpikesToMatch_physiology_tbt(spikes,photo_tbt);
% 
% % Bin spikes into PSTHs for behavior alignments
% photo_tbt=addGoodUnitsAsFields(photo_tbt,spikes,2,binsize,bsmooth,true,suppressFigs);
% 
% % Save various alignments to behavior events
% % Cue
% event='cue';
% timeWindow=[-1 16]; % in seconds from cue onset
% subDir='cue';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(photo_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'cueAligned',dataout,alignComp,phys_timepointsComp);
% 
% % Uncued reach
% event='all_reachBatch';
% timeWindow=[5 16]; % in seconds from cue onset
% subDir='uncued_reach';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(photo_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'uncuedReach',dataout,alignComp,phys_timepointsComp);
% 
% % Cue preceded by no reach and followed by no reach (for at least 4 sec)
% event='cue_noReach';
% timeWindow=[-1 16]; % in seconds from cue onset
% subDir='cue_noReach';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(photo_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'cueNoReach',dataout,alignComp,phys_timepointsComp);
% 
% % Cue followed by success
% event='cue_followedby_success';
% timeWindow=[-1 16]; % in seconds from cue onset
% subDir='cue_followedby_success';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(photo_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'cueFollowedBySuccess',dataout,alignComp,phys_timepointsComp);
% 
% % Cued success
% event='success_fromPerchOrWheel';
% timeWindow=[0+cueOffset 3]; % in seconds from cue onset
% subDir='cued_success';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(photo_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'cuedSuccess',dataout,alignComp,phys_timepointsComp);
% 
% % Cued failure
% event='misses_and_pelletMissing_and_drop';
% timeWindow=[0+cueOffset 3]; % in seconds from cue onset
% subDir='cued_failure';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(photo_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'cuedFailure',dataout,alignComp,phys_timepointsComp);
% 
% % Cued drop
% event='drop_fromPerchOrWheel';
% timeWindow=[0+cueOffset 3]; % in seconds from cue onset
% subDir='cued_drop';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(photo_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'cuedDrop',dataout,alignComp,phys_timepointsComp);
% 
% % Cued miss
% event='misses_and_pelletMissing';
% timeWindow=[0+cueOffset 3]; % in seconds from cue onset
% subDir='cued_miss';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(photo_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'cuedMiss',dataout,alignComp,phys_timepointsComp);
% 
% % Uncued success
% event='success_fromPerchOrWheel';
% timeWindow=[3 16]; % in seconds from cue onset
% subDir='uncued_success';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(photo_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'uncuedSuccess',dataout,alignComp,phys_timepointsComp);
% 
% % Uncued failure
% event='misses_and_pelletMissing_and_drop';
% timeWindow=[3 16]; % in seconds from cue onset
% subDir='uncued_failure';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(photo_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'uncuedFailure',dataout,alignComp,phys_timepointsComp);
% 
% % Uncued drop
% event='drop_fromPerchOrWheel';
% timeWindow=[3 16]; % in seconds from cue onset
% subDir='uncued_drop';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(photo_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'uncuedDrop',dataout,alignComp,phys_timepointsComp);
% 
% % Uncued miss
% event='misses_and_pelletMissing';
% timeWindow=[3 16]; % in seconds from cue onset
% subDir='uncued_miss';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(photo_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'uncuedMiss',dataout,alignComp,phys_timepointsComp);
% 
% % Cued failure and no reaching afterward
% event='failure_noSuccessBeforeAndNoReachingAfter';
% timeWindow=[0+cueOffset 3]; % in seconds from cue onset
% subDir='cued_failure_then_noReach';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(photo_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'cuedFailureThenNoReach',dataout,alignComp,phys_timepointsComp);
% 
% % Uncued failure and no reaching afterward
% event='failure_noSuccessBeforeAndNoReachingAfter';
% timeWindow=[3 16]; % in seconds from cue onset
% subDir='uncued_failure_then_noReach';
% [~,dataout,~,alignComp,~,~,~,phys_timepointsComp]=plotPhysiologyResult(photo_tbt,behavior_tbt,[],event,['unit' num2str(assign)],'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
% saveStuff(saveDir,subDir,[num2str(assign) 'onCh' num2str(unit_on_channel)],saveName,'uncuedFailureThenNoReach',dataout,alignComp,phys_timepointsComp);
% 
% end
% 
% function saveStuff(saveDir,subDir,assign,saveName,alignmentName,dataout,alignComp,phys_timepointsComp)
% 
% if ~exist([saveDir sep subDir],'dir')
%     mkdir([saveDir sep subDir]);
% end
% save([saveDir sep subDir sep 'unit' assign '_' saveName '_' alignmentName '.mat'],'dataout','alignComp','phys_timepointsComp');
% 
% end