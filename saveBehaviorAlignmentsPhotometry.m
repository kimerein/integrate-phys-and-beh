function saveBehaviorAlignmentsPhotometry(phys_tbt,behavior_tbt,saveDir,saveName,which_photoch)

suppressFigs=true;
cueOffset=-0.16; % in sec

% Save various alignments to behavior events
% Cue
event='cue';
timeWindow=[-1 16]; % in seconds from cue onset
subDir='cue';
[~,dataout,~,alignComp,~,~,phys_timepointsComp]=plotPhotometryResult(phys_tbt,behavior_tbt,[],event,which_photoch,'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
saveStuff(saveDir,subDir,'',saveName,'cueAligned',dataout,alignComp,phys_timepointsComp);

% Uncued reach
event='all_reachBatch';
timeWindow=[5 16]; % in seconds from cue onset
subDir='uncued_reach';
[~,dataout,~,alignComp,~,~,phys_timepointsComp]=plotPhotometryResult(phys_tbt,behavior_tbt,[],event,which_photoch,'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
saveStuff(saveDir,subDir,'',saveName,'uncuedReach',dataout,alignComp,phys_timepointsComp);

% Cue preceded by no reach and followed by no reach (for at least 4 sec)
event='cue_noReach';
timeWindow=[-1 16]; % in seconds from cue onset
subDir='cue_noReach';
[~,dataout,~,alignComp,~,~,phys_timepointsComp]=plotPhotometryResult(phys_tbt,behavior_tbt,[],event,which_photoch,'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
saveStuff(saveDir,subDir,'',saveName,'cueNoReach',dataout,alignComp,phys_timepointsComp);

% Cue followed by success
event='cue_followedby_success';
timeWindow=[-1 16]; % in seconds from cue onset
subDir='cue_followedby_success';
[~,dataout,~,alignComp,~,~,phys_timepointsComp]=plotPhotometryResult(phys_tbt,behavior_tbt,[],event,which_photoch,'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
saveStuff(saveDir,subDir,'',saveName,'cueFollowedBySuccess',dataout,alignComp,phys_timepointsComp);

% Cued success
event='success_fromPerchOrWheel';
timeWindow=[0+cueOffset 3]; % in seconds from cue onset
subDir='cued_success';
[~,dataout,~,alignComp,~,~,phys_timepointsComp]=plotPhotometryResult(phys_tbt,behavior_tbt,[],event,which_photoch,'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
saveStuff(saveDir,subDir,'',saveName,'cuedSuccess',dataout,alignComp,phys_timepointsComp);

% Cued failure
event='misses_and_pelletMissing_and_drop';
timeWindow=[0+cueOffset 3]; % in seconds from cue onset
subDir='cued_failure';
[~,dataout,~,alignComp,~,~,phys_timepointsComp]=plotPhotometryResult(phys_tbt,behavior_tbt,[],event,which_photoch,'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
saveStuff(saveDir,subDir,'',saveName,'cuedFailure',dataout,alignComp,phys_timepointsComp);

% Cued drop
event='drop_fromPerchOrWheel';
timeWindow=[0+cueOffset 3]; % in seconds from cue onset
subDir='cued_drop';
[~,dataout,~,alignComp,~,~,phys_timepointsComp]=plotPhotometryResult(phys_tbt,behavior_tbt,[],event,which_photoch,'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
saveStuff(saveDir,subDir,'',saveName,'cuedDrop',dataout,alignComp,phys_timepointsComp);

% Cued miss
event='misses_and_pelletMissing';
timeWindow=[0+cueOffset 3]; % in seconds from cue onset
subDir='cued_miss';
[~,dataout,~,alignComp,~,~,phys_timepointsComp]=plotPhotometryResult(phys_tbt,behavior_tbt,[],event,which_photoch,'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
saveStuff(saveDir,subDir,'',saveName,'cuedMiss',dataout,alignComp,phys_timepointsComp);

% Uncued success
event='success_fromPerchOrWheel';
timeWindow=[3 16]; % in seconds from cue onset
subDir='uncued_success';
[~,dataout,~,alignComp,~,~,phys_timepointsComp]=plotPhotometryResult(phys_tbt,behavior_tbt,[],event,which_photoch,'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
saveStuff(saveDir,subDir,'',saveName,'uncuedSuccess',dataout,alignComp,phys_timepointsComp);

% Uncued failure
event='misses_and_pelletMissing_and_drop';
timeWindow=[3 16]; % in seconds from cue onset
subDir='uncued_failure';
[~,dataout,~,alignComp,~,~,phys_timepointsComp]=plotPhotometryResult(phys_tbt,behavior_tbt,[],event,which_photoch,'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
saveStuff(saveDir,subDir,'',saveName,'uncuedFailure',dataout,alignComp,phys_timepointsComp);

% Uncued drop
event='drop_fromPerchOrWheel';
timeWindow=[3 16]; % in seconds from cue onset
subDir='uncued_drop';
[~,dataout,~,alignComp,~,~,phys_timepointsComp]=plotPhotometryResult(phys_tbt,behavior_tbt,[],event,which_photoch,'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
saveStuff(saveDir,subDir,'',saveName,'uncuedDrop',dataout,alignComp,phys_timepointsComp);

% Uncued miss
event='misses_and_pelletMissing';
timeWindow=[3 16]; % in seconds from cue onset
subDir='uncued_miss';
[~,dataout,~,alignComp,~,~,phys_timepointsComp]=plotPhotometryResult(phys_tbt,behavior_tbt,[],event,which_photoch,'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
saveStuff(saveDir,subDir,'',saveName,'uncuedMiss',dataout,alignComp,phys_timepointsComp);

% Cued failure and no reaching afterward
event='failure_noSuccessBeforeAndNoReachingAfter';
timeWindow=[0+cueOffset 3]; % in seconds from cue onset
subDir='cued_failure_then_noReach';
[~,dataout,~,alignComp,~,~,phys_timepointsComp]=plotPhotometryResult(phys_tbt,behavior_tbt,[],event,which_photoch,'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
saveStuff(saveDir,subDir,'',saveName,'cuedFailureThenNoReach',dataout,alignComp,phys_timepointsComp);

% Uncued failure and no reaching afterward
event='failure_noSuccessBeforeAndNoReachingAfter';
timeWindow=[3 16]; % in seconds from cue onset
subDir='uncued_failure_then_noReach';
[~,dataout,~,alignComp,~,~,phys_timepointsComp]=plotPhotometryResult(phys_tbt,behavior_tbt,[],event,which_photoch,'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
saveStuff(saveDir,subDir,'',saveName,'uncuedFailureThenNoReach',dataout,alignComp,phys_timepointsComp);

% Cued reach when pellet missing
event='pelletmissingreach_reachStarts';
timeWindow=[0+cueOffset 3]; % in seconds from cue onset
subDir='cued_pelletMissing';
[~,dataout,~,alignComp,~,~,phys_timepointsComp]=plotPhotometryResult(phys_tbt,behavior_tbt,[],event,which_photoch,'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
saveStuff(saveDir,subDir,'',saveName,'cuedReachPelletMissing',dataout,alignComp,phys_timepointsComp);

% Uncued reach when pellet missing
event='pelletmissingreach_reachStarts';
timeWindow=[3 16]; % in seconds from cue onset
subDir='uncued_pelletMissing';
[~,dataout,~,alignComp,~,~,phys_timepointsComp]=plotPhotometryResult(phys_tbt,behavior_tbt,[],event,which_photoch,'cueZone_onVoff','first',timeWindow,[],suppressFigs); close all;
saveStuff(saveDir,subDir,'',saveName,'uncuedReachPelletMissing',dataout,alignComp,phys_timepointsComp);

end

function saveStuff(saveDir,subDir,assign,saveName,alignmentName,dataout,alignComp,phys_timepointsComp)

if ~exist([saveDir sep subDir],'dir')
    mkdir([saveDir sep subDir]);
end
save([saveDir sep subDir sep 'ch' assign '_' saveName '_' alignmentName '.mat'],'dataout','alignComp','phys_timepointsComp');

end
