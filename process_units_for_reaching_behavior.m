% function process_units_for_reaching_behavior()
% for processing unit alignments to behavior

clear all

onVPN=false; %true; % if want to skip reading in spikes from raw data
goodUnitLabel=2; 
% for discarding units dead or moved away
dsinds=225;
percentThresh=5;
timeStretchThresh=60*10; % in seconds
plotInference=false;
channelSpacing=20; % in microns
skipIfBehaviorAlignmentsAlreadyExist=true; % if true, will skip any units for which a behavior alignment already exists in the cue sub-folder
% CAN COMMENT OUT SOME BEHAVIOR ALIGNMENTS IN
% saveBehaviorAlignmentsSingleNeuron.m IF DON'T WANT TO REPOPULATE
% EVERYTHING
skipUnitDetails=false; % if true, will skip populating the unit details
skipUnitDetailsUnlessNoQCfig=true; % skips populating unit details unless QC fig does not yet exist
% BUT NOTE WILL STILL REPOPULATE UNIT DETAILS
% second row is Matlab index, first row is depth on probe, where 32 is most
% dorsall, 1 is most ventral
% for A1x32Edge
% was using 20 micron spacing between chs
chDepthMapping=[1   21; ...
                2   24; ...
                3   22; ...
                4   23; ...
                5   20; ...
                6   26; ...
                7   19; ...
                8   28; ...
                9   18; ...
                10  30; ...
                11  17; ...
                12  32; ...
                13  25; ...
                14  27; ...
                15  29; ...
                16  31; ...
                17  2; ...
                18  4; ...
                19  6; ...
                20  8; ...
                21  1; ...
                22  5; ...
                23  3; ...
                24  10; ...
                25  7; ...
                26  12; ...
                27  9; ...
                28  13; ...
                29  11; ...
                30  15; ...
                31  14; ...
                32  16];

%% 1. Get locations of data: aligned phys events and spikes, raw phys data, probe location and depth, etc.
% Cell array mapping
% Date, mouse, raw data physiology file name and directory, probe depth,
% recording location, location of trial by trial data, spikes location,
% location of SU aligned to behavior, raw data binary name

% ultimately, ingest this from .csv file
dataTable='C:\Users\sabatini\Downloads\Spike sorting analysis - data as of 20221105.csv';
data_loc_array=table2cell(readtable(dataTable,'Format','%s%s%s%u%s%s%s%s%s%u%u%s'));

% data_loc_array=cell(2,6);
% i=1; data_loc_array{i,1}='20210311'; data_loc_array{i,2}='Oct_B'; 
% data_loc_array{i,3}='Z:\MICROSCOPE\Kim\WHISPER recs\Oct_B\20210311\phys\OctB2__g0'; % raw data phys location, or empty
% data_loc_array{i,4}=2300; % ventral depth of recording probe tip in microns relative to bregma
% data_loc_array{i,5}='pDMS-tail';
% data_loc_array{i,6}='Z:\MICROSCOPE\Kim\WHISPER recs\Oct_B\20210311\tbt';
% data_loc_array{i,7}='Z:\MICROSCOPE\Kim\WHISPER recs\Oct_B\20210311\spike output';
% data_loc_array{i,8}='Z:\MICROSCOPE\Kim\WHISPER recs\Oct_B\20210311\SU aligned to behavior';
% data_loc_array{i,9}='OctB2__g0_t0.nidq.bin';
% data_loc_array{i,10}=2000; % dorsal-most depth that is posterior striatum
% data_loc_array{i,11}=2600; % ventral-most depth that is posterior striatum
% data_loc_array{i,12}='Nkx2.1'; % if opto-tagging, which tag
% i=2; data_loc_array{i,1}='20210311'; data_loc_array{i,2}='dLight1'; 
% data_loc_array{i,3}='Z:\MICROSCOPE\Kim\WHISPER recs\dLight1\20210311\phys\dLight1__g0'; % raw data phys location, or empty
% data_loc_array{i,4}=2300; % ventral depth in microns relative to bregma
% data_loc_array{i,5}='pDMS-tail';
% data_loc_array{i,6}='Z:\MICROSCOPE\Kim\WHISPER recs\dLight1\20210311\tbt';
% data_loc_array{i,7}='Z:\MICROSCOPE\Kim\WHISPER recs\dLight1\20210311\phys\dLight1__g0';
% data_loc_array{i,8}='Z:\MICROSCOPE\Kim\WHISPER recs\dLight1\20210311\SU aligned to behavior';
% data_loc_array{i,9}='dLight1__g0_t0.nidq.bin';
% data_loc_array{i,10}=2000; % dorsal-most depth that is posterior striatum
% data_loc_array{i,11}=2600; % ventral-most depth that is posterior striatum
% data_loc_array{i,12}='none'; % if opto-tagging, which tag

%% 2. Auto-populate SU QC and behavior alignments -- for 84 sessions, running all of this takes about 12 hrs
% Discarding trials where unit dead or moved away
for i=1:size(data_loc_array,1)
    % load tbt (trial by trial) data
    [physiology_tbt,beh2_tbt,behavior_tbt,photometry_tbt]=loadReachingExptPhysData(data_loc_array(i,:));
    % get spikes
    dd=dir(data_loc_array{i,7});
    disp(['Processing ' data_loc_array{i,7}]);
%     try
    for j=1:length(dd)
        if ~isempty(regexp(dd(j).name,'spikes'))
            % load spikes
            a=load([data_loc_array{i,7} sep dd(j).name]);
            if ~isfield(a,'spikes')
                disp(['cannot find spikes in ' data_loc_array{i,7} sep dd(j).name]);
                continue
            end
            spikes=a.spikes;
            switch dd(j).name
                case 'spikes.mat'
                    trodeChsForSpikes=[1 2 3 4];
                case 'spikes2.mat'
                    trodeChsForSpikes=[5 6 7 8];
                case 'spikes3.mat'
                    trodeChsForSpikes=[9 10 11];
                case 'spikes4.mat'
                    trodeChsForSpikes=[13 14];
                case 'spikes5.mat'
                    trodeChsForSpikes=[15 16 17 18];
                case 'spikes6.mat'
                    trodeChsForSpikes=[19 20];
                case 'spikes7.mat'
                    trodeChsForSpikes=[21 22 23 24];
                case 'spikes8.mat'
                    trodeChsForSpikes=[23 25 26 27];
                case 'spikes9.mat'
                    trodeChsForSpikes=[28 29 30 31];
                case 'spikes_sorted.mat'
                    trodeChsForSpikes=[1 2 3 4];
                case 'spikes2_sorted.mat'
                    trodeChsForSpikes=[5 6 7 8];
                case 'spikes3_sorted.mat'
                    trodeChsForSpikes=[9 10 11];
                case 'spikes4_sorted.mat'
                    trodeChsForSpikes=[13 14];
                case 'spikes5_sorted.mat'
                    trodeChsForSpikes=[15 16 17 18];
                case 'spikes6_sorted.mat'
                    trodeChsForSpikes=[19 20];
                case 'spikes7_sorted.mat'
                    trodeChsForSpikes=[21 22 23 24];
                case 'spikes8_sorted.mat'
                    trodeChsForSpikes=[23 25 26 27];
                case 'spikes9_sorted.mat'
                    trodeChsForSpikes=[28 29 30 31];
                otherwise 
                    error('do not recognize name of spikes mat file');
            end
            % find good units
            gu=find(spikes.labels(:,2)==goodUnitLabel);
            if skipUnitDetailsUnlessNoQCfig==true
                % are all SU_QC figs already populated? if yes, continue
                allsue=zeros(1,length(gu));
                for k=1:length(gu)
                    currU=gu(k);
                    currAssign=spikes.labels(currU,1);
                    allsue(k)=SU_QC_file_exists(data_loc_array{i,8}, currAssign, trodeChsForSpikes);
                end
                if all(allsue==1)
                    disp(['Already made SU QC figs for this spikes']);
                else
                    doBecauseMissing=true;
                end
            else
                doBecauseMissing=true;
            end
            % make opto-tagged alignment
            if skipUnitDetails==false || doBecauseMissing==true
                [~,tbtspikes]=organizeSpikesToMatch_physiology_tbt(spikes,physiology_tbt);
            end
            % if no good units, just continue
            if isempty(gu)
                continue
            end
            if skipUnitDetails==false || doBecauseMissing==true
                optoAligned_phys_tbt=alignToOpto(addGoodUnitsAsFields(physiology_tbt,tbtspikes,2,1,false,true,true));
                optoAligned_phys_tbt=checkWaveformsDuringOpto(optoAligned_phys_tbt,tbtspikes);
                % save opto alignment
                if ~exist([data_loc_array{i,8} sep 'opto_aligned'],'dir')
                    mkdir([data_loc_array{i,8} sep 'opto_aligned']);
                end
                save([data_loc_array{i,8} sep 'opto_aligned' sep 'phys_tbt_for_' dd(j).name],'optoAligned_phys_tbt');
                clear tbtspikes
            end
            % unit by unit
            % and other alignments
            for k=1:length(gu)
                currU=gu(k);
                currAssign=spikes.labels(currU,1);
                % find ch where waveform biggest
                amp=nan(1,size(spikes.waveforms,3));
                for l=1:size(spikes.waveforms,3)
                    amp(l)=abs(min(reshape(mean(spikes.waveforms(spikes.assigns==currAssign,:,l),1,'omitnan'),1,size(spikes.waveforms,2)),[],2,'omitnan'));
                end
                [~,si]=sort(amp);
                % sort trode chs for units
                if any(si>length(trodeChsForSpikes))
                    disp(['problem indexing spike channel for unit ' num2str(currAssign) ' in ' data_loc_array{i,7} sep dd(j).name]);
                    si=si(si<=length(trodeChsForSpikes));
                end
                trodeChsForSpikes=trodeChsForSpikes(si);
                forplot_trodeChs=trodeChsForSpikes;
                if length(forplot_trodeChs)<4
                    forplot_trodeChs=[forplot_trodeChs ones(1,4-length(forplot_trodeChs))*forplot_trodeChs(1)];
                end
                % 
                if skipIfBehaviorAlignmentsAlreadyExist==true
                    if cue_alignment_file_exists([data_loc_array{i,8} sep 'cue'], currAssign, trodeChsForSpikes(end))
                        skipBehAlign=true;
                    else
                        skipBehAlign=false;
                    end
                else
                    skipBehAlign=false;
                end
                if skipUnitDetails==false
                    addTag=isThisUnitOptoTagged(optoAligned_phys_tbt,currAssign,data_loc_array{i,12});
                end
                % check whether QC figure already exists for this unit
                [sue,qc_fname]=SU_QC_file_exists(data_loc_array{i,8}, currAssign, trodeChsForSpikes(end));
                if skipUnitDetails==true
                    % need to get addTag from SU_QC fig
                    if sue==true
                        r=regexp(qc_fname,'_');
                        addTag=qc_fname(r(1)+1:r(2)-1);
                    else
                        error('Cannot skip unit details if SU_QC file does not yet exist');
                    end
                end
                if skipBehAlign==false || sue==false
                    if ~sue && onVPN==false
                        % make behavior alignments and single unit QC figs
                        unit_data=saveBehaviorAlignmentsSingleNeuron(physiology_tbt,spikes,currAssign,beh2_tbt,forplot_trodeChs,trodeChsForSpikes(end),data_loc_array{i,8},addTag,data_loc_array{i,3},data_loc_array{i,9});
                    else
                        if sue==false
                            % make QC fig without reading in raw data
                            unit_data=saveBehaviorAlignmentsSingleNeuron(physiology_tbt,spikes,currAssign,beh2_tbt,forplot_trodeChs,trodeChsForSpikes(end),data_loc_array{i,8},addTag,data_loc_array{i,3},[]);
                        else
                            % just redo behavior alignments, skipping QC fig
                            spikes.skipQC=true;
                            saveBehaviorAlignmentsSingleNeuron(physiology_tbt,spikes,currAssign,beh2_tbt,forplot_trodeChs,trodeChsForSpikes(end),data_loc_array{i,8},addTag,data_loc_array{i,3},[]);
                            unit_data=[data_loc_array{i,8} sep qc_fname];
                        end
                    end
                else
                    unit_data=[data_loc_array{i,8} sep qc_fname];
                end
                dud=true;
                if skipUnitDetailsUnlessNoQCfig==true
                    if sue==true
                        dud=false;
                    else
                        dud=true;
                    end
                end
                if skipUnitDetails==true
                    dud=false;
                end
                if dud==true
                    % discard trials where unit dead or moved away
                    [physiology_tbt.(['unit' num2str(currAssign) '_dontUseTrials']),physiology_tbt.(['unit' num2str(currAssign) '_meanFR'])]=inferUnitStability(unit_data,physiology_tbt,dsinds,percentThresh,timeStretchThresh,plotInference);
                    % get depth of unit
                    physiology_tbt.(['unit' num2str(currAssign) '_depth'])=getUnitDepth_forWHISPER(trodeChsForSpikes(end), chDepthMapping, channelSpacing, data_loc_array{i,4});
                    if physiology_tbt.(['unit' num2str(currAssign) '_depth'])>=data_loc_array{i,10} && physiology_tbt.(['unit' num2str(currAssign) '_depth'])<=data_loc_array{i,11}
                        physiology_tbt.(['unit' num2str(currAssign) '_inStructure'])=true;
                    else
                        physiology_tbt.(['unit' num2str(currAssign) '_inStructure'])=false;
                    end
                    % get waveform features
                    [physiology_tbt.(['unit' num2str(currAssign) '_halfwidth']), physiology_tbt.(['unit' num2str(currAssign) '_peakToTrough']), ~, physiology_tbt.(['unit' num2str(currAssign) '_avWaveforms'])]=getSUWaveformFeatures(filtspikes_without_sweeps(spikes,0,'assigns',currAssign),spikes.params.Fs);
                    physiology_tbt.(['unit' num2str(currAssign) '_amp'])=amp(si(end));
                    % use unit details to classify as type of striatum unit
                    [physiology_tbt.(['unit' num2str(currAssign) '_isFS']),physiology_tbt.(['unit' num2str(currAssign) '_isTAN']),physiology_tbt.(['unit' num2str(currAssign) '_isSPN']),physiology_tbt.(['unit' num2str(currAssign) '_isLowFRThin'])]=classifyStriatumUnits(physiology_tbt.(['unit' num2str(currAssign) '_halfwidth']),physiology_tbt.(['unit' num2str(currAssign) '_meanFR']));
                    % save unit details for this unit only
                    % will be used later to group units
                    if ~exist([data_loc_array{i,8} sep 'unit_details'],'dir')
                        mkdir([data_loc_array{i,8} sep 'unit_details']);
                    end
                    saveUnitDets([data_loc_array{i,8} sep 'unit_details'],physiology_tbt,currAssign,trodeChsForSpikes(end));
                    saveOptoTagDets([data_loc_array{i,8} sep 'opto_aligned'],optoAligned_phys_tbt,currAssign,trodeChsForSpikes(end));
                end
            end
            % save unit details for these spikes
            if skipUnitDetails==false
                if ~exist([data_loc_array{i,8} sep 'unit_details'],'dir')
                    mkdir([data_loc_array{i,8} sep 'unit_details']);
                end
                save([data_loc_array{i,8} sep 'unit_details' sep 'phys_tbt_for_' dd(j).name],'physiology_tbt');
            end
        end
    end
%     catch
%         disp('caught error');
%     end
end

%% 3. Load data locations
dd=cell(1,size(data_loc_array,1));
for i=1:size(data_loc_array,1)
    % load locations of SU data aligned to behavior
    % e.g., 'Z:\MICROSCOPE\Kim\WHISPER recs\Mar_3\20210721\SU aligned to behavior';
    dd{i}=data_loc_array{i,8};
end

%% 4. Make figures -- about 6 min to load 84 sessions of unit data
% choose type of response to plot
response_to_plot='uncued_failure'; % can be any of the directories created in saveBehaviorAlignmentsSingleNeuron.m

% doUnitTest.m is used to test whether to include unit in this plot
% will include unit if unitdets match the following
% [inStructure isFS isTAN isSPN isLowFRThin]
plotUnitCriteria=[-100 0 0 1 0]; % -100 is a wildcard, else 0 (false) and 1 (true)
getCriteriaForUnitsToPlot(plotUnitCriteria);
% read in some units
dd_more=cell(1,length(dd)); 
for i=1:length(dd)
    dd_more{i}=[dd{i} sep response_to_plot];
end
whichUnitsToGrab='_'; % '_' for all units, or can be something like 'D1tagged'
Response=getAndSaveResponse(dd_more,whichUnitsToGrab,settingsForStriatumUnitPlots,[]);
% for opto-tagging comparing tagged v untagged
% [tagged_Response,D1orD2taggingExpt,putAlignPeakAt]=getAndSaveResponse(dd_more,'D1tagged',settingsForStriatumUnitPlots,[]);
% [untagged_Response]=getAndSaveResponse(dd_more,'__',settingsForStriatumUnitPlots,putAlignPeakAt);
% untagged_Response=takeOnlyUnitsFromSess(untagged_Response,unique(tagged_Response.fromWhichSess));

% exclude units with too few trials
trial_n_cutoff=3;
Response=excludeTooFewTrials(Response,trial_n_cutoff,true);

% plot some stuff
plotVariousSUResponsesAlignedToBeh('meanAcrossUnits',Response,1); % down sample factor is last arg
out=plotVariousSUResponsesAlignedToBeh('scatterInTimeWindows',Response,[-3 0.25],[0.25 4]); % time windows relative to alignment companion peak
out=plotVariousSUResponsesAlignedToBeh('modulationIndex',Response,[-3 0.25],[0.25 4]); % time windows relative to alignment companion peak
plotVariousSUResponsesAlignedToBeh('ZscoredAndSmoothed',Response,[-3 0.25],[0.25 4]); % time windows relative to alignment companion peak
plotVariousSUResponsesAlignedToBeh('populationVectorCorrelation',Response,0.25,[0.25 4]); % time bins step, then slice at time window
plotVariousSUResponsesAlignedToBeh('trialVectorCorrelation',Response,0.25,[0.25 4]); % time bins step, then slice at time window

% compare same unit responses to different events
% line up same units using excluded field
out=plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',Response1,Response2,'modulationIndex',[-3 0.25],[0.25 4],'ColorLabel',Response3,[-1 -0.5],[-0.4 1]); 
% 'ColorLabel' then Response3 must be after other args for plot type
% optionally can pass in different windows for Response3 calculation (or
% don't include those args to use same windows as for Response1 and
% Response2)
plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',cue_D1tagged,cued_success_D1tagged,'modulationIndex',[-1 -0.5],[0.25 1],'ColorLabel',cued_failure_D1tagged);
plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',Response1,Response2,'meanAcrossUnits',1);
plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',cued_success_D1tagged,cued_failure_D1tagged,'scatterInTimeWindows',[-3 0.25],[0.25 4],'ColorLabel',cue_D1tagged,[-1 -0.5],[-0.4 1]);

out=plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',cued_success_Response,cued_failure_Response,'modulationIndex',[-3 0.25],[0.25 4],'ColorLabel',uncued_success_Response);
failure_minus_success=out.response2minus1;
out=plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',cued_success_Response,uncued_success_Response,'modulationIndex',[-3 0.25],[0.25 4],'ColorLabel',cued_failure_Response);
uncued_minus_cued=out.response2minus1;

out=plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',cued_success_D1tagged,cued_failure_D1tagged,'modulationIndex',[-1 -0.5],[0.25 1],'ColorLabel',uncued_success_D1tagged);
failure_minus_success=out.response2minus1;
out=plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',cued_success_D1tagged,uncued_success_D1tagged,'modulationIndex',[-1 -0.5],[0.25 1],'ColorLabel',cued_failure_D1tagged);
uncued_minus_cued=out.response2minus1;
figure(); scatter(uncued_minus_cued,failure_minus_success);
line([0 0],[-2 2]); line([-2 2],[0 0]);
xlabel('uncued minus cued'); ylabel('failure minus success');
[r,p]=corrcoef(uncued_minus_cued(~isnan(uncued_minus_cued) & ~isnan(failure_minus_success)),failure_minus_success(~isnan(uncued_minus_cued) & ~isnan(failure_minus_success)));
out=plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',cued_success_A2atagged,cued_failure_A2atagged,'modulationIndex',[-1 -0.5],[0.25 1],'ColorLabel',uncued_success_A2atagged);
failure_minus_success=out.response2minus1;
out=plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',cued_success_A2atagged,uncued_success_A2atagged,'modulationIndex',[-1 -0.5],[0.25 1],'ColorLabel',cued_failure_A2atagged);
uncued_minus_cued=out.response2minus1;
figure(); scatter(uncued_minus_cued,failure_minus_success);
line([0 0],[-2 2]); line([-2 2],[0 0]);
xlabel('uncued minus cued'); ylabel('failure minus success');
[r,p]=corrcoef(uncued_minus_cued(~isnan(uncued_minus_cued) & ~isnan(failure_minus_success)),failure_minus_success(~isnan(uncued_minus_cued) & ~isnan(failure_minus_success)));
title('A2a');

out=plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',cued_success_D1tagged,cued_failure_D1tagged,'modulationIndex',[-3 -0.5],[0.25 4],'ColorLabel',uncued_success_D1tagged);
failure_minus_success=out.response2minus1;
out=plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',cued_success_D1tagged,uncued_success_D1tagged,'modulationIndex',[-3 -0.5],[0.25 4],'ColorLabel',cued_failure_D1tagged);
uncued_minus_cued=out.response2minus1;
figure(); scatter(uncued_minus_cued,failure_minus_success);
line([0 0],[-2 2]); line([-2 2],[0 0]);
xlabel('uncued minus cued'); ylabel('failure minus success'); title('D1');
out=plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',cued_success_A2atagged,cued_failure_A2atagged,'modulationIndex',[-3 -0.5],[0.25 4],'ColorLabel',uncued_success_A2atagged);
failure_minus_success=out.response2minus1;
out=plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',cued_success_A2atagged,uncued_success_A2atagged,'modulationIndex',[-3 -0.5],[0.25 4],'ColorLabel',cued_failure_A2atagged);
uncued_minus_cued=out.response2minus1;
figure(); scatter(uncued_minus_cued,failure_minus_success);
line([0 0],[-2 2]); line([-2 2],[0 0]);
xlabel('uncued minus cued'); ylabel('failure minus success'); title('A2a');

plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',cue_D1tagged,cued_success_D1tagged,'modulationIndex',[-1 -0.5],[-0.4 1],'ColorLabel',cued_failure_D1tagged);
plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',cue_A2atagged,cued_success_A2atagged,'modulationIndex',[-1 -0.5],[-0.4 1],'ColorLabel',cued_failure_A2atagged);

% scriptToOrganizeD1vD2unitResponses is old and probably won't work
% getCriteriaForUnitsToPlot('alltrue');
% scriptToOrganizeD1vD2unitResponses_wrapper(dd);

plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',cued_success_D1tagged,cued_failure_D1tagged,'modulationIndex',[-3 0.25],[0.5 4],'ColorLabel',cue_D1tagged,[-1 -0.5],[0 0.25]);

plotVariousSUResponsesAlignedToBeh('modRatioPSTH',cued_success_Response,cued_failure_Response,cueNoReach_Response,[-3 0.25],[0.5 2.5]);

% D1mods=plotVariousSUResponsesAlignedToBeh('scatter3D',cued_success_D1tagged,cued_failure_D1tagged,cueNoReach_D1tagged,[-3 0.25],[0.5 4]);
% D1mods_forCue=plotVariousSUResponsesAlignedToBeh('scatter3D',cued_success_D1tagged,cued_failure_D1tagged,cueNoReach_D1tagged,[-1 -0.5],[-0.4 1]);
% D1mods.modIndex3=D1mods_forCue.modIndex3;
% % D1meanFR=D1mods_forCue.meanFR3;
% A2amods=plotVariousSUResponsesAlignedToBeh('scatter3D',cued_success_A2atagged,cued_failure_A2atagged,cueNoReach_A2atagged,[-3 0.25],[0.5 4]);
% A2amods_forCue=plotVariousSUResponsesAlignedToBeh('scatter3D',cued_success_A2atagged,cued_failure_A2atagged,cueNoReach_A2atagged,[-1 -0.5],[-0.4 1]);
% A2amods.modIndex3=A2amods_forCue.modIndex3;
% % A2ameanFR=A2amods_forCue.meanFR3;
% SVMforD1vA2a(D1mods,A2amods);
trial_n_cutoff=0;
D1mods=plotVariousSUResponsesAlignedToBeh('scatter3D',excludeTooFewTrials(cued_success_D1tagged,trial_n_cutoff,true),excludeTooFewTrials(cued_failure_D1tagged,trial_n_cutoff,true),excludeTooFewTrials(cueNoReach_D1tagged,trial_n_cutoff,true),[-3 0.25],[0.5 2.5]);
D1mods_forCue=plotVariousSUResponsesAlignedToBeh('scatter3D',excludeTooFewTrials(cued_success_D1tagged,trial_n_cutoff,true),excludeTooFewTrials(cued_failure_D1tagged,trial_n_cutoff,true),excludeTooFewTrials(cueNoReach_D1tagged,trial_n_cutoff,true),[-1 -0.5],[-0.4 1]);
D1mods.modIndex3=D1mods_forCue.modIndex3;
% D1meanFR=D1mods_forCue.meanFR3;
A2amods=plotVariousSUResponsesAlignedToBeh('scatter3D',excludeTooFewTrials(cued_success_A2atagged,trial_n_cutoff,true),excludeTooFewTrials(cued_failure_A2atagged,trial_n_cutoff,true),excludeTooFewTrials(cueNoReach_A2atagged,trial_n_cutoff,true),[-3 0.25],[0.5 2.5]);
A2amods_forCue=plotVariousSUResponsesAlignedToBeh('scatter3D',excludeTooFewTrials(cued_success_A2atagged,trial_n_cutoff,true),excludeTooFewTrials(cued_failure_A2atagged,trial_n_cutoff,true),excludeTooFewTrials(cueNoReach_A2atagged,trial_n_cutoff,true),[-1 -0.5],[-0.4 1]);
A2amods.modIndex3=A2amods_forCue.modIndex3;
% A2ameanFR=A2amods_forCue.meanFR3;
SVMforD1vA2a(D1mods,A2amods);