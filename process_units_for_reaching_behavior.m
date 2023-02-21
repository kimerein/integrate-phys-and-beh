% function process_units_for_reaching_behavior()
%% for processing unit alignments to behavior

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
dataTable='C:\Users\sabatini\Downloads\Spike sorting analysis - Combined phys and photo.csv';
data_loc_array=table2cell(readtable(dataTable,'Format','%s%s%s%u%s%s%s%s%s%u%u%s%u%s%s','Delimiter', ',', 'HeaderLines', 0, 'ReadVariableNames', true));

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
% data_loc_array{i,13}=false; % flag if penetration was not through str
% data_loc_array{i,14}='green_ch'; % which photo ch to use for photometry
% data_loc_array{i,15}='Z:\MICROSCOPE\Kim\WHISPER recs\dLight1\20210311\photometry aligned to beh'; % where to save beh-aligned photometry data

%% 2. Auto-populate SU QC and behavior alignments -- for 84 sessions, running all of this takes about 12 hrs
% Discarding trials where unit dead or moved away
errorInDirectories={};
for i=1:size(data_loc_array,1)
    if strcmp(data_loc_array{i,6},'no_tbt')
        continue
    end
    % load tbt (trial by trial) data
    [physiology_tbt,beh2_tbt,behavior_tbt,photometry_tbt]=loadReachingExptPhysData(data_loc_array(i,:));
    % save photometry alignments to behavior, if photometry was acquired
    if ~isempty(photometry_tbt) && ~strcmp(data_loc_array{i,14},'photo_not_working') && ~strcmp(data_loc_array{i,14},'not_finished')
        mkdir(data_loc_array{i,15});
        saveBehaviorAlignmentsPhotometry(photometry_tbt,behavior_tbt,data_loc_array{i,15},'',data_loc_array{i,14}); 
        disp(['Photometry for ' data_loc_array{i,6} ': Done']);
    end
    if strcmp(data_loc_array{i,7},'no_spikes')
        % no spikes recorded
        continue
    end
    if data_loc_array{i,13}==1
        % recording is not in structure
        disp(['Skipping  ' data_loc_array{i,6} ' because not in structure']);
        continue
    end
    % get spikes if recorded physiology
    dd=dir(data_loc_array{i,7});
    disp(['Processing ' data_loc_array{i,7}]);
    %try
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
                    doBecauseMissing=false;
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
    disp(['Physiology for ' data_loc_array{i,7} ': Done']);
%     catch
%         disp(['caught error while processing ' data_loc_array{i,6}]);
%         errorInDirectories{length(errorInDirectories)+1}=data_loc_array{i,6};
%     end
end

%% 3. Load data locations
dd=cell(1,size(data_loc_array,1));
for i=1:size(data_loc_array,1)
    % load locations of SU data aligned to behavior
    % e.g., 'Z:\MICROSCOPE\Kim\WHISPER recs\Mar_3\20210721\SU aligned to behavior';
    dd{i}=data_loc_array{i,8};
end

%% 3.5 For Mac, change pointers to files on server
% /Volumes/Neurobio/MICROSCOPE/Kim/WHISPER recs
if ismac()
    for i=1:length(dd)
        temp=dd{i};
        temp=['/Volumes/Neurobio/' temp(4:end)];
        f=regexp(temp,'\');
        temp(f)=sep;
        dd{i}=temp;
    end
end

%% 4.0 Get photometry
dd_photo=cell(1,size(data_loc_array,1));
for i=1:size(data_loc_array,1)
    % load locations of photo data aligned to behavior
    % e.g., 'Z:\MICROSCOPE\Kim\WHISPER recs\Mar_3\20210721\photometry aligned to behavior';
    dd_photo{i}=data_loc_array{i,15};
end
% choose type of response to plot
response_to_plot='uncued_pelletMissing'; % can be any of the directories created in saveBehaviorAlignmentsSingleNeuron.m
dd_photo_more=cell(1,length(dd_photo)); 
for i=1:length(dd_photo)
    dd_photo_more{i}=[dd_photo{i} sep response_to_plot];
end
Response_photo=getAndSaveResponse(dd_photo_more,'_',settingsForStriatumPhotometryPlots,[]);

%% 4.02 Get behavior d-primes for these days
settingsForBeh=reachExpt_analysis_settings('display settings');
% override some of these settings
settingsForBeh.check_for_human=0; settingsForBeh.tryForFiles={};
Response_beh=getAndSaveBeh(data_loc_array(:,6),settingsForBeh);

%% 4.03 Plot some photometry results
sub=subResponse(Response_photo,'fromWhichSess',find(Response_beh.dprimes<0.5));
plotVariousSUResponsesAlignedToBeh('meanAcrossUnits',sub,1);
% Get photo peaks
photolocs=getPhotoPeaks(sub,'dip',[0 3]);
figure(); histogram(photolocs,25);

%% 4. Make figures -- about 6 min to load 84 sessions of unit data
% choose type of response to plot
response_to_plot='cued_pelletMissing'; % can be any of the directories created in saveBehaviorAlignmentsSingleNeuron.m

% doUnitTest.m is used to test whether to include unit in this plot
% will include unit if unitdets match the following
% [inStructure isFS isTAN isSPN isLowFRThin]
plotUnitCriteria=[1 0 0 1 0]; % -100 is a wildcard, else 0 (false) and 1 (true)
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
trial_n_cutoff=6;
Response=excludeTooFewTrials(Response,trial_n_cutoff,true);

% plot some stuff
plotVariousSUResponsesAlignedToBeh('meanAcrossUnits',Response,1); % down sample factor is last arg
out=plotVariousSUResponsesAlignedToBeh('scatterInTimeWindows',Response,[-3 0.25],[0.25 4]); % time windows relative to alignment companion peak
out=plotVariousSUResponsesAlignedToBeh('modulationIndex',Response,[-3 0.25],[0.25 4]); % time windows relative to alignment companion peak
plotVariousSUResponsesAlignedToBeh('ZscoredAndSmoothed',Response,[-3 0.25],[0.25 4]); % time windows relative to alignment companion peak
plotVariousSUResponsesAlignedToBeh('populationVectorCorrelation',Response,0.25,[0.25 4]); % time bins step, then slice at time window
plotVariousSUResponsesAlignedToBeh('trialVectorCorrelation',Response,0.25,[0.25 4]); % time bins step, then slice at time window

% get sig responses
out=plotVariousSUResponsesAlignedToBeh('getResponsive',Response,[-3 0.25],[0.5 2.5]);
Response.isSig=out.isResponsive.isSig';
Response.pvals=out.isResponsive.pvals';
Response.positiveMod=out.isResponsive.positiveMod';
Response.sustained=out.isResponsive.sustained';

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
plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',cued_success,Response,'coloredPSTH',[-3 -0.5],[0 3],'ColorLabel',cue_Response);
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

% plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',cuedSuccess_Response,uncued_success_Response,'scatterInTimeWindows',[-3 0.25],[0.5 4],'ColorLabel',cue_Response,[-1 -0.5],[0 1]);
% plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',cue_Response,uncued_success_Response,'scatterInTimeWindows','timeWindowsForResponse1',[-1 -0.5],[0 1],'timeWindowsForResponse2',[-3 0.25],[0.5 4]);
plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',excludeTooFewTrials(cued_success_Response,trial_n_cutoff,false),excludeTooFewTrials(uncued_success_Response,trial_n_cutoff,false),'scatterInTimeWindows',[-3 -0.1],[0.5 4],'ColorLabel',excludeTooFewTrials(cue_Response,trial_n_cutoff,false),[-1 -0.5],[0 1]);
plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',excludeTooFewTrials(cued_failure_Response,trial_n_cutoff,false),excludeTooFewTrials(uncued_failure_Response,trial_n_cutoff,false),'scatterInTimeWindows',[-3 -0.1],[0.5 4],'ColorLabel',excludeTooFewTrials(cue_Response,trial_n_cutoff,false),[-1 -0.5],[0 1]);

% BEST SO FAR
trial_n_cutoff=3;
plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',excludeTooFewTrials(cued_success_Response,trial_n_cutoff,false),excludeTooFewTrials(uncued_success_Response,trial_n_cutoff,false),'scatterInTimeWindows',[-3 -0.1],[0.5 4],'ColorLabel',excludeTooFewTrials(cue_Response,trial_n_cutoff,false),[-1 -0.3],[-0.3 1]);
plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',excludeTooFewTrials(cued_failure_Response,trial_n_cutoff,false),excludeTooFewTrials(uncued_failure_Response,trial_n_cutoff,false),'scatterInTimeWindows',[-3 -0.1],[0.5 4],'ColorLabel',excludeTooFewTrials(cue_Response,trial_n_cutoff,false),[-1 -0.3],[-0.3 1]);

%% LDA analysis
whichSess=82; % 81 has a fair number of cells

% downSampBy=25; % downsamp 60 ms bins by this much
% takeNPointsAfterEvent=3;

downSampBy=4; % downsamp 60 ms bins by this much
takeNPointsAfterEvent=15;
takeNPointsBeforeEvent=7;

response_to_plot1='cued_success';
response_to_plot2='uncued_success';
LDA_analysis(whichSess,downSampBy,takeNPointsAfterEvent,takeNPointsBeforeEvent,response_to_plot1,response_to_plot2,dd,'SVM');

%% Put together tensors
whichSess=[81 60 75 5 11 22 31 47 62 63 65 66 67 68 70 82];
whichSess=[81 75 11 22 47 62 63 65 67 82]; % with some cue responsive neuron(s), e.g., 47
whichSess=[81 75];
downSampBy=5; % downsamp 60 ms bins by this much
takeNPointsAfterEvent=12;
takeNPointsBeforeEvent=0;
tensor1=[]; allLabels1=[];
for i=1:length(whichSess)
    [tensor, allLabels, timepoints_for_tensor]=getTensorsForOneSession(whichSess(i), downSampBy, takeNPointsAfterEvent, takeNPointsBeforeEvent, dd);
    if ~isempty(tensor1)
        [tensor,allLabels]=catAndExpandTensor(tensor1, tensor, allLabels1, allLabels);
    end
    tensor1=tensor; allLabels1=allLabels;
    close all;
end
disp(size(tensor))
save(['C:\Users\sabatini\Documents\tensor.mat'],'tensor');
save(['C:\Users\sabatini\Documents\allLabels.mat'],'allLabels');
save(['C:\Users\sabatini\Documents\timepoints_for_tensor.mat'],'timepoints_for_tensor');

%% GLM analysis
for i=441:450
    if data_loc_array{i,13}==1
        continue
    end
    if strcmp(data_loc_array{i,8},'no_tbt')
        continue
    end
    GLM_analysis(i,data_loc_array,10);
    close all;
end

% Read in unit names
plotUnitCriteria=[1 0 0 1 0]; getCriteriaForUnitsToPlot(plotUnitCriteria); dd_more=cell(1,length(dd)); 
for i=1:length(dd)
    dd_more{i}=[dd{i} sep 'cue']; % just a placeholder to read in names
end
whichUnitsToGrab='_'; 
[unitbyunit_names.names,unitbyunit_names.excluded,unitbyunit_names.D1orD2taggingExpt,unitbyunit_names.D1taggedCells,unitbyunit_names.A2ataggedCells,unitbyunit_names.fromWhichSess,unitbyunit_names.fromWhichUnit]=alignToReadInUnitNames(dd_more,whichUnitsToGrab,settingsForStriatumUnitPlots);
% Read in GLM coefficients
for i=1:length(dd)
    dd_more{i}=[dd{i} sep 'matglm']; % just a placeholder to read in names
end
[all_glm_coef,unitnames_glm,fromWhichSess_glm]=readInGLMCoefs(dd_more,settingsForStriatumUnitPlots);
% Outliers for matlab glm
outli=[312 313 319 320 322 324 374 376 377 378 379 380 381 382 383 386 387 388 389];
all_glm_coef=all_glm_coef(~ismember(1:size(all_glm_coef,1),outli),:); unitnames_glm=unitnames_glm(~ismember(1:length(unitnames_glm),outli)); fromWhichSess_glm=fromWhichSess_glm(~ismember(1:length(fromWhichSess_glm),outli));
indexGLMcellsIntoUnitNames=getNamesIndexIntoNamesList(unitnames_glm,unitbyunit_names);
load('Z:\MICROSCOPE\Kim\WHISPER recs\Mar_6\20210625\SU aligned to behavior\forglm\output\neuron1_glm.mat')
[ts,allco]=plotGLMcoef(all_glm_coef,0,feature_names,10*0.01,9,'mean'); title('python glm');
load('Z:\MICROSCOPE\Kim\WHISPER recs\Mar_6\20210625\SU aligned to behavior\matglm\features_for_glm.mat');
[ts,allco]=plotGLMcoef(all_glm_coef,[],fnames,10*0.01,9,'mean'); title('mat glm');
whichCoefToUse=[4 5 6]; studyGLMcoef(all_glm_coef,ts,whichCoefToUse);

%% Get significant responses 
% Get significance from trial by trial
a=load('Z:\MICROSCOPE\Kim\20221107 figures for Sys Club talk\allunits_cueResponse_trialbytrial.mat'); allTrials_cue_Response=a.allTrials_cue_Response;
a=load('Z:\MICROSCOPE\Kim\20221107 figures for Sys Club talk\allunits_cuedSuccess_trialbytrial.mat'); allTrials_cuedSuccess=a.Response;
a=load('Z:\MICROSCOPE\Kim\20221107 figures for Sys Club talk\allunits_cuedFailure_trialbytrial.mat'); allTrials_cuedFailure=a.Response;
a=load('Z:\MICROSCOPE\Kim\20221107 figures for Sys Club talk\allunits_uncuedSuccess_trialbytrial.mat'); allTrials_uncuedSuccess=a.Response;
a=load('Z:\MICROSCOPE\Kim\20221107 figures for Sys Club talk\allunits_uncuedFailure_trialbytrial.mat'); allTrials_uncuedFailure=a.Response;
responsesForSig={'allTrials_cue_Response','allTrials_cuedSuccess','allTrials_cuedFailure','allTrials_uncuedSuccess','allTrials_uncuedFailure'}; % responses aligned to various behavior events
responseTimeWindows={[-1 -0.3 -0.3 1],[-3 0 1 2],[-3 0 2 3],[-3 0 1 2],[-3 0 2 3]};
out=plotVariousSUResponsesAlignedToBeh('matchUnitsAcrossResponses',allTrials_cue_Response,allTrials_cuedSuccess,allTrials_cuedFailure,allTrials_uncuedSuccess,allTrials_uncuedFailure);
allTrials_cue_Response=out.Response1; allTrials_cuedSuccess=out.Response2; allTrials_cuedFailure=out.Response3; allTrials_uncuedSuccess=out.Response4; allTrials_uncuedFailure=out.Response5;
for i=1:length(responsesForSig) % get sig responses
    temp=responseTimeWindows{i};
    out=plotVariousSUResponsesAlignedToBeh('getResponsive',eval(responsesForSig{i}),temp(1:2),temp(3:4));
    if i==1 % init
        isSig=nan(length(out.isResponsive.isSig),length(responsesForSig));
        pvals=nan(length(out.isResponsive.isSig),length(responsesForSig));
        positiveMod=nan(length(out.isResponsive.isSig),length(responsesForSig));
        inWindow1Mean=nan(length(out.isResponsive.isSig),length(responsesForSig));
        inWindow2Mean=nan(length(out.isResponsive.isSig),length(responsesForSig));
        allWindowMean=nan(length(out.isResponsive.isSig),length(responsesForSig));
        allWindowVar=nan(length(out.isResponsive.isSig),length(responsesForSig));
    end
    isSig(:,i)=out.isResponsive.isSig; pvals(:,i)=out.isResponsive.pvals; positiveMod(:,i)=out.isResponsive.positiveMod; inWindow1Mean(:,i)=out.isResponsive.inWindow1Mean; inWindow2Mean(:,i)=out.isResponsive.inWindow2Mean;
    allWindowMean(:,i)=out.isResponsive.allWindowMean; allWindowVar(:,i)=out.isResponsive.allWindowVar;
end
anyIsSig=any(isSig==1,2);

%% Set up data matrix
% Units X conditions (alignments to beh events) X time
a=load('Z:\MICROSCOPE\Kim\20230205 all SU alignments\all trials averaged not downsampled\cue.mat'); cue_Response=a.Response; 
a=load('Z:\MICROSCOPE\Kim\20230205 all SU alignments\all trials averaged not downsampled\cued_success.mat'); cued_success_Response=a.Response;  
a=load('Z:\MICROSCOPE\Kim\20230205 all SU alignments\all trials averaged not downsampled\cued_failure.mat'); cued_failure_Response=a.Response; 
a=load('Z:\MICROSCOPE\Kim\20230205 all SU alignments\all trials averaged not downsampled\uncued_success.mat'); uncued_success_Response=a.Response; 
a=load('Z:\MICROSCOPE\Kim\20230205 all SU alignments\all trials averaged not downsampled\uncued_failure.mat'); uncued_failure_Response=a.Response;
% a=load('Z:\MICROSCOPE\Kim\20230205 all SU alignments\all trials averaged not downsampled\cue_noReach.mat'); cue_noReach_Response=a.Response;
%a=load('Z:\MICROSCOPE\Kim\20221129 lab meeting\responses unit by unit\uncued_reach.mat'); uncued_reach_Response=a.Response;
trial_n_cutoff=0;
out=plotVariousSUResponsesAlignedToBeh('matchUnitsAcrossResponses',excludeTooFewTrials(cue_Response,trial_n_cutoff,false),excludeTooFewTrials(cued_success_Response,trial_n_cutoff,false),excludeTooFewTrials(cued_failure_Response,trial_n_cutoff,false),excludeTooFewTrials(uncued_success_Response,trial_n_cutoff,false),excludeTooFewTrials(uncued_failure_Response,trial_n_cutoff,false));
cue_Response=out.Response1; cued_success_Response=out.Response2; cued_failure_Response=out.Response3; uncued_success_Response=out.Response4; uncued_failure_Response=out.Response5;
% out=plotVariousSUResponsesAlignedToBeh('matchUnitsAcrossResponses',excludeTooFewTrials(cue_noReach_Response,trial_n_cutoff,false),excludeTooFewTrials(cued_success_Response,trial_n_cutoff,false),excludeTooFewTrials(cued_failure_Response,trial_n_cutoff,false),excludeTooFewTrials(uncued_success_Response,trial_n_cutoff,false),excludeTooFewTrials(uncued_failure_Response,trial_n_cutoff,false));
% cue_noReach_Response=out.Response1;  cued_success_Response=out.Response2; cued_failure_Response=out.Response3; uncued_success_Response=out.Response4; uncued_failure_Response=out.Response5;
D1tag=cued_success_Response.D1tag(cued_success_Response.excluded==0); A2atag=cued_success_Response.A2atag(cued_success_Response.excluded==0); 
save('Z:\MICROSCOPE\Kim\20221129 lab meeting\responses unit by unit\for_data_matrix_D1vA2a.mat','D1tag','A2atag');

% takePointsBeforeZero=5; %15;
% takePointsAfterZero=70;

takePointsBeforeZero=30; %15;
takePointsAfterZero=450;

dataMatrix=setUpDataMatrix(cue_Response,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,takePointsBeforeZero,takePointsAfterZero);
% dataMatrix=setUpDataMatrix(cue_noReach_Response,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,takePointsBeforeZero,takePointsAfterZero);
dataMatrix(dataMatrix<0)=0; % no firing rates below 0
% take just outcome alignments
% dataMatrix=dataMatrix(:,:,[2:5]);
% dataMatrix_cueNoReach=setUpDataMatrix(cue_noReach_Response,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,takePointsBeforeZero,takePointsAfterZero);
% dataMatrix_cueNoReach(dataMatrix_cueNoReach<0)=0; % no firing rates below 0

% load('Z:\MICROSCOPE\Kim\20221129 lab meeting\dataMatrix.mat')
% dataMatrix=dataMatrix(:,6:end-35,:);
% dataMatrix(:,1:end-6,3)=dataMatrix(:,7:end,3);
% dataMatrix(:,1:end-6,5)=dataMatrix(:,7:end,5);

% Peak dopamine at 0.83 sec
% Dopamine dip at 1.5685 sec
dataMatrix=dataMatrix(:,31:end,:);
dataMatrix(:,1:end-73,3)=dataMatrix(:,74:end,3);
dataMatrix(:,1:end-73,5)=dataMatrix(:,74:end,5);

% clear newDataMatrix
% for i=1:size(dataMatrix,3)
%     temp=reshape(dataMatrix(:,:,i),size(dataMatrix(:,:,i),1),size(dataMatrix(:,:,i),2));
%     %temp=downSampMatrix(reshape(dataMatrix(:,:,i),size(dataMatrix(:,:,i),1),size(dataMatrix(:,:,i),2)),10);
%     for j=1:size(temp,1)
%         temp(j,:)=smoothdata(temp(j,:),'gaussian',10); %42);
%     end
% %     newDataMatrix(:,:,i)=abs(diff(temp,1,2)); 
%     newDataMatrix(:,:,i)=temp; 
% end

% PCA, CCA, etc.
% CCA: Find orthogonal dimensions of max covariance between X=[] and Y=[]
plotN=6;
boot=1; % num iterations for bootstrap
principaledCA(dataMatrix,{'units','time','conditions'},plotN,boot);

% zscore_cueR=cue_Response.unitbyunit_y; zscore_cueR(zscore_cueR<0.0001)=0;
% zscore_cueR=(zscore_cueR)./repmat(std(zscore_cueR(:,500:1000),[],2,'omitnan'),1,size(zscore_cueR,2));
% figure(); plot(nanmean(zscore_cueR,1));
% cuedcellresponse=nanmean(zscore_cueR(:,286:328),2)-nanmean(zscore_cueR(:,210:282),2);
% cuedcellresponse(cuedcellresponse>35)=35;
% figure(); histogram(cuedcellresponse,500);
figure(); plot(nanmean(cued_failure_Response.unitbyunit_x,1),nanmean(cued_failure_Response.unitbyunit_y(idx==2 & cuedcellresponse>0.5,:),1),'Color','b');
hold on; plot(nanmean(cued_failure_Response.aligncomp_x,1),nanmean(cued_failure_Response.aligncomp_y(idx==2 & cuedcellresponse>0.5,:),1),'Color','c'); title('failure cued');
figure(); plot(nanmean(cued_success_Response.unitbyunit_x,1),nanmean(cued_success_Response.unitbyunit_y(idx==2 & cuedcellresponse>0.5,:),1),'Color','b');
hold on; plot(nanmean(cued_success_Response.aligncomp_x,1),nanmean(cued_success_Response.aligncomp_y(idx==2 & cuedcellresponse>0.5,:),1),'Color','c'); title('success cued');
figure(); plot(nanmean(cued_failure_Response.unitbyunit_x,1),nanmean(cued_failure_Response.unitbyunit_y(idx==2 & cuedcellresponse<0.5,:),1),'Color','b');
hold on; plot(nanmean(cued_failure_Response.aligncomp_x,1),nanmean(cued_failure_Response.aligncomp_y(idx==2 & cuedcellresponse<0.5,:),1),'Color','c'); title('failure medium');
figure(); plot(nanmean(cued_success_Response.unitbyunit_x,1),nanmean(cued_success_Response.unitbyunit_y(idx==2 & cuedcellresponse<0.5,:),1),'Color','b');
hold on; plot(nanmean(cued_success_Response.aligncomp_x,1),nanmean(cued_success_Response.aligncomp_y(idx==2 & cuedcellresponse<0.5,:),1),'Color','c'); title('success medium');
figure(); plot(nanmean(cued_failure_Response.unitbyunit_x,1),nanmean(cued_failure_Response.unitbyunit_y(idx==2 & cuedcellresponse<0,:),1),'Color','b');
hold on; plot(nanmean(cued_failure_Response.aligncomp_x,1),nanmean(cued_failure_Response.aligncomp_y(idx==2 & cuedcellresponse<0,:),1),'Color','c'); title('failure uncued');
figure(); plot(nanmean(cued_success_Response.unitbyunit_x,1),nanmean(cued_success_Response.unitbyunit_y(idx==2 & cuedcellresponse<0,:),1),'Color','b');
hold on; plot(nanmean(cued_success_Response.aligncomp_x,1),nanmean(cued_success_Response.aligncomp_y(idx==2 & cuedcellresponse<0,:),1),'Color','c'); title('success uncued');