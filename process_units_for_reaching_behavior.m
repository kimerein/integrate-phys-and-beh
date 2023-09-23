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
skipUnitDetails=true; % if true, will skip populating the unit details
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
response_to_plot='uncued_success'; % can be any of the directories created in saveBehaviorAlignmentsSingleNeuron.m
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
response_to_plot='cue'; % can be any of the directories created in saveBehaviorAlignmentsSingleNeuron.m

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
takeNPointsAfterEvent=17; 
takeNPointsBeforeEvent=0; %takeNPointsBeforeEvent=10; 
tensor1=[]; allLabels1=[];
whichSess=336; %335; 
for i=1:length(whichSess)
    [tensor, allLabels, timepoints_for_tensor]=getTensorsForOneSession(whichSess(i), downSampBy, takeNPointsAfterEvent, takeNPointsBeforeEvent, dd);
    if ~isempty(tensor1)
        [tensor,allLabels]=catAndExpandTensor(tensor1, tensor, allLabels1, allLabels);
    end
    tensor1=tensor; allLabels1=allLabels;
    close all;
end
tensor=tensor1; allLabels=allLabels1;
disp(size(tensor))
save(['C:\Users\sabatini\Documents\currtens\tensor.mat'],'tensor');
save(['C:\Users\sabatini\Documents\currtens\allLabels.mat'],'allLabels');
save(['C:\Users\sabatini\Documents\currtens\timepoints_for_tensor.mat'],'timepoints_for_tensor');

%% Select a neuron type
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\matlab glm training set\combine mat and python glms\consensus_idx_from_glm_when_normByGLMcoefIntegral.mat');
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\matlab glm training set\unitnames_glm.mat');
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\matlab glm training set\fromWhichSess_glm.mat');
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\python glm training set\py_all_glm_coef_butIndexedIntoMatCoefs.mat');
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\python glm training set\py_metrics_butIndexedIntoMatCoefs.mat');
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\matlab glm training set\combine mat and python glms\allZeroCoeffsAfterOutcome.mat');
whichCellTypeToTake=[1 2]; % 1 is success-continuing
if length(whichSess)==1
    % match these neurons to their type classification
    currnames=unitnames_glm(fromWhichSess_glm==whichSess);
    if length(currnames)==size(tensor,1)
        % neurons match
%         tensor=tensor(ismember(idx_from_glm(fromWhichSess_glm==whichSess),whichCellTypeToTake),:,:);
        cuethresh=1;
        tensorPart1=nanmean(tensor(idx_from_glm(fromWhichSess_glm==whichSess)==1 & cueco(fromWhichSess_glm==whichSess)>cuethresh,:,:),1);
        tensorPart2=nanmean(tensor(idx_from_glm(fromWhichSess_glm==whichSess)==2 & cueco(fromWhichSess_glm==whichSess)>cuethresh,:,:),1);
        tensorPart3=nanmean(tensor(idx_from_glm(fromWhichSess_glm==whichSess)==1 & cueco(fromWhichSess_glm==whichSess)<=cuethresh,:,:),1);
        tensorPart4=nanmean(tensor(idx_from_glm(fromWhichSess_glm==whichSess)==2 & cueco(fromWhichSess_glm==whichSess)<=cuethresh,:,:),1);
%         tensor=cat(1,tensorPart1,tensorPart2);
        tensor=[tensorPart1; tensorPart2; tensorPart3; tensorPart4];
%         morefortensor=cat(1,tensorPart1,tensorPart2);
%         tensor=cat(1,tensor,morefortensor);
        figure(); imagesc(nanmean(tensor(:,:,allLabels==0),3)); title('cued success');
        figure(); imagesc(nanmean(tensor(:,:,allLabels==1),3)); title('uncued success');
        figure(); imagesc(nanmean(tensor(:,:,allLabels==2),3)); title('cued failure');
        figure(); imagesc(nanmean(tensor(:,:,allLabels==3),3)); title('uncued failure');
        [meanOfAll,meanReal,meanShuffle]=trialTypeDecode(tensor,allLabels,timepoints_for_tensor); % Plot trial type decode
    elseif length(currnames)<size(tensor,1)
        disp(['length of currnames: ' num2str(length(currnames))]);
        disp(['length of neurons in tensor: ' num2str(size(tensor,1))]);
        error('neurons do not match');
        disp('Cutting tensor');
        %disp(currnames); 
        temp=table2cell(readtable('C:\Users\sabatini\Documents\currtens\tens_unit_names.csv','Format','%s','Delimiter', ',', 'HeaderLines', 0, 'ReadVariableNames', true));
        whichToTake=matchUpNames(temp,currnames);
        tensor=tensor(whichToTake==1,:,:);
    else
        disp(['length of currnames: ' num2str(length(currnames))]);
        disp(['length of neurons in tensor: ' num2str(size(tensor,1))]);
        disp('Not a match');
    end
end
disp(['Neurons of selected type: ' num2str(size(tensor,1))]);
cueco=nansum(py_all_glm_coef(:,[19:22]),2); % look for greater than 2
disp([max(cueco(fromWhichSess_glm==whichSess))]);
disp([cueco(fromWhichSess_glm==whichSess) py_metrics.cuecoef_over1sec(fromWhichSess_glm==whichSess) idx_from_glm(fromWhichSess_glm==whichSess)]);
save(['C:\Users\sabatini\Documents\currtens\tensor.mat'],'tensor');

%% Project single session tensor onto best TCA
TCAtooktimeafterzero=450*0.01; 
% interpolate current tensor to match TCA times, second dimension is time
interptens=[];
for i=1:size(tensor,1)
    for j=1:size(tensor,3)
        interptens(i,:,j)=interp1(timepoints_for_tensor,tensor(i,:,j),linspace(0,TCAtooktimeafterzero,450));
    end
end
trialWeightsOntoTypeTimeFactors=projectOntoCPdecomp('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM training set\TCA decomp using GLM training set\allconditions_cpmodel.mat',interptens,[1 2 3 4 5 6]);
tensor=trialWeightsOntoTypeTimeFactors;
save(['C:\Users\sabatini\Documents\currtens\tensor.mat'],'tensor');

%% GLM analysis
for i=446:450 %407:450 %351:381 %441:450
    if data_loc_array{i,13}==1
        continue
    end
    if strcmp(data_loc_array{i,8},'no_tbt')
        continue
    end
    % Make cue training set as combination of reach types
    makeComboTrainingSet({'cued_success','cued_failure','uncued_success','uncued_failure'},...
                         {'cuedSuccess','cuedFailure','uncuedSuccess','uncuedFailure'},'cue','cueAligned',data_loc_array{i,8});
    % Set up and fit GLM
    GLM_analysis(i,data_loc_array,10);
    close all;
end

% Read in unit names
plotUnitCriteria=[1 0 0 1 0]; getCriteriaForUnitsToPlot(plotUnitCriteria); dd_more=cell(1,length(dd)); 
for i=1:length(dd)
    dd_more{i}=[dd{i} sep 'cued_reach']; % just a placeholder to read in names
end
whichUnitsToGrab='_'; 
[unitbyunit_names.names,unitbyunit_names.excluded,unitbyunit_names.D1orD2taggingExpt,unitbyunit_names.D1taggedCells,unitbyunit_names.A2ataggedCells,unitbyunit_names.fromWhichSess,unitbyunit_names.fromWhichUnit]=alignToReadInUnitNames(dd_more,whichUnitsToGrab,settingsForStriatumUnitPlots);
% Read in GLM coefficients
for i=1:length(dd)
    dd_more{i}=[dd{i} sep 'matglm_trainingSet']; % just a placeholder to read in names
end
[all_glm_coef,unitnames_glm,fromWhichSess_glm,pva_glm]=readInGLMCoefs(dd_more,settingsForStriatumUnitPlots);
% Outliers for matlab glm
% outli=[312 313 319 320 322 324 374 376 377 378 379 380 381 382 383 386 387 388 389];
% all_glm_coef=all_glm_coef(~ismember(1:size(all_glm_coef,1),outli),:); unitnames_glm=unitnames_glm(~ismember(1:length(unitnames_glm),outli)); fromWhichSess_glm=fromWhichSess_glm(~ismember(1:length(fromWhichSess_glm),outli));
% f=find(cued_success_Response.excluded==0); fmore=find(unitbyunit_names.excluded==0); ftormv=find(~ismember(fmore,f));
% unitbyunit_names.names=unitbyunit_names.names(~ismember(1:length(unitbyunit_names.names),ftormv));
% unitbyunit_names.excluded=cued_success_Response.excluded;
getridisnan=find(isnan(indexGLMcellsIntoUnitNames));
all_glm_coef=all_glm_coef(~ismember(1:size(all_glm_coef,1),getridisnan),:); unitnames_glm=unitnames_glm(~ismember(1:length(unitnames_glm),getridisnan)); fromWhichSess_glm=fromWhichSess_glm(~ismember(1:length(fromWhichSess_glm),getridisnan));
indexGLMcellsIntoUnitNames=getNamesIndexIntoNamesList(unitnames_glm,unitbyunit_names);
% load('Z:\MICROSCOPE\Kim\WHISPER recs\Mar_6\20210625\SU aligned to behavior\forglm\output\neuron1_glm.mat')
% [ts,allco]=plotGLMcoef(all_glm_coef,0,feature_names,10*0.01,9,'mean'); title('python glm');
% load('Z:\MICROSCOPE\Kim\WHISPER recs\Mar_6\20210625\SU aligned to behavior\matglm\features_for_glm.mat');
% [ts,allco]=plotGLMcoef(all_glm_coef,[],fnames,10*0.01,9,'mean'); title('mat glm');
% [ts,allco]=plotGLMcoef(all_glm_coef(idx(indexGLMcellsIntoUnitNames)==1,:),[],fnames,10*0.01,9,'mean'); title('mat glm');
load('Z:\MICROSCOPE\Kim\WHISPER recs\Mar_2\20210805\SU aligned to behavior\matglm_trainingSet\neuron1_glm_coef.mat');
load('Z:\MICROSCOPE\Kim\WHISPER recs\Mar_2\20210805\SU aligned to behavior\matglm_trainingSet\features_for_glm.mat');
load('Z:\MICROSCOPE\Kim\WHISPER recs\Mar_2\20210805\SU aligned to behavior\matglm_trainingSet\shifts.mat');
[ts,allco]=plotGLMcoef(coef,[],fnames,10*0.01,nansum(shifts<0),'mean',false,[]); title('mat glm');
whichCoefToUse=[4 5 6]; studyGLMcoef(all_glm_coef,ts,whichCoefToUse);
metrics=getMetricsForAllGLMcoef(all_glm_coef,[],fnames,10*0.01,nansum(shifts<0));
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\matlab glm training set\unitnames_glm.mat'); mat_unitnames_glm=unitnames_glm;
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\python glm training set\unitnames_glm.mat');
indexPyGLMintoMatGLM=getNamesIndexIntoNamesList(unitnames_glm,mat_unitnames_glm);
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\matlab glm training set\metrics.mat'); mat_metrics=metrics; backup_mat_metrics=metrics;
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\matlab glm training set\all_glm_coef.mat');  mat_all_glm_coef=all_glm_coef;
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\matlab glm training set\fromWhichSess_glm.mat');
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\matlab glm training set\pvals_glm.mat');
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\python glm training set\metrics.mat'); f=fieldnames(mat_metrics);
for i=1:length(f)
    temp=mat_metrics.(f{i});
    temp(indexPyGLMintoMatGLM)=metrics.(f{i});
    mat_metrics.(f{i})=temp;
end
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\python glm training set\all_glm_coef.mat'); py_all_glm_coef=mat_all_glm_coef; py_all_glm_coef(indexPyGLMintoMatGLM,:)=all_glm_coef;
doingGLMfigures(py_all_glm_coef,mat_metrics,mat_all_glm_coef,backup_mat_metrics,fromWhichSess_glm,pva_glm);
% relate glm classification to TCA classification
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\matlab glm training set\combine mat and python glms\consensus_idx_from_glm_when_normByGLMcoefIntegral.mat');
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\test set\cued_success_Response.mat');
idx_g=nan(size(cued_success_Response.unitbyunit_x,1),1);
idx_g(indexGLMcellsIntoUnitNames(~isnan(indexGLMcellsIntoUnitNames)))=idx_from_glm(~isnan(indexGLMcellsIntoUnitNames));
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM training set\TCA decomp using GLM training set\idx_groupLabelsFromTCA.mat');
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM training set\TCA decomp using GLM training set\allcell_PCs.mat');
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\matlab glm training set\combine mat and python glms\figures\tsne.mat');
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\matlab glm training set\combine mat and python glms\figures\tsne_whichNanned.mat');
Ywithnans=nan(length(Ynanned),1); Ywithnans(Ynanned==0)=Y(:,1);
Y_tomatchTCA=nan(size(cued_success_Response.unitbyunit_x,1),1); 
Y_tomatchTCA(indexGLMcellsIntoUnitNames(~isnan(indexGLMcellsIntoUnitNames)))=Ywithnans(~isnan(indexGLMcellsIntoUnitNames),1);
figure(); scatter(Y_tomatchTCA,allcell_PCs.score(:,1)); forvio{1}=allcell_PCs.score(Y_tomatchTCA<0,1); forvio{2}=allcell_PCs.score(Y_tomatchTCA>=0,1);
figure(); violin(forvio);
[r,p]=corrcoef(Y_tomatchTCA(~isnan(Y_tomatchTCA)),allcell_PCs.score(~isnan(Y_tomatchTCA),1));
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
% a=load('Z:\MICROSCOPE\Kim\20230205 all SU alignments\all trials averaged not downsampled\cue.mat'); cue_Response=a.Response; 
a=load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM test set\excluded trials where opto during cue\cued_success_Response.mat'); cued_success_Response=a.cued_success_Response;  
a=load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM test set\excluded trials where opto during cue\cued_failureNotDrop_Response.mat'); cued_failure_Response=a.cued_failureNotDrop_Response; 
a=load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM test set\excluded trials where opto during cue\uncued_success_Response.mat'); uncued_success_Response=a.uncued_success_Response; 
a=load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM test set\excluded trials where opto during cue\uncued_failureNotDrop_Response.mat'); uncued_failure_Response=a.uncued_failureNotDrop_Response;
% a=load('Z:\MICROSCOPE\Kim\20230205 all SU alignments\all trials averaged not downsampled\cue_noReach.mat'); cue_noReach_Response=a.Response;
%a=load('Z:\MICROSCOPE\Kim\20221129 lab meeting\responses unit by unit\uncued_reach.mat'); uncued_reach_Response=a.Response;

% Remove units with firing rates too low
% toolow=zeros(size(nanmean(cued_success_Response.unitbyunit_y,2)));
% toolow=nanmean(cued_success_Response.unitbyunit_y,2)<0.1;
% trmv=cued_success_Response.excluded;
% trmv(toolow==1)=1;
% trmv=logical(trmv);
% cued_success_Response=removeUnitFromResponse(cued_success_Response,trmv);
% cued_failure_Response=removeUnitFromResponse(cued_failure_Response,trmv);
% uncued_failure_Response=removeUnitFromResponse(uncued_failure_Response,trmv);
% uncued_success_Response=removeUnitFromResponse(uncued_success_Response,trmv);

trial_n_cutoff=0;
out=plotVariousSUResponsesAlignedToBeh('matchUnitsAcrossResponses',excludeTooFewTrials(cued_success_Response,trial_n_cutoff,false),excludeTooFewTrials(cued_success_Response,trial_n_cutoff,false),excludeTooFewTrials(cued_failure_Response,trial_n_cutoff,false),excludeTooFewTrials(uncued_success_Response,trial_n_cutoff,false),excludeTooFewTrials(uncued_failure_Response,trial_n_cutoff,false));
% out=plotVariousSUResponsesAlignedToBeh('matchUnitsAcrossResponses',excludeTooFewTrials(cue_Response,trial_n_cutoff,false),excludeTooFewTrials(cued_success_Response,trial_n_cutoff,false),excludeTooFewTrials(cued_failure_Response,trial_n_cutoff,false),excludeTooFewTrials(uncued_success_Response,trial_n_cutoff,false),excludeTooFewTrials(uncued_failure_Response,trial_n_cutoff,false));
% cue_Response=out.Response1; 
cued_success_Response=out.Response2; cued_failure_Response=out.Response3; uncued_success_Response=out.Response4; uncued_failure_Response=out.Response5;
% out=plotVariousSUResponsesAlignedToBeh('matchUnitsAcrossResponses',excludeTooFewTrials(cue_noReach_Response,trial_n_cutoff,false),excludeTooFewTrials(cued_success_Response,trial_n_cutoff,false),excludeTooFewTrials(cued_failure_Response,trial_n_cutoff,false),excludeTooFewTrials(uncued_success_Response,trial_n_cutoff,false),excludeTooFewTrials(uncued_failure_Response,trial_n_cutoff,false));
% cue_noReach_Response=out.Response1;  cued_success_Response=out.Response2; cued_failure_Response=out.Response3; uncued_success_Response=out.Response4; uncued_failure_Response=out.Response5;
D1tag=cued_success_Response.D1tag(cued_success_Response.excluded==0); A2atag=cued_success_Response.A2atag(cued_success_Response.excluded==0); 
% save('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\training\for_data_matrix_D1vA2a.mat','D1tag','A2atag');

% takePointsBeforeZero=5; %15;
% takePointsAfterZero=70;

takePointsBeforeZero=30; %15;
takePointsAfterZero=450;

dataMatrix=setUpDataMatrix(cued_success_Response,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,takePointsBeforeZero,takePointsAfterZero);
% dataMatrix=setUpDataMatrix(cue_Response,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,takePointsBeforeZero,takePointsAfterZero);
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

clear newDataMatrix
for i=1:size(dataMatrix,3)
    temp=reshape(dataMatrix(:,:,i),size(dataMatrix(:,:,i),1),size(dataMatrix(:,:,i),2));
    %temp=downSampMatrix(reshape(dataMatrix(:,:,i),size(dataMatrix(:,:,i),1),size(dataMatrix(:,:,i),2)),10);
    for j=1:size(temp,1)
        temp(j,:)=smoothdata(temp(j,:),'gaussian',42);%10);
    end
%     newDataMatrix(:,:,i)=abs(diff(temp,1,2)); 
    newDataMatrix(:,:,i)=temp; 
end

% PCA, CCA, etc.
% CCA: Find orthogonal dimensions of max covariance between X=[] and Y=[]
plotN=6;
boot=1; % num iterations for bootstrap
[groupLabelsFromTCA,cuez]=principaledCA(newDataMatrix,{'units','time','conditions'},plotN,boot,'high');

% zscore_cueR=cue_Response.unitbyunit_y; zscore_cueR(zscore_cueR<0.0001)=0;
% zscore_cueR=(zscore_cueR)./repmat(std(zscore_cueR(:,500:1000),[],2,'omitnan'),1,size(zscore_cueR,2));
% figure(); plot(nanmean(zscore_cueR,1));
% cuedcellresponse=nanmean(zscore_cueR(:,286:328),2)-nanmean(zscore_cueR(:,210:282),2);
% cuedcellresponse(cuedcellresponse>35)=35;
% figure(); histogram(cuedcellresponse,500);

% figure(); plot(nanmean(cued_failure_Response.unitbyunit_x,1),nanmean(cued_failure_Response.unitbyunit_y(idx==2 & cuedcellresponse>0.5,:),1),'Color','b');
% hold on; plot(nanmean(cued_failure_Response.aligncomp_x,1),nanmean(cued_failure_Response.aligncomp_y(idx==2 & cuedcellresponse>0.5,:),1),'Color','c'); title('failure cued');
% figure(); plot(nanmean(cued_success_Response.unitbyunit_x,1),nanmean(cued_success_Response.unitbyunit_y(idx==2 & cuedcellresponse>0.5,:),1),'Color','b');
% hold on; plot(nanmean(cued_success_Response.aligncomp_x,1),nanmean(cued_success_Response.aligncomp_y(idx==2 & cuedcellresponse>0.5,:),1),'Color','c'); title('success cued');
% figure(); plot(nanmean(cued_failure_Response.unitbyunit_x,1),nanmean(cued_failure_Response.unitbyunit_y(idx==2 & cuedcellresponse<0.5,:),1),'Color','b');
% hold on; plot(nanmean(cued_failure_Response.aligncomp_x,1),nanmean(cued_failure_Response.aligncomp_y(idx==2 & cuedcellresponse<0.5,:),1),'Color','c'); title('failure medium');
% figure(); plot(nanmean(cued_success_Response.unitbyunit_x,1),nanmean(cued_success_Response.unitbyunit_y(idx==2 & cuedcellresponse<0.5,:),1),'Color','b');
% hold on; plot(nanmean(cued_success_Response.aligncomp_x,1),nanmean(cued_success_Response.aligncomp_y(idx==2 & cuedcellresponse<0.5,:),1),'Color','c'); title('success medium');
% figure(); plot(nanmean(cued_failure_Response.unitbyunit_x,1),nanmean(cued_failure_Response.unitbyunit_y(idx==2 & cuedcellresponse<0,:),1),'Color','b');
% hold on; plot(nanmean(cued_failure_Response.aligncomp_x,1),nanmean(cued_failure_Response.aligncomp_y(idx==2 & cuedcellresponse<0,:),1),'Color','c'); title('failure uncued');
% figure(); plot(nanmean(cued_success_Response.unitbyunit_x,1),nanmean(cued_success_Response.unitbyunit_y(idx==2 & cuedcellresponse<0,:),1),'Color','b');
% hold on; plot(nanmean(cued_success_Response.aligncomp_x,1),nanmean(cued_success_Response.aligncomp_y(idx==2 & cuedcellresponse<0,:),1),'Color','c'); title('success uncued');

% analyzeProbabilityOfOnAfterOutcome(dd,[0 2],[],[]);
% plotSU_contextAndOutcome('Z:\MICROSCOPE\Kim\WHISPER recs\Mar_1\20210803\SU aligned to behavior',failure_off);

%% Plot unit summaries according to groupLabelsFromTCA
% if ~exist('cuez','var')
%     temp=normalizeDataMatrix(newDataMatrix,[2 3],'sd'); cuez=reshape(max(temp(:,1:150,1),[],2,'omitnan'),size(newDataMatrix,1),1);
%     cuez=log(cuez); cuez(cuez<-2)=-2; figure(); histogram(cuez,200); 
% end
% a=load('Z:\MICROSCOPE\Kim\20230205 all SU alignments\all trials averaged not downsampled\uncued_reach.mat'); uncuedReach_Response=a.Response;  
% a=load('Z:\MICROSCOPE\Kim\20230205 all SU alignments\all trials averaged not downsampled\cue_noReach.mat'); cue_noReach_Response=a.Response;
% a=load('Z:\MICROSCOPE\Kim\20230205 all SU alignments\all trials averaged not downsampled\cued_reach.mat'); cuedReach_Response=a.Response;
% % fixing cuedReach response because missing units compared to others
% load('Z:\MICROSCOPE\Kim\20230205 all SU alignments\all trials averaged not downsampled\unitbyunit_names.mat');
% load('Z:\MICROSCOPE\Kim\20230205 all SU alignments\all trials averaged not downsampled\unitbyunitnames_cuedreach.mat')
% cuedReach_Response=matchExcludedBasedOnUnitNames(unitbyunitnames_cuedreach,unitbyunit_names,cuedReach_Response);

% out=plotVariousSUResponsesAlignedToBeh('matchUnitsAcrossResponses',cue_noReach_Response,cued_success_Response,[],[],[]);
% cue_noReach_Response=out.Response1;  
% cuedReach_Response=cued_success_Response; cuedReach_Response.unitbyunit_y=(cued_success_Response.unitbyunit_y+cued_failure_Response.unitbyunit_y(:,1:size(cued_success_Response.unitbyunit_y,2)))./2; 
% uncuedReach_Response=uncued_success_Response; uncuedReach_Response.unitbyunit_y=(uncued_success_Response.unitbyunit_y+uncued_failure_Response.unitbyunit_y(:,1:size(uncued_success_Response.unitbyunit_y,2)))./2; 
% cuez=getCueTunedUnits(cue_noReach_Response,uncuedReach_Response,'vs_uncued_reach','max'); % method 3rd arg can be 'vs_uncued_reach' or 'cue_vs_baseline' or 'justcue'

clear r
% load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\training\TCA\idx_groupLabelsFromTCA.mat'); temp=idx; idx(temp==2)=1; idx(temp==1)=2; % labels flipped for this TCA
% load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\tensor regression\rank 2\idx.mat'); 
% load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\tensor regression\BEST SO FAR\currtens\idx.mat'); 
% load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\tensor regression\BEST SO FAR\currtens\isHighWeight_train.mat');
% load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\tensor regression\BEST SO FAR\currtens\isVeryHighWeight_train.mat');

% load('C:\Users\sabatini\Documents\currtens BEST SO FAR\factor_0.mat');
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\tensor regression\BEST SO FAR\currtens THIS IS GOOD\factor_0.mat');
idx=factor(:,1)>=factor(:,2);
isHighWeight_train=factor(:,1)>0.0262 | factor(:,2)>0.0207; % top 100 neurons for each factor
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\matlab glm training set\combine mat and python glms\consensus_idx_from_glm_when_normByGLMcoefIntegral.mat');

usingGLMidx=false;
if usingGLMidx==true

    % CONSIDER THROWING OUT CUED TRIALS WHERE REACH SEEMS
    % PREEMPTIVE!!!!!!!!!!!!!!!!!!!!!!!!

    load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\unitbyunit_names_to_match_cued_success_Response.mat');
    load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\matlab glm training set\combine mat and python glms\consensus_idx_from_glm_when_normByGLMcoefIntegral.mat');
    load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM test set\excluded trials where opto during cue\cued_success_Response.mat');
    load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\python glm training set\py_all_glm_coef_butIndexedIntoMatCoefs.mat');
    load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\python glm training set\py_metrics_butIndexedIntoMatCoefs.mat');
    load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\matlab glm training set\unitnames_glm.mat');

%     temp=py_metrics.allSucc_sustained+py_metrics.cXsucc_sustained>0.01;
% %     temp=ones(size(py_metrics.allDrop_sustained))==1; %-py_metrics.cXsucc_sustained>0;
%     [n,x]=histcounts(py_metrics.allSucc_sustained+py_metrics.cXsucc_sustained>0.01,-1-0.005:0.01:1+0.005);
%     [n,x]=cityscape_hist(n,x); figure(); plot(x,n,'Color','k'); ylim([0 50]);
%     idx_from_glm(~temp)=2; % failure-continuing
%     idx_from_glm(temp)=1; % success-continuing
%     f=find(py_metrics.allSucc_sustained+py_metrics.cXsucc_sustained>0.01);
%     temp=py_metrics.cXsucc_sustained-py_metrics.cXfail_sustained;
%     idx_from_glm=temp;
    
    conditionNow=idx_from_glm==2; 
    names_of_units=unitnames_glm(conditionNow & ~isnan(indexGLMcellsIntoUnitNames));
    indexGLMcellsIntoUnitNames=getNamesIndexIntoNamesList(unitnames_glm,unitbyunit_names); 
    sessions_of_units=cued_success_Response.fromWhichSess(indexGLMcellsIntoUnitNames(conditionNow & ~isnan(indexGLMcellsIntoUnitNames)));
    cued_success_Response.carrythroughindsintowhich_units=nan(size(cued_success_Response.fromWhichSess));
    cued_success_Response.carrythroughindsintowhich_units(indexGLMcellsIntoUnitNames(conditionNow & ~isnan(indexGLMcellsIntoUnitNames)))=1:length(names_of_units);
    py_metrics.activeMoreBeforeCuedReach=nansum(py_all_glm_coef(:,[1:16 427:427+20 498:498+20 569:569+20]),2);
    py_metrics.activeMoreBeforeAnyReach=nansum(py_all_glm_coef(:,[214:214+20 285:285+20 356:356+20]),2);
    py_metrics.cuecoef_over1sec=nansum(py_all_glm_coef(:,[20:30]),2);
    py_metrics.cuecoef_overall=nansum(py_all_glm_coef(:,[19:71]),2);
    py_metrics.cuecoef_at1sec=nansum(py_all_glm_coef(:,[30:33]),2);
    py_metrics.cXsucc_over1sec=nansum(py_all_glm_coef(:,[427:427+19]),2);
    py_metrics.allSucc_over1sec=nansum(py_all_glm_coef(:,[214:214+19]),2);
    py_metrics.cXsucc_overpoint5sec=nansum(py_all_glm_coef(:,[427:427+10]),2);
    py_metrics.allSucc_overpoint5sec=nansum(py_all_glm_coef(:,[214:214+10]),2);    
    idx=nan(size(cued_success_Response.unitbyunit_x,1),1);
    idx(indexGLMcellsIntoUnitNames(~isnan(indexGLMcellsIntoUnitNames)))=idx_from_glm(~isnan(indexGLMcellsIntoUnitNames)); cued_success_Response.idx=idx;
%     whichGLMinds=[286:286+71 356:356+71 499:499+71 568:568+71];
    whichGLMinds=[1:71];
    cued_success_Response=addMetricsToResponse(cued_success_Response,py_metrics,py_all_glm_coef,indexGLMcellsIntoUnitNames,whichGLMinds);

%     temp=cued_success_Response.glmcoef_index21+cued_success_Response.glmcoef_index22+cued_success_Response.glmcoef_index23+cued_success_Response.glmcoef_index24+cued_success_Response.glmcoef_index25+cued_success_Response.glmcoef_index26;
%     grp1=temp>0.1 & idx==1;
%     grp2=temp>0.1 & idx==2;
% %     grp1=~(temp>0) & cued_success_Response.idx==1;
% %     grp2=~(temp>0) & cued_success_Response.idx==2;
%     idx(:)=nan;
%     idx(grp2)=2; 
%     idx(grp1)=1;
%     cued_success_Response.idx=idx;

    r{1}=cued_success_Response;
    load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM test set\excluded trials where opto during cue\cued_failure_Response.mat'); cued_failure_Response.idx=idx; r{2}=cued_failure_Response;
    load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM test set\excluded trials where opto during cue\uncued_failure_Response.mat'); uncued_failure_Response.idx=idx; r{3}=uncued_failure_Response;
    load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM test set\excluded trials where opto during cue\uncued_success_Response.mat'); uncued_success_Response.idx=idx; r{4}=uncued_success_Response;
%     load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM test set\excluded trials where opto during cue\cued_drop_Response.mat'); r{5}=cued_drop_Response;
%     load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM test set\cued_failureNotDrop_Response.mat'); r{6}=cued_failureNotDrop_Response;
%     load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM test set\excluded trials where opto during cue\uncued_drop_Response.mat'); r{6}=uncued_drop_Response;
%     load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM test set\uncued_failureNotDrop_Response.mat'); r{8}=uncued_failureNotDrop_Response;
%     load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\all trials\cued_reach_Response.mat'); r{9}=cued_reach_Response;
%     load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\all trials\uncued_reach_Response.mat'); r{10}=uncued_reach_Response;
%     load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM test set\cued_failure_noReach_Response.mat'); r{11}=cued_failure_noReach_Response;
%     load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM test set\uncued_failure_noReach_Response.mat'); r{12}=uncued_failure_noReach_Response;
%     load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM test set\all_success_Response.mat'); r{13}=all_success_Response;
%     load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM test set\all_failure_Response.mat'); r{14}=all_failure_Response;
else
    load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\unitbyunit_names_to_match_cued_success_Response.mat');
    load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM test set\excluded trials where opto during cue\cued_success_Response.mat');
    load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\tensor regression\BEST SO FAR\currtens\excludedInCurrent_cued_success_Response.mat');
    nownotexcl=find(cued_success_Response.excluded==0); wasnotexcl=find(excludedInCurrent_cued_success_Response==0);

    cued_success_Response.isHighWeight=nan(size(cued_success_Response.ns));
    cued_success_Response.isVeryHighWeight=nan(size(cued_success_Response.ns));
    cued_success_Response.idx=nan(size(cued_success_Response.ns)); 
    cued_success_Response.idx(ismember(nownotexcl,wasnotexcl))=idx; 
    cued_success_Response.isHighWeight(ismember(nownotexcl,wasnotexcl))=isHighWeight_train';
    cued_success_Response.isVeryHighWeight(ismember(nownotexcl,wasnotexcl))=isVeryHighWeight_train;
    
    % cued_success_Response.idx=idx;
    load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\python glm training set\py_all_glm_coef_butIndexedIntoMatCoefs.mat');
    load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\python glm training set\py_metrics_butIndexedIntoMatCoefs.mat');
    load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\python glm training set\allr2scores.mat');
    load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\matlab glm training set\unitnames_glm.mat');
    indexGLMcellsIntoUnitNames=getNamesIndexIntoNamesList(unitnames_glm,unitbyunit_names); 
    load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\matlab glm training set\combine mat and python glms\maxNormed_consensus_glm_coef.mat');
    maxnormedglmcoef=nan(size(cued_success_Response.unitbyunit_x,1),size(nama_consensus_all_glm_coef,2));
    maxnormedglmcoef(indexGLMcellsIntoUnitNames(~isnan(indexGLMcellsIntoUnitNames)),:)=nama_consensus_all_glm_coef(~isnan(indexGLMcellsIntoUnitNames),:);
    py_metrics.allr2scores=allr2scores;
    py_metrics.idx_from_glm=idx_from_glm;
    py_metrics.activeMoreBeforeCuedReach=nansum(py_all_glm_coef(:,[1:16 427:427+20 498:498+20 569:569+20]),2);
    py_metrics.activeMoreBeforeAnyReach=nansum(py_all_glm_coef(:,[214:214+20 285:285+20 356:356+20]),2);
    py_metrics.cuecoef_over1sec=nansum(py_all_glm_coef(:,[20:30]),2);
    py_metrics.cuecoef_overall=nansum(py_all_glm_coef(:,[19:71]),2);
    py_metrics.cuecoef_at1sec=nansum(py_all_glm_coef(:,[30:33]),2);
    py_metrics.cXsucc_over1sec=nansum(py_all_glm_coef(:,[427:427+19]),2);
    py_metrics.allSucc_over1sec=nansum(py_all_glm_coef(:,[214:214+19]),2);
    py_metrics.cXsucc_overpoint5sec=nansum(py_all_glm_coef(:,[427:427+10]),2);
    py_metrics.allSucc_overpoint5sec=nansum(py_all_glm_coef(:,[214:214+10]),2);
    py_metrics.precuebaseline=nansum(py_all_glm_coef(:,[1:1+10]),2);

    % [num2str((1:640)') feature_names]
    % Tensor regression temporal factors started at time 0, wrt reach,
    % continues for 4.5 secs
    % GLM coefficients according to shifts are from -2 to 5 sec wrt event
    % GLM inds for 0 to 2 sec: shift names 0 to 20, 20 to 40 index
    % GLM inds for 2 to 4.5 sec: shift names 21 to 45, 41 to 65 index
    % GLM inds for 2 to 5 sec: shift names 21 to 50, 41 to 70 index
    % all success starts at 214
    % sus1to5sec
    py_metrics.allfail_sus1to5sec=(nanmean(py_all_glm_coef(:,[305+10:305+50]),2)+nanmean(py_all_glm_coef(:,[376+10:376+50]),2))/2;
    py_metrics.alldrop_sus1to5sec=nanmean(py_all_glm_coef(:,[305+10:305+50]),2);
    py_metrics.allmiss_sus1to5sec=nanmean(py_all_glm_coef(:,[376+10:376+50]),2);
    py_metrics.cXfail_sus1to5sec=(nanmean(py_all_glm_coef(:,[518+10:518+50]),2)+nanmean(py_all_glm_coef(:,[589+10:589+50]),2))/2;
    py_metrics.cXdrop_sus1to5sec=nanmean(py_all_glm_coef(:,[518+10:518+50]),2);
    py_metrics.cXmiss_sus1to5sec=nanmean(py_all_glm_coef(:,[589+10:589+50]),2);
    py_metrics.allsucc_sus1to5sec=nanmean(py_all_glm_coef(:,[234+10:234+50]),2);
    py_metrics.cXsucc_sus1to5sec=nanmean(py_all_glm_coef(:,[447+10:447+50]),2);
    py_metrics.allsucc_sus0to5sec=nanmean(py_all_glm_coef(:,[234:234+50]),2);
    py_metrics.cXsucc_sus0to5sec=nanmean(py_all_glm_coef(:,[447:447+50]),2);
    py_metrics.allfail_minus2to0sec=(nanmean(py_all_glm_coef(:,[305-20:305]),2)+nanmean(py_all_glm_coef(:,[376-20:376]),2))/2;
    py_metrics.cXfail_minus2to0sec=(nanmean(py_all_glm_coef(:,[518-20:518]),2)+nanmean(py_all_glm_coef(:,[589-20:589]),2))/2;
    py_metrics.allsucc_minus2to0sec=nanmean(py_all_glm_coef(:,[234-20:234]),2);
    py_metrics.cXsucc_minus2to0sec=nanmean(py_all_glm_coef(:,[447-20:447]),2);
    py_metrics.reachCoefAv=nanmean(py_all_glm_coef(:,[640:710]),2);

    firsthalfie=nanmean(py_all_glm_coef(:,[230:234]),2);
    secondhalfie=nanmean(py_all_glm_coef(:,[234:284]),2);
    py_metrics.beforeaftersuccess_modulation_index=(secondhalfie-firsthalfie)./(abs(secondhalfie)+abs(firsthalfie));
    
    firsthalfie=nanmean(py_all_glm_coef(:,[301:305]),2);
    secondhalfie=nanmean(py_all_glm_coef(:,[305:355]),2);
    py_metrics.beforeafterdrop_modulation_index=(secondhalfie-firsthalfie)./(abs(secondhalfie)+abs(firsthalfie));

    firsthalfie=nanmean(py_all_glm_coef(:,[372:376]),2);
    secondhalfie=nanmean(py_all_glm_coef(:,[376:426]),2);
    py_metrics.beforeaftermiss_modulation_index=(secondhalfie-firsthalfie)./(abs(secondhalfie)+abs(firsthalfie));
   
    py_metrics.beforeafterfailure_modulation_index=(py_metrics.beforeafterdrop_modulation_index+py_metrics.beforeaftermiss_modulation_index)/2;

    firsthalfie=nanmean(py_all_glm_coef(:,[234:254]),2);
    secondhalfie=nanmean(py_all_glm_coef(:,[255:284]),2);
    py_metrics.allsuccess_modulation_index=(secondhalfie-firsthalfie)./(abs(secondhalfie)+abs(firsthalfie));
    py_metrics.allsuccess_change=(secondhalfie-firsthalfie);
    firsthalfie=nanmean(py_all_glm_coef(:,[447:467]),2);
    secondhalfie=nanmean(py_all_glm_coef(:,[468:497]),2);
    py_metrics.cXsuccess_modulation_index=(secondhalfie-firsthalfie)./(abs(secondhalfie)+abs(firsthalfie));
    firsthalfie=nanmean(py_all_glm_coef(:,[305:325]),2);
    secondhalfie=nanmean(py_all_glm_coef(:,[326:355]),2);
    py_metrics.alldrop_modulation_index=(secondhalfie-firsthalfie)./(abs(secondhalfie)+abs(firsthalfie));
    firsthalfie=nanmean(py_all_glm_coef(:,[376:396]),2);
    secondhalfie=nanmean(py_all_glm_coef(:,[397:426]),2);
    py_metrics.allmiss_modulation_index=(secondhalfie-firsthalfie)./(abs(secondhalfie)+abs(firsthalfie));
    py_metrics.allmiss_change=(secondhalfie-firsthalfie);
    py_metrics.allfailure_modulation_index=(py_metrics.alldrop_modulation_index+py_metrics.allmiss_modulation_index)/2;
    firsthalfie=nanmean(py_all_glm_coef(:,[518:538]),2);
    secondhalfie=nanmean(py_all_glm_coef(:,[539:568]),2);
    py_metrics.cXdrop_modulation_index=(secondhalfie-firsthalfie)./(abs(secondhalfie)+abs(firsthalfie));
    firsthalfie=nanmean(py_all_glm_coef(:,[589:609]),2);
    secondhalfie=nanmean(py_all_glm_coef(:,[610:639]),2);
    py_metrics.cXmiss_modulation_index=(secondhalfie-firsthalfie)./(abs(secondhalfie)+abs(firsthalfie));
    py_metrics.cXfailure_modulation_index=(py_metrics.cXdrop_modulation_index+py_metrics.cXmiss_modulation_index)/2;
    py_metrics.combosuccess_modulation_index=(py_metrics.allsuccess_modulation_index+py_metrics.cXsuccess_modulation_index)/2;
    py_metrics.combofailure_modulation_index=(py_metrics.allfailure_modulation_index+py_metrics.cXfailure_modulation_index)/2;

    whichGLMinds=[1:71];
%     whichGLMinds=[];
    cued_success_Response=addMetricsToResponse(cued_success_Response,py_metrics,py_all_glm_coef,indexGLMcellsIntoUnitNames,whichGLMinds);
    cued_success_Response.consensus_idx=zeros(size(cued_success_Response.idx));
    cued_success_Response.consensus_idx(cued_success_Response.idx_from_glm==1 & cued_success_Response.idx==1)=1; % 1 is succ-continuing, 0 is fail-continuing

    clearcell=abs(cued_success_Response.allSucc_sustained-cued_success_Response.allFail_sustained)>0.0005;
    figure(); scatter(cued_success_Response.combosuccess_modulation_index(clearcell & cued_success_Response.allr2scores>0 & cued_success_Response.idx_from_glm==1),cued_success_Response.combofailure_modulation_index(clearcell & cued_success_Response.allr2scores>0 & cued_success_Response.idx_from_glm==1),[],'b');
    hold on; scatter(cued_success_Response.combosuccess_modulation_index(clearcell & cued_success_Response.allr2scores>0 & cued_success_Response.idx_from_glm==2),cued_success_Response.combofailure_modulation_index(clearcell & cued_success_Response.allr2scores>0 & cued_success_Response.idx_from_glm==2),[],'r');
    scatter(cued_success_Response.cXsuccess_modulation_index(clearcell & cued_success_Response.allr2scores>0 & cued_success_Response.idx_from_glm==1),cued_success_Response.cXfailure_modulation_index(clearcell & cued_success_Response.allr2scores>0 & cued_success_Response.idx_from_glm==1),[],'b');
    hold on; scatter(cued_success_Response.cXsuccess_modulation_index(clearcell & cued_success_Response.allr2scores>0 & cued_success_Response.idx_from_glm==2),cued_success_Response.cXfailure_modulation_index(clearcell & cued_success_Response.allr2scores>0 & cued_success_Response.idx_from_glm==2),[],'r');
    scatter(cued_success_Response.allsuccess_modulation_index(clearcell & cued_success_Response.allr2scores>0 & cued_success_Response.idx_from_glm==1),cued_success_Response.allfailure_modulation_index(clearcell & cued_success_Response.allr2scores>0 & cued_success_Response.idx_from_glm==1),[],'b');
    hold on; scatter(cued_success_Response.allsuccess_modulation_index(clearcell & cued_success_Response.allr2scores>0 & cued_success_Response.idx_from_glm==2),cued_success_Response.allfailure_modulation_index(clearcell & cued_success_Response.allr2scores>0 & cued_success_Response.idx_from_glm==2),[],'r');
    xie=-1:0.01:1; yie=-0.75*exp(-1.1*(xie-0.6))+1.5; hold all; plot(xie,yie);
    title('GLM classification');

    figure(); scatter(cued_success_Response.combosuccess_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.consensus_idx==1),cued_success_Response.combofailure_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.consensus_idx==1),[],'b');
    hold on; scatter(cued_success_Response.combosuccess_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.consensus_idx==0),cued_success_Response.combofailure_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.consensus_idx==0),[],'r');
    figure(); scatter(cued_success_Response.allsuccess_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.consensus_idx==1),cued_success_Response.allfailure_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.consensus_idx==1),[],'b');
    hold on; scatter(cued_success_Response.allsuccess_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.consensus_idx==0),cued_success_Response.allfailure_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.consensus_idx==0),[],'r');
    scatter(cued_success_Response.cXsuccess_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.consensus_idx==1),cued_success_Response.cXfailure_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.consensus_idx==1),[],'b');
    hold on; scatter(cued_success_Response.cXsuccess_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.consensus_idx==0),cued_success_Response.cXfailure_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.consensus_idx==0),[],'r');
    figure(); scatter(cued_success_Response.combosuccess_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.consensus_idx==1),cued_success_Response.combofailure_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.consensus_idx==1),[],'k');
    hold on; scatter(cued_success_Response.combosuccess_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.consensus_idx==0),cued_success_Response.combofailure_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.consensus_idx==0),[],'r');
    
%     scatter(cued_success_Response.combosuccess_modulation_index(cued_success_Response.cXsuccess_modulation_index>cued_success_Response.allsuccess_modulation_index & cued_success_Response.isHighWeight==1 & cued_success_Response.consensus_idx==1),cued_success_Response.combofailure_modulation_index(cued_success_Response.cXsuccess_modulation_index>cued_success_Response.allsuccess_modulation_index & cued_success_Response.isHighWeight==1 & cued_success_Response.consensus_idx==1),[],'b','filled');
%     hold on; scatter(cued_success_Response.combosuccess_modulation_index(cued_success_Response.cXfailure_modulation_index>cued_success_Response.allfailure_modulation_index & cued_success_Response.isHighWeight==1 & cued_success_Response.consensus_idx==0),cued_success_Response.combofailure_modulation_index(cued_success_Response.cXfailure_modulation_index>cued_success_Response.allfailure_modulation_index & cued_success_Response.isHighWeight==1 & cued_success_Response.consensus_idx==0),[],'r','filled');
    
    scatter(cued_success_Response.allsuccess_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.consensus_idx==1),cued_success_Response.allfailure_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.consensus_idx==1),[],'k',"square");
    hold on; scatter(cued_success_Response.allsuccess_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.consensus_idx==0),cued_success_Response.allfailure_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.consensus_idx==0),[],'r',"square");
    scatter(cued_success_Response.cXsuccess_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.consensus_idx==1),cued_success_Response.cXfailure_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.consensus_idx==1),[],'k',"diamond");
    hold on; scatter(cued_success_Response.cXsuccess_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.consensus_idx==0),cued_success_Response.cXfailure_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.consensus_idx==0),[],'r',"diamond");
    xie=-1:0.01:1; yie=-0.75*exp(-1.1*(xie-0.6))+1.5; hold all; plot(xie,yie);

    % What I used for Fig 5 supp
    cmap=[0, 0.75, 0.75; 0.4940, 0.1840, 0.5560];
    figure(); s=scatter(cued_success_Response.combosuccess_modulation_index(cued_success_Response.idx_from_glm==2),cued_success_Response.combofailure_modulation_index(cued_success_Response.idx_from_glm==2),100,cmap(2,:),'filled');
    s.MarkerFaceAlpha=0.4; s.MarkerEdgeAlpha=0.4;
    hold on;
    s=scatter(cued_success_Response.combosuccess_modulation_index(cued_success_Response.idx_from_glm==1),cued_success_Response.combofailure_modulation_index(cued_success_Response.idx_from_glm==1),100,cmap(1,:),'filled');
    s.MarkerFaceAlpha=0.4; s.MarkerEdgeAlpha=0.4;
    s=scatter(cued_success_Response.combosuccess_modulation_index(cued_success_Response.idx==2 & cued_success_Response.isHighWeight==1),cued_success_Response.combofailure_modulation_index(cued_success_Response.idx==2 & cued_success_Response.isHighWeight==1),100,cmap(2,:),'filled');
    s.MarkerFaceAlpha=0.4; s.MarkerEdgeAlpha=0.4;
    hold on;
    s=scatter(cued_success_Response.combosuccess_modulation_index(cued_success_Response.idx==1 & cued_success_Response.isHighWeight==1),cued_success_Response.combofailure_modulation_index(cued_success_Response.idx==1 & cued_success_Response.isHighWeight==1),100,cmap(1,:),'filled');
    s.MarkerFaceAlpha=0.4; s.MarkerEdgeAlpha=0.4;
    %figure(); scatter(cued_success_Response.combosuccess_modulation_index(cued_success_Response.allr2scores>0 & cued_success_Response.idx_from_glm==1),cued_success_Response.combofailure_modulation_index(cued_success_Response.allr2scores>0 & cued_success_Response.idx_from_glm==1),[],'b');
    %hold on; scatter(cued_success_Response.combosuccess_modulation_index(cued_success_Response.allr2scores>0 & cued_success_Response.idx_from_glm==2),cued_success_Response.combofailure_modulation_index(cued_success_Response.allr2scores>0 & cued_success_Response.idx_from_glm==2),[],'r');
    figure(); scatter(cued_success_Response.cXsuccess_modulation_index(cued_success_Response.allr2scores>0 & cued_success_Response.idx_from_glm==1),cued_success_Response.cXfailure_modulation_index(cued_success_Response.allr2scores>0 & cued_success_Response.idx_from_glm==1),[],'b');
    hold on; scatter(cued_success_Response.cXsuccess_modulation_index(cued_success_Response.allr2scores>0 & cued_success_Response.idx_from_glm==2),cued_success_Response.cXfailure_modulation_index(cued_success_Response.allr2scores>0 & cued_success_Response.idx_from_glm==2),[],'r');
    scatter(cued_success_Response.allsuccess_modulation_index(cued_success_Response.allr2scores>0 & cued_success_Response.idx_from_glm==1),cued_success_Response.allfailure_modulation_index(cued_success_Response.allr2scores>0 & cued_success_Response.idx_from_glm==1),[],'b');
    hold on; scatter(cued_success_Response.allsuccess_modulation_index(cued_success_Response.allr2scores>0 & cued_success_Response.idx_from_glm==2),cued_success_Response.allfailure_modulation_index(cued_success_Response.allr2scores>0 & cued_success_Response.idx_from_glm==2),[],'r');
    xie=-1:0.01:1; yie=-0.75*exp(-1.1*(xie-0.6))+1.5; hold all; plot(xie,yie);
    title('GLM classification');
    %figure(); scatter(cued_success_Response.combosuccess_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.idx==0),cued_success_Response.combofailure_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.idx==0),[],'b');
    %hold on; scatter(cued_success_Response.combosuccess_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.idx==1),cued_success_Response.combofailure_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.idx==1),[],'r');
    figure(); scatter(cued_success_Response.cXsuccess_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.idx==0),cued_success_Response.cXfailure_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.idx==0),[],'b');
    hold on; scatter(cued_success_Response.cXsuccess_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.idx==1),cued_success_Response.cXfailure_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.idx==1),[],'r');
    scatter(cued_success_Response.allsuccess_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.idx==0),cued_success_Response.allfailure_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.idx==0),[],'b');
    hold on; scatter(cued_success_Response.allsuccess_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.idx==1),cued_success_Response.allfailure_modulation_index(cued_success_Response.isHighWeight==1 & cued_success_Response.idx==1),[],'r');
    xie=-1:0.01:1; yie=-0.75*exp(-1.1*(xie-0.6))+1.5; hold all; plot(xie,yie);
    title('Tensor reg classification');
%     cued_success_Response.consensus_idx=nan(size(cued_success_Response.idx));
%     cued_success_Response.consensus_idx(cued_success_Response.idx_from_glm==1 | cued_success_Response.idx==0 & ~(cued_success_Response.idx_from_glm==2 | cued_success_Response.idx==1))=1; % 1 is succ-continuing, 0 is fail-continuing
%     cued_success_Response.consensus_idx(cued_success_Response.idx_from_glm==2 | cued_success_Response.idx==1 & ~(cued_success_Response.idx_from_glm==1 | cued_success_Response.idx==0))=2;
    % CHANGE LABELS FOR IDX TO MATCH OTHERS
    temp=cued_success_Response.idx; cued_success_Response.idx(temp==0)=1; cued_success_Response.idx(temp==1)=2;
    %figure(); scatter(cued_success_Response.combosuccess_modulation_index(cued_success_Response.consensus_idx==1),cued_success_Response.combofailure_modulation_index(cued_success_Response.consensus_idx==1),[],'b');
    %hold on; scatter(cued_success_Response.combosuccess_modulation_index(cued_success_Response.consensus_idx==2),cued_success_Response.combofailure_modulation_index(cued_success_Response.consensus_idx==2),[],'r');
    figure(); scatter(cued_success_Response.allsuccess_modulation_index(cued_success_Response.consensus_idx==1),cued_success_Response.allfailure_modulation_index(cued_success_Response.consensus_idx==1),[],'b');
    hold on; scatter(cued_success_Response.allsuccess_modulation_index(cued_success_Response.consensus_idx==2),cued_success_Response.allfailure_modulation_index(cued_success_Response.consensus_idx==2),[],'r');
    scatter(cued_success_Response.cXsuccess_modulation_index(cued_success_Response.consensus_idx==1),cued_success_Response.cXfailure_modulation_index(cued_success_Response.consensus_idx==1),[],'b');
    hold on; scatter(cued_success_Response.cXsuccess_modulation_index(cued_success_Response.consensus_idx==2),cued_success_Response.cXfailure_modulation_index(cued_success_Response.consensus_idx==2),[],'r');
    cued_success_Response.idx_from_boundary=nan(size(cued_success_Response.idx));
    cued_success_Response.idx_from_boundary(cued_success_Response.combosuccess_modulation_index>0.45)=1;
    cued_success_Response.idx_from_boundary(cued_success_Response.combosuccess_modulation_index<-0.52)=2;
%     cued_success_Response.idx_far_from_boundary=nan(size(cued_success_Response.idx));
%     cued_success_Response.idx_far_from_boundary(cued_success_Response.combosuccess_modulation_index>0.2 & cued_success_Response.combosuccess_modulation_index~=1)=1;
%     cued_success_Response.idx_far_from_boundary(cued_success_Response.combosuccess_modulation_index<=-0.2 & cued_success_Response.combosuccess_modulation_index~=-1)=2;
    cued_success_Response.idx_diagonal_boundary=nan(size(cued_success_Response.idx));
    cued_success_Response.idx_diagonal_boundary(:)=2;
    cued_success_Response.idx_diagonal_boundary((cued_success_Response.combosuccess_modulation_index-cued_success_Response.combofailure_modulation_index)>0.2)=1;
    cued_success_Response.combo_boundary=nan(size(cued_success_Response.idx));
    cued_success_Response.combo_boundary(cued_success_Response.idx_from_boundary==1 & cued_success_Response.idx_diagonal_boundary==1)=1;
    cued_success_Response.combo_boundary(cued_success_Response.idx_from_boundary==2 & cued_success_Response.idx_diagonal_boundary==2)=2;   
    cued_success_Response.consensus_idx=nan(size(cued_success_Response.idx));
    cued_success_Response.consensus_idx(cued_success_Response.combo_boundary==1 | (cued_success_Response.idx==1 & cued_success_Response.isHighWeight==1) & (cued_success_Response.combo_boundary~=2 & cued_success_Response.idx~=2))=1; 
    cued_success_Response.consensus_idx(cued_success_Response.combo_boundary==2 | (cued_success_Response.idx==2 & cued_success_Response.isHighWeight==1) & (cued_success_Response.combo_boundary~=1 & cued_success_Response.idx~=1))=2;

    % Working on final boundary functions
%     figure(); scatter(cued_success_Response.allsuccess_modulation_index(cued_success_Response.consensus_idx==1),cued_success_Response.allfailure_modulation_index(cued_success_Response.consensus_idx==1),[],'b');
%     hold on; scatter(cued_success_Response.allsuccess_modulation_index(cued_success_Response.consensus_idx==2),cued_success_Response.allfailure_modulation_index(cued_success_Response.consensus_idx==2),[],'r');
%     hold on; x1=-pi+0.01:0.01:-0.01; plot(x1,-(0.2*csc(3.25*(x1+0.2))+0.7)); xlim([-1 1]); ylim([-1 1]);
%     figure(); scatter(cued_success_Response.cXsuccess_modulation_index(cued_success_Response.consensus_idx==1),cued_success_Response.cXfailure_modulation_index(cued_success_Response.consensus_idx==1),[],'b');
%     hold on; scatter(cued_success_Response.cXsuccess_modulation_index(cued_success_Response.consensus_idx==2),cued_success_Response.cXfailure_modulation_index(cued_success_Response.consensus_idx==2),[],'r');
%     x2=0.01:0.01:pi-0.01; hold on; plot(x2,-(0.1*csc(3*x2)-1.1)); xlim([-1 1]); ylim([-1 1]);

%     % THIS IS THE BOUNDARY FUNCTION:
%     % xie=-1:0.01:1; yie=-0.75*exp(-1.1*(xie-0.6))+1.5; hold all; plot(xie,yie);
%     % Use it to define cell types, based on glm coefs
%     ythresh=-0.75*exp(-1.1.*(cued_success_Response.combosuccess_modulation_index-0.6))+1.5; belowthresh1=cued_success_Response.combofailure_modulation_index<ythresh; abovethresh1=cued_success_Response.combofailure_modulation_index>=ythresh;
%     ythresh=-0.75*exp(-1.1.*(cued_success_Response.allsuccess_modulation_index-0.6))+1.5; belowthresh2=cued_success_Response.allfailure_modulation_index<ythresh; abovethresh2=cued_success_Response.allfailure_modulation_index>=ythresh;
%     ythresh=-0.75*exp(-1.1.*(cued_success_Response.cXsuccess_modulation_index-0.6))+1.5; belowthresh3=cued_success_Response.cXfailure_modulation_index<ythresh; abovethresh3=cued_success_Response.cXfailure_modulation_index>=ythresh;
%     cued_success_Response.idx=nan(size(cued_success_Response.idx)); 
% %     cued_success_Response.idx(belowthresh1==1 | belowthresh2==1 | belowthresh3==1 | (cued_success_Response.isHighWeight==1 & cued_success_Response.consensus_idx==0))=2; % 2 is failure-continuing
% %     cued_success_Response.idx((abovethresh1==1 & abovethresh2==1 & abovethresh3==1) | (cued_success_Response.isHighWeight==1 & cued_success_Response.consensus_idx==1) | ((abovethresh1==1 | abovethresh2==1 | abovethresh3==1) & ~(belowthresh1==1 | belowthresh2==1 | belowthresh3==1)))=1; % 1 is success-continuing
%     cued_success_Response.idx((belowthresh1==1 | belowthresh2==1 | belowthresh3==1))=2; % 2 is failure-continuing
%     cued_success_Response.idx(((abovethresh1==1 & abovethresh2==1 & abovethresh3==1) | ((abovethresh1==1 | abovethresh2==1 | abovethresh3==1) & ~(belowthresh1==1 | belowthresh2==1 | belowthresh3==1))))=1; % 1 is success-continuing
% %     cued_success_Response.idx((cued_success_Response.isHighWeight==1 & cued_success_Response.consensus_idx==0))=2; % 2 is failure-continuing
% %     cued_success_Response.idx((cued_success_Response.isHighWeight==1 & cued_success_Response.consensus_idx==1))=1; % 1 is success-continuing
%     cued_success_Response.idx(cued_success_Response.combosuccess_modulation_index==-1 & cued_success_Response.allsuccess_modulation_index==-1 & cued_success_Response.cXsuccess_modulation_index==-1)=nan;
%     idx=cued_success_Response.idx;

    % Define cued units
    cuedamount=defineCuedUnits(cued_success_Response);
    cued_success_Response.cuedamount=cuedamount;

    r{1}=cued_success_Response;
%     load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM test set\excluded trials where opto during cue\cued_failureNotDrop_Response.mat'); r{2}=cued_failureNotDrop_Response;
%     load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM test set\excluded trials where opto during cue\uncued_failureNotDrop_Response.mat'); r{3}=uncued_failureNotDrop_Response;
    load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM test set\excluded trials where opto during cue\cued_failure_Response.mat'); r{2}=cued_failure_Response;
    load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM test set\excluded trials where opto during cue\uncued_failure_Response.mat'); r{3}=uncued_failure_Response;
    load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM test set\excluded trials where opto during cue\uncued_success_Response.mat'); r{4}=uncued_success_Response;
    load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM test set\excluded trials where opto during cue\cue_Response.mat'); cue_Response=cue_Response; r{5}=cue_Response;
%     load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM test set\only trials where opto during cue\cued_drop_Response.mat'); r{5}=cued_drop_Response;
%     load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM test set\cued_failureNotDrop_Response.mat'); r{6}=cued_failureNotDrop_Response;
%     load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM test set\only trials where opto during cue\uncued_drop_Response.mat'); r{6}=uncued_drop_Response;
%     load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM test set\uncued_failureNotDrop_Response.mat'); r{8}=uncued_failureNotDrop_Response;
%     load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\all trials\cued_reach_Response.mat'); r{9}=cued_reach_Response;
%     load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\all trials\uncued_reach_Response.mat'); r{10}=uncued_reach_Response;
%     load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM test set\cued_failure_noReach_Response.mat'); r{11}=cued_failure_noReach_Response;
%     load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM test set\uncued_failure_noReach_Response.mat'); r{12}=uncued_failure_noReach_Response;
%     load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM test set\all_success_Response.mat'); r{13}=all_success_Response;
%     load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM test set\all_failure_Response.mat'); r{14}=all_failure_Response;

%     load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\test set\cued_success_Response.mat'); cued_success_Response.idx=idx; r{1}=cued_success_Response;
%     load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\test set\cued_failure_Response.mat'); cued_failure_Response.idx=idx; r{2}=cued_failure_Response;
%     load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\test set\uncued_failure_Response.mat'); uncued_failure_Response.idx=idx; r{3}=uncued_failure_Response;
%     load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\test set\uncued_success_Response.mat'); uncued_success_Response.idx=idx; r{4}=uncued_success_Response;
%     load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\test set\cued_drop_Response.mat'); r{5}=cued_drop_Response;
%     load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\test set\cued_failureNotDrop_Response.mat'); r{6}=cued_failureNotDrop_Response;
%     load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\test set\uncued_drop_Response.mat'); r{7}=uncued_drop_Response;
%     load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\test set\uncued_failureNotDrop_Response.mat'); r{8}=uncued_failureNotDrop_Response;
%     load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\all trials\cued_reach_Response.mat'); r{9}=cued_reach_Response;
%     load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\all trials\uncued_reach_Response.mat'); r{10}=uncued_reach_Response;
%     load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\test set\cued_failure_noReach_Response.mat'); r{11}=cued_failure_noReach_Response;
%     load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\test set\uncued_failure_noReach_Response.mat'); r{12}=uncued_failure_noReach_Response;
end
r=matchAllUnits(r);
cued_success_Response=r{1};
cued_failure_Response=r{2};
uncued_failure_Response=r{3};
uncued_success_Response=r{4};
cue_Response=r{5};
% cued_drop_Response=r{5};
% cued_failureNotDrop_Response=r{6};
% uncued_drop_Response=r{6};
% uncued_failureNotDrop_Response=r{8};
% cued_reach_Response=r{9};
% uncued_reach_Response=r{10};
% cued_failure_noReach_Response=r{11};
% uncued_failure_noReach_Response=r{12};
% all_success_Response=r{13};
% all_failure_Response=r{14};
% groupLabelsFromTCA=cued_success_Response.idx;

% Exclude non-SPN units, i.e., firing rate > 4 Hz
% nonSPNs=[762 797 1541]; 
% trmv=zeros(length(cued_success_Response.excluded),1); trmv(nonSPNs)=1; trmv=logical(trmv);
% cued_success_Response=removeUnitFromResponse(cued_success_Response,trmv);
% cued_failure_Response=removeUnitFromResponse(cued_failure_Response,trmv);
% uncued_failure_Response=removeUnitFromResponse(uncued_failure_Response,trmv);
% uncued_success_Response=removeUnitFromResponse(uncued_success_Response,trmv);
% cued_drop_Response=removeUnitFromResponse(cued_drop_Response,trmv);
% cued_failureNotDrop_Response=removeUnitFromResponse(cued_failureNotDrop_Response,trmv);
% uncued_drop_Response=removeUnitFromResponse(uncued_drop_Response,trmv);
% uncued_failureNotDrop_Response=removeUnitFromResponse(uncued_failureNotDrop_Response,trmv);
% cued_reach_Response=removeUnitFromResponse(cued_reach_Response,trmv);
% uncued_reach_Response=removeUnitFromResponse(uncued_reach_Response,trmv);
% cued_failure_noReach_Response=removeUnitFromResponse(cued_failure_noReach_Response,trmv);
% uncued_failure_noReach_Response=removeUnitFromResponse(uncued_failure_noReach_Response,trmv);

usingGLMidx=false;
if usingGLMidx==true
    % Remove units with firing rates too low??
    toolow=zeros(size(nanmean(cued_success_Response.unitbyunit_y,2)));
%     toolow=nanmean(cued_success_Response.unitbyunit_y,2)+nanmean(uncued_success_Response.unitbyunit_y,2)+nanmean(cued_failure_Response.unitbyunit_y,2)+nanmean(uncued_failure_Response.unitbyunit_y,2)<2; 
    
    % remove all units with nan classification
    f=find(cued_success_Response.excluded==0); 
    trmv=cued_success_Response.excluded;
    trmv(f(isnan(cued_success_Response.idx)))=1; 
    trmv(toolow==1)=1;
    trmv=logical(trmv);
    cued_success_Response=removeUnitFromResponse(cued_success_Response,trmv);
    cued_failure_Response=removeUnitFromResponse(cued_failure_Response,trmv);
    uncued_failure_Response=removeUnitFromResponse(uncued_failure_Response,trmv);
    uncued_success_Response=removeUnitFromResponse(uncued_success_Response,trmv);
%     cued_drop_Response=removeUnitFromResponse(cued_drop_Response,trmv);
%     cued_failureNotDrop_Response=removeUnitFromResponse(cued_failureNotDrop_Response,trmv);
%     uncued_drop_Response=removeUnitFromResponse(uncued_drop_Response,trmv);
%     uncued_failureNotDrop_Response=removeUnitFromResponse(uncued_failureNotDrop_Response,trmv);
%     cued_reach_Response=removeUnitFromResponse(cued_reach_Response,trmv);
%     uncued_reach_Response=removeUnitFromResponse(uncued_reach_Response,trmv);
%     cued_failure_noReach_Response=removeUnitFromResponse(cued_failure_noReach_Response,trmv);
%     uncued_failure_noReach_Response=removeUnitFromResponse(uncued_failure_noReach_Response,trmv);
%     all_success_Response=removeUnitFromResponse(all_success_Response,trmv);
%     all_failure_Response=removeUnitFromResponse(all_failure_Response,trmv);
end

% Average firing rates could try excluding all with trial_n_cutoff=5 or 8
% or 9 or 12 or more (but only for succ-continuing)
plotUnitSummariesAfterTCAlabels(cued_success_Response.consensus_idx,[],cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,[],'justAvs','justAvs');
% plotUnitSummariesAfterTCAlabels(cued_success_Response.idx,[],cued_success_Response,cued_drop_Response,uncued_success_Response,uncued_drop_Response,[],'justAvs','justAvs');
% plotUnitSummariesAfterTCAlabels(cued_success_Response.idx,[],cued_success_Response,cued_failureNotDrop_Response,uncued_success_Response,uncued_failureNotDrop_Response,[],'justAvs','justAvs');
% plotUnitSummariesAfterTCAlabels(cued_success_Response.idx,[],cued_success_Response,cued_failure_noReach_Response,uncued_success_Response,uncued_failure_noReach_Response,[],'justAvs','justAvs');

% fi=find(~isnan(cued_success_Response.carrythroughindsintowhich_units)); 
% i=1; 
% figure(); plot(cued_drop_Response.unitbyunit_x(fi(i),:),smoothdata(cued_drop_Response.unitbyunit_y(fi(i),:),'gaussian',10)); 
% hold on; plot(cued_drop_Response.aligncomp_x(fi(i),:),cued_drop_Response.aligncomp_y(fi(i),:),'Color','b');
% disp(['Doing ' names_of_units{i} ' from session ' num2str(sessions_of_units(i))]);

%%%%%%%%%%%%% Inspect units one by one
% plotGLMcoef(consensus_glm_coef(f(156),:),0,feature_names,10*0.01,nansum(shifts<0),'mean',false,[]); title('python glm');
% i=find(ismember(names_of_units,unitnames_glm{f(156)}));
% figure(); plot(uncued_drop_Response.unitbyunit_x(i,:),smoothdata(uncued_drop_Response.unitbyunit_y(i,:),'gaussian',10),'Color','r');
% hold on; plot(uncued_drop_Response.aligncomp_x(i,:),uncued_drop_Response.aligncomp_y(i,:),'Color','b');
% figure(); plot(uncued_success_Response.unitbyunit_x(i,:),smoothdata(uncued_success_Response.unitbyunit_y(i,:),'gaussian',10),'Color','g');
% hold on; plot(uncued_success_Response.aligncomp_x(i,:),uncued_success_Response.aligncomp_y(i,:),'Color','b');
% disp(['Doing ' names_of_units{i} ' from session ' num2str(sessions_of_units(i))]);

%%%%%%%%%%%% Plot cue response of subpopulation
% f=find(cued_success_Response.excluded==0); 
% trmv=cued_success_Response.excluded;
% trmv(f(~(cued_success_Response.idx==1 & cued_success_Response.isHighWeight==1)))=1;
% trmv=logical(trmv);
% subCueResponse=removeUnitFromResponse(cue_Response,trmv);
% plotVariousSUResponsesAlignedToBeh('meanAcrossUnits',subCueResponse,1);

%% TUNING OF PERSISTENT ACTIVITY
clear r
% load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\training\TCA\idx_groupLabelsFromTCA.mat');
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\test set\cued_success_Response.mat'); cued_success_Response.idx=idx; r{1}=cued_success_Response;
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\test set\cued_failure_Response.mat'); cued_failure_Response.idx=cued_success_Response.idx; r{2}=cued_failure_Response;
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\test set\uncued_failure_Response.mat'); uncued_failure_Response.idx=idx; r{3}=uncued_failure_Response;
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\test set\uncued_success_Response.mat'); uncued_success_Response.idx=idx; r{4}=uncued_success_Response;
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\all trials\cued_reach_Response.mat'); r{5}=cued_reach_Response;
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\all trials\uncued_reach_Response.mat'); r{6}=uncued_reach_Response;
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\test set\cued_drop_Response.mat'); r{7}=cued_drop_Response;
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\test set\uncued_drop_Response.mat'); r{8}=uncued_drop_Response;
r=matchAllUnits(r);

cued_success_Response=r{1};
cued_failure_Response=r{2};
uncued_failure_Response=r{3};
uncued_success_Response=r{4};
cued_reach_Response=r{5};
uncued_reach_Response=r{6};
cued_drop_Response=r{7};
uncued_drop_Response=r{8};
groupLabelsFromTCA=cued_success_Response.idx;

% Exclude non-SPN units, i.e., firing rate > 4 Hz
% nonSPNs=[762 797 1541]; trmv=zeros(length(cued_success_Response.excluded),1); trmv(nonSPNs)=1; trmv=logical(trmv);
% cued_success_Response=removeUnitFromResponse(cued_success_Response,trmv);
% cued_failure_Response=removeUnitFromResponse(cued_failure_Response,trmv);
% uncued_failure_Response=removeUnitFromResponse(uncued_failure_Response,trmv);
% uncued_success_Response=removeUnitFromResponse(uncued_success_Response,trmv);
% cued_reach_Response=removeUnitFromResponse(cued_reach_Response,trmv);
% uncued_reach_Response=removeUnitFromResponse(uncued_reach_Response,trmv);
% cued_drop_Response=removeUnitFromResponse(cued_drop_Response,trmv);
% uncued_drop_Response=removeUnitFromResponse(uncued_drop_Response,trmv);

% CUED
% cuedReach_Response=cued_success_Response;
% cuedReach_Response.unitbyunit_y=(cued_success_Response.unitbyunit_y+cued_failure_Response.unitbyunit_y)./2;
% cuez=getCueTunedUnits(cued_reach_Response,uncued_reach_Response,'cue_vs_baseline_no_index','mean',1,[7 16],[-2 0],[7 16],[-2 0]); % method 3rd arg can be 'vs_uncued_reach' or 'cue_vs_baseline' or 'justcue'
% cuez=getCueTunedUnits(cued_reach_Response,uncued_reach_Response,'cue_vs_baseline_no_index','mean',1,[7 16],[-1 0.5],[7 16],[-1 0.5]); % method 3rd arg can be 'vs_uncued_reach' or 'cue_vs_baseline' or 'justcue'
%%% for the sustained reach-related activity, take cells that turn on
%%% within 1 sec window of the arm outstretched
% grp 2
cuez=getCueTunedUnits(cued_reach_Response,uncued_reach_Response,'cue_vs_baseline_no_index','mean',1,[7 16],[-0.1 0],[7 16],[-1 0.5]); % method 3rd arg can be 'vs_uncued_reach' or 'cue_vs_baseline' or 'justcue'
plotUnitSummariesAfterTCAlabels(cued_success_Response.idx,cuez,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,[],'cued','tuning');
% grp 1
cuez=getCueTunedUnits(cued_reach_Response,uncued_reach_Response,'cue_vs_baseline_no_index','mean',1,[7 16],[-1 0],[7 16],[-1 0.5]); % method 3rd arg can be 'vs_uncued_reach' or 'cue_vs_baseline' or 'justcue'
plotUnitSummariesAfterTCAlabels(cued_success_Response.idx,cuez,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,[],'cued','tuning');
cuez=getCueTunedUnits(cued_reach_Response,uncued_reach_Response,'cue_vs_baseline_no_index','mean',1,[7 16],[-1 -0.5],[7 16],[-1 0.5]); % method 3rd arg can be 'vs_uncued_reach' or 'cue_vs_baseline' or 'justcue'
plotUnitSummariesAfterTCAlabels(cued_success_Response.idx,cuez,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,[],'cued','tuning');

% UNCUED
% cuez=getCueTunedUnits(uncued_reach_Response,uncuedReach_Response,'cue_vs_baseline_no_index','mean',1,[7 16],[-2 0],[7 16],[-2 0]); % method 3rd arg can be 'vs_uncued_reach' or 'cue_vs_baseline' or 'justcue'
% plotUnitSummariesAfterTCAlabels(groupLabelsFromTCA,cuez,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,[],'uncued');

% Plot DIFFERENT IN CUED V UNCUED
% get rid of cells that NEVER turn on -- they are not useful for this
% comparison
trial_n_cutoff=10;
out=plotVariousSUResponsesAlignedToBeh('matchUnitsAcrossResponses',excludeTooFewTrials(cued_reach_Response,trial_n_cutoff,false),excludeTooFewTrials(cued_success_Response,trial_n_cutoff,false),excludeTooFewTrials(cued_failure_Response,trial_n_cutoff,false),excludeTooFewTrials(uncued_success_Response,trial_n_cutoff,false),excludeTooFewTrials(uncued_failure_Response,trial_n_cutoff,false));
cued_reach_Response=out.Response1; cued_success_Response=out.Response2; cued_failure_Response=out.Response3; uncued_success_Response=out.Response4; uncued_failure_Response=out.Response5;
r{1}=cued_reach_Response; r{2}=uncued_reach_Response; r=matchAllUnits(r); uncued_reach_Response=r{2};
doesntSpike=~any([downSampMatrix(cued_success_Response.unitbyunit_y(:,201:500),100) downSampMatrix(cued_failure_Response.unitbyunit_y(:,1:2000),100) downSampMatrix(uncued_success_Response.unitbyunit_y(:,1:2000),100) downSampMatrix(uncued_failure_Response.unitbyunit_y(:,1:2000),100)]>1,2); 
newexcl=cued_success_Response.excluded; newexcl(doesntSpike)=1; newexcl=logical(newexcl);
cued_success_Response=removeUnitFromResponse(cued_success_Response,newexcl);
cued_failure_Response=removeUnitFromResponse(cued_failure_Response,newexcl);
uncued_success_Response=removeUnitFromResponse(uncued_success_Response,newexcl);
uncued_failure_Response=removeUnitFromResponse(uncued_failure_Response,newexcl);
cued_reach_Response=removeUnitFromResponse(cued_reach_Response,newexcl);
uncued_reach_Response=removeUnitFromResponse(uncued_reach_Response,newexcl);
cued_drop_Response=removeUnitFromResponse(cued_drop_Response,newexcl);
uncued_drop_Response=removeUnitFromResponse(uncued_drop_Response,newexcl);
% cuez=getCueTunedUnits(cued_reach_Response,uncued_reach_Response,'vs_uncued_reach_no_index','mean',1,{[-5 -4],[9 12.5]},[-3.3 0],{[-5 -4],[9 12.5]},[-3.3 0]); % method 3rd arg can be 'vs_uncued_reach' or 'cue_vs_baseline' or 'justcue' or 'vs_uncued_reach_no_index'
% cuez=getCueTunedUnits(cued_reach_Response,uncued_reach_Response,'justcue_v_justuncue','mean',1,[4 12],[-3.3 0.5],[4 12],[-3.3 0.5]); % method 3rd arg can be 'vs_uncued_reach' or 'cue_vs_baseline' or 'justcue' or 'vs_uncued_reach_no_index'
% cuez=getCueTunedUnits(cued_reach_Response,uncued_reach_Response,'justcue_v_justuncue','mean',1,[4 12.5],[-2 0],[4 12.5],[-2 0]); % method 3rd arg can be 'vs_uncued_reach' or 'cue_vs_baseline' or 'justcue' or 'vs_uncued_reach_no_index'
% cuez=getCueTunedUnits(cue_noReach_Response,uncued_reach_Response,'vs_uncued_reach_no_index','mean',1,[9 12.5],[-0.37 1.5],[9 12.5],[-2 0]); % method 3rd arg can be 'vs_uncued_reach' or 'cue_vs_baseline' or 'justcue' or 'vs_uncued_reach_no_index'
% cuez=getCueTunedUnits(cued_reach_Response,uncued_reach_Response,'justcue_v_justuncue','mean',1,[7 16],[-0.1 0],[7 16],[-0.1 0]);
% cuez=getCueTunedUnits(cued_reach_Response,uncued_reach_Response,'justcue_v_justuncue','mean',1,[7 16],[-0.1 0],[7 16],[-0.1 0]);
cuez=getCueTunedUnits(cued_reach_Response,uncued_reach_Response,'justcue_v_justuncue','mean',1,[7 16],[-1 -0.2],[7 16],[-1 -0.2]);
plotUnitSummariesAfterTCAlabels(cued_success_Response.idx,cuez,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,[],'uncuedOverCued','tuning');

plotUnitSummariesAfterTCAlabels(cued_success_Response.idx,cued_success_Response.cXmiss_sustained-cued_success_Response.allMiss_sustained,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,[],'uncuedOverCued','tuning');
plotUnitSummariesAfterTCAlabels(cued_success_Response.idx,cued_success_Response.postCueAmp_over1sec,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,[],'uncuedOverCued','tuning');
plotUnitSummariesAfterTCAlabels(cued_success_Response.idx,cued_success_Response.cuecoef_over1sec,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,[],'uncuedOverCued','tuning');

plotUnitSummariesAfterTCAlabels(cued_success_Response.idx,cued_success_Response.allSucc_sustained+cued_success_Response.cXsucc_sustained,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,[],'uncuedOverCued','tuning');
plotUnitSummariesAfterTCAlabels(cued_success_Response.idx,cued_success_Response.cXsucc_over1sec-cued_success_Response.allSucc_over1sec,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,[],'uncuedOverCued','tuning');
plotUnitSummariesAfterTCAlabels(cued_success_Response.idx,cued_success_Response.cuecoef_over1sec,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,[],'uncuedOverCued','tuning');

plotUnitSummariesAfterTCAlabels(cued_success_Response.idx,(cued_success_Response.glmcoef_index21+cued_success_Response.glmcoef_index22+cued_success_Response.glmcoef_index23+cued_success_Response.glmcoef_index24+cued_success_Response.glmcoef_index25+cued_success_Response.glmcoef_index26),cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,[],'uncuedOverCued','tuning');
plotUnitSummariesAfterTCAlabels(cued_success_Response.idx,(cued_success_Response.glmcoef_index21+cued_success_Response.glmcoef_index22+cued_success_Response.glmcoef_index23+cued_success_Response.glmcoef_index24+cued_success_Response.glmcoef_index25+cued_success_Response.glmcoef_index26)-...
                                                          (cued_success_Response.glmcoef_index20+cued_success_Response.glmcoef_index19+cued_success_Response.glmcoef_index18+cued_success_Response.glmcoef_index17+cued_success_Response.glmcoef_index16+cued_success_Response.glmcoef_index16),cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,[],'uncuedOverCued','tuning');

plotUnitSummariesAfterTCAlabels(cued_success_Response.idx,cued_success_Response.isHighWeight,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,[],'uncuedOverCued','tuning');

plotUnitSummariesAfterTCAlabels(cued_success_Response.idx,cued_success_Response.cXfail_sustained,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,[],'uncuedOverCued','tuning');

% % 20230914
% FOR CUE TUNING SPREAD OUT
binsForTuning{1}=[-10 -0.01 10]; binsForTuning{2}=[-10 -0.01 10];
tuningOutput=plotUnitSummariesAfterTCAlabels(cued_success_Response.consensus_idx,cued_success_Response.cXfail_sus1to5sec-cued_success_Response.allfail_sus1to5sec,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,[],'uncuedOverCued','tuning',binsForTuning);
plotTuningOutputScatter(tuningOutput,'grp1_fail','grp1_fail_uncue',1,[2 5]);
binsForTuning{1}=[-0.01 -0.0001 10]; binsForTuning{2}=[-0.01 -0.0001 10];
tuningOutput=plotUnitSummariesAfterTCAlabels(cued_success_Response.consensus_idx,cued_success_Response.cXfail_sus1to5sec-cued_success_Response.allfail_sus1to5sec,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,[],'uncuedOverCued','tuning',binsForTuning);
plotTuningOutputScatter(tuningOutput,'grp1_fail','grp1_fail_uncue',1,[2 5]);
binsForTuning{1}=[-0.0001 0.0001 10]; binsForTuning{2}=[-0.0001 0.0001 10];
tuningOutput=plotUnitSummariesAfterTCAlabels(cued_success_Response.consensus_idx,cued_success_Response.cXfail_sus1to5sec-cued_success_Response.allfail_sus1to5sec,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,[],'uncuedOverCued','tuning',binsForTuning);
plotTuningOutputScatter(tuningOutput,'grp1_fail','grp1_fail_uncue',1,[2 5]);
binsForTuning{1}=[-0.0001 0.0001 10]; binsForTuning{2}=[-0.0001 0.0001 10];
tuningOutput=plotUnitSummariesAfterTCAlabels(cued_success_Response.consensus_idx,cued_success_Response.cXfail_sus1to5sec-cued_success_Response.allfail_sus1to5sec,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,[],'uncuedOverCued','tuning',binsForTuning);
plotTuningOutputScatter(tuningOutput,'grp1_fail','grp1_fail_uncue',2,[2 5]);

% binsForTuning{1}=[-10 -0.0001 10]; binsForTuning{2}=[-10 -0.0001 10];
% tuningOutput=plotUnitSummariesAfterTCAlabels(cued_success_Response.consensus_idx,cued_success_Response.cXfail_sus1to5sec-cued_success_Response.allfail_sus1to5sec,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,[],'uncuedOverCued','tuning',binsForTuning);
binsForTuning{1}=[-10 0.5 10]; binsForTuning{2}=[-10 0.5 10];
tuningOutput=plotUnitSummariesAfterTCAlabels(cued_success_Response.consensus_idx,cued_success_Response.cuecoef_over1sec,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,[],'uncuedOverCued','tuning',binsForTuning);
[allgp1_cuedfailFR_cuedir2,allgp1_uncuedfailFR_cuedir2]=plotTuningOutputScatter(tuningOutput,'grp1_fail','grp1_fail_uncue',2,[2 5]);
[allgp1_cuedfailFR_cuedir1,allgp1_uncuedfailFR_cuedir1]=plotTuningOutputScatter(tuningOutput,'grp1_fail','grp1_fail_uncue',1,[2 5]);
[allgp1_cuedsuccFR_cuedir2,allgp1_uncuedsuccFR_cuedir2]=plotTuningOutputScatter(tuningOutput,'grp1_succ','grp1_succ_uncue',2,[2 5]);
[allgp1_cuedsuccFR_cuedir1,allgp1_uncuedsuccFR_cuedir1]=plotTuningOutputScatter(tuningOutput,'grp1_succ','grp1_succ_uncue',1,[2 5]);
binsForTuning{1}=[-10 0 10]; binsForTuning{2}=[-10 0 10];
tuningOutput=plotUnitSummariesAfterTCAlabels(cued_success_Response.consensus_idx,cued_success_Response.cXsucc_sus1to5sec,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,[],'uncuedOverCued','tuning',binsForTuning);
[allgp2_cuedfailFR_cuedir2,allgp2_uncuedfailFR_cuedir2]=plotTuningOutputScatter(tuningOutput,'grp2_fail','grp2_fail_uncue',2,[2 5]);
[allgp2_cuedfailFR_cuedir1,allgp2_uncuedfailFR_cuedir1]=plotTuningOutputScatter(tuningOutput,'grp2_fail','grp2_fail_uncue',1,[2 5]);
[allgp2_cuedsuccFR_cuedir2,allgp2_uncuedsuccFR_cuedir2]=plotTuningOutputScatter(tuningOutput,'grp2_succ','grp2_succ_uncue',2,[2 5]);
[allgp2_cuedsuccFR_cuedir1,allgp2_uncuedsuccFR_cuedir1]=plotTuningOutputScatter(tuningOutput,'grp2_succ','grp2_succ_uncue',1,[2 5]);
% % plotUnitSummariesAfterTCAlabels(cued_success_Response.consensus_idx,cued_success_Response.cXfail_sus1to5sec-cued_success_Response.allfail_sus1to5sec,cued_success_Response,cued_failure_Response,uncued_success_Response,cue_Response,[],'uncuedOverCued','tuning',binsForTuning);
% % or consider just showing stats for all units
% binsForTuning{1}=[-10 -9 10]; binsForTuning{2}=[-10 -9 10];
% % binsForTuning{1}=[-10 -0.00005 10]; binsForTuning{2}=[-10 -0.00005 10];
% tuningOutput=plotUnitSummariesAfterTCAlabels(cued_success_Response.consensus_idx,cued_success_Response.cXfail_sus1to5sec,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,[],'uncuedOverCued','tuning',binsForTuning);
% [allgp1_cuedfailFR,allgp1_uncuedfailFR]=plotTuningOutputScatter(tuningOutput,'grp1_fail','grp1_fail_uncue',2,[2 5]);
% [allgp2_cuedfailFR,allgp2_uncuedfailFR]=plotTuningOutputScatter(tuningOutput,'grp2_fail','grp2_fail_uncue',2,[2 5]);
% plotUnitSummariesAfterTCAlabels(cued_success_Response.consensus_idx,cued_success_Response.cXfail_sus1to5sec,cued_success_Response,cued_failure_Response,uncued_success_Response,cue_Response,[],'uncuedOverCued','tuning',binsForTuning);

% all units
binsForTuning{1}=[-10 -9 10]; binsForTuning{2}=[-10 -9 10];
tuningOutput=plotUnitSummariesAfterTCAlabels(cued_success_Response.consensus_idx,cued_success_Response.cXfail_sus1to5sec,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,[],'uncuedOverCued','tuning',binsForTuning);
[allgp1_cuedfailFR,allgp1_uncuedfailFR]=plotTuningOutputScatter(tuningOutput,'grp1_fail','grp1_fail_uncue',2,[2 5]);
[allgp2_cuedfailFR,allgp2_uncuedfailFR]=plotTuningOutputScatter(tuningOutput,'grp2_fail','grp2_fail_uncue',2,[2 5]);
plotTuningOutputScatter(tuningOutput,'grp1_succ_uncue','grp1_fail_uncue',2,[2 5]);
[allgp1_cuedsuccFR,allgp1_uncuedsuccFR]=plotTuningOutputScatter(tuningOutput,'grp1_succ','grp1_succ_uncue',2,[2 5]);
[allgp2_cuedsuccFR,allgp2_uncuedsuccFR]=plotTuningOutputScatter(tuningOutput,'grp2_succ','grp2_succ_uncue',2,[2 5]);
plotTuningOutputScatter(tuningOutput,'grp2_succ','grp2_fail',2,[2 5]);

% % binsForTuning{1}=[-10 0 10]; binsForTuning{2}=[-10 0 10];
% binsForTuning{1}=[-10 -9 10]; binsForTuning{2}=[-10 -9 10];
% tuningOutput=plotUnitSummariesAfterTCAlabels(cued_success_Response.consensus_idx,cued_success_Response.cXsucc_sus0to5sec,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,[],'uncuedOverCued','tuning',binsForTuning); % first n is light green (lower range), second n is dark red (higher range)
% [allgp2_cuedsuccFR,allgp2_uncuedsuccFR]=plotTuningOutputScatter(tuningOutput,'grp2_succ','grp2_succ_uncue',2,[2 5]);
% [allgp1_cuedsuccFR,allgp1_uncuedsuccFR]=plotTuningOutputScatter(tuningOutput,'grp1_succ','grp1_succ_uncue',2,[2 5]);
% 
% % sort by cued
% binsForTuning{1}=[-10 1 10]; binsForTuning{2}=[-10 1 10];
% tuningOutput=plotUnitSummariesAfterTCAlabels(cued_success_Response.consensus_idx,cued_success_Response.cuecoef_over1sec,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,[],'uncuedOverCued','tuning',binsForTuning);
% % cuez=getCueTunedUnits(uncuedReach_Response,cuedReach_Response,'justcue_v_justuncue','mean',1,[4 12],[-2 0],[4 12],[-2 0]); 
% % plotUnitSummariesAfterTCAlabels(groupLabelsFromTCA,cuez,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,[],'uncuedOverCued');

% not randomly mixed selectivity
figure(); 
yaxis=abs(cued_success_Response.cXsucc_sus0to5sec-cued_success_Response.allsucc_sus0to5sec)-abs(cued_success_Response.cXfail_sus1to5sec-cued_success_Response.allfail_sus1to5sec);
xaxis=(cued_success_Response.cXsucc_sus0to5sec+cued_success_Response.allsucc_sus0to5sec)-(cued_success_Response.cXfail_sus1to5sec+cued_success_Response.allfail_sus1to5sec);
scatter(xaxis,yaxis,[],'k'); hold on; 
scatter(xaxis(cued_success_Response.consensus_idx==1),yaxis(cued_success_Response.consensus_idx==1),[],'b');
scatter(xaxis(cued_success_Response.consensus_idx==2),yaxis(cued_success_Response.consensus_idx==2),[],'r');
% before outcome
% figure(); 
% yaxis=abs(cued_success_Response.cXsucc_minus2to0sec-cued_success_Response.allsucc_minus2to0sec)-abs(cued_success_Response.cXfail_minus2to0sec-cued_success_Response.allfail_minus2to0sec);
% xaxis=(cued_success_Response.cXsucc_minus2to0sec+cued_success_Response.allsucc_minus2to0sec)-(cued_success_Response.cXfail_minus2to0sec+cued_success_Response.allfail_minus2to0sec);
% scatter(xaxis,yaxis,[],'k'); hold on; 
% scatter(xaxis(cued_success_Response.consensus_idx==1),yaxis(cued_success_Response.consensus_idx==1),[],'b');
% scatter(xaxis(cued_success_Response.consensus_idx==2),yaxis(cued_success_Response.consensus_idx==2),[],'r');
% shuffle fail wrt succ coeffs
yaxis_part1=abs(cued_success_Response.cXsucc_sus0to5sec-cued_success_Response.allsucc_sus0to5sec);
yaxis_part2=abs(cued_success_Response.cXfail_sus1to5sec-cued_success_Response.allfail_sus1to5sec);
xaxis_part1=(cued_success_Response.cXsucc_sus0to5sec+cued_success_Response.allsucc_sus0to5sec);
xaxis_part2=(cued_success_Response.cXfail_sus1to5sec+cued_success_Response.allfail_sus1to5sec);
yaxis_part2=yaxis_part2(randperm(length(yaxis_part2)));
xaxis_part2=xaxis_part2(randperm(length(xaxis_part2)));
yaxis_shuff=yaxis_part1-yaxis_part2;
xaxis_shuff=xaxis_part1-xaxis_part2;
% relative to cue
figure(); 
yaxis=abs(cued_success_Response.cXsucc_sus0to5sec-cued_success_Response.allsucc_sus0to5sec)-abs(cued_success_Response.cXfail_sus1to5sec-cued_success_Response.allfail_sus1to5sec);
xaxis=cued_success_Response.cuecoef_over1sec;
scatter(xaxis,yaxis,[],'k'); hold on; 
scatter(xaxis(cued_success_Response.consensus_idx==1),yaxis(cued_success_Response.consensus_idx==1),[],'b');
scatter(xaxis(cued_success_Response.consensus_idx==2),yaxis(cued_success_Response.consensus_idx==2),[],'r');

%% decode trial type
% axis x is activity of gp1 units (use window 2 to 5 sec)
% axis y is activity of gp2 units (use window 2 to 5 sec)
figure(); nBoot=100; nUnits=85;
takeThese_gp1=nan(nBoot,nUnits); takeThese_gp2=nan(nBoot,nUnits);
for i=1:nBoot
    takeThese_gp1(i,:)=randsample(length(allgp1_cuedsuccFR),nUnits); takeThese_gp2(i,:)=randsample(length(allgp2_cuedsuccFR),nUnits);
end
temp1=nan(nBoot,nUnits); temp2=nan(nBoot,nUnits);
for i=1:nBoot
    temp1(i,:)=nanmean(allgp1_cuedsuccFR(takeThese_gp1(i,:))); temp2(i,:)=nanmean(allgp2_cuedsuccFR(takeThese_gp2(i,:)));
end
scatter(temp1,temp2,[],'g'); hold on; scatter(nanmean(temp1),nanmean(temp2),[],'g','filled'); cuedsuccmeanx=nanmean(temp1); cuedsuccmeany=nanmean(temp2);
temp1=nan(nBoot,nUnits); temp2=nan(nBoot,nUnits);
for i=1:nBoot
    temp1(i,:)=nanmean(allgp1_cuedfailFR(takeThese_gp1(i,:))); temp2(i,:)=nanmean(allgp2_cuedfailFR(takeThese_gp2(i,:)));
end
scatter(temp1,temp2,[],'r'); scatter(nanmean(temp1),nanmean(temp2),[],'r','filled'); cuedfailmeanx=nanmean(temp1); cuedfailmeany=nanmean(temp2);
temp1=nan(nBoot,nUnits); temp2=nan(nBoot,nUnits);
for i=1:nBoot
    temp1(i,:)=nanmean(allgp1_uncuedsuccFR(takeThese_gp1(i,:))); temp2(i,:)=nanmean(allgp2_uncuedsuccFR(takeThese_gp2(i,:)));
end
scatter(temp1,temp2,[],'b'); scatter(nanmean(temp1),nanmean(temp2),[],'b','filled'); uncuedsuccmeanx=nanmean(temp1); uncuedsuccmeany=nanmean(temp2);
temp1=nan(nBoot,nUnits); temp2=nan(nBoot,nUnits);
for i=1:nBoot
    temp1(i,:)=nanmean(allgp1_uncuedfailFR(takeThese_gp1(i,:))); temp2(i,:)=nanmean(allgp2_uncuedfailFR(takeThese_gp2(i,:)));
end
scatter(temp1,temp2,[],'y'); scatter(nanmean(temp1),nanmean(temp2),[],'y','filled'); uncuedfailmeanx=nanmean(temp1); uncuedfailmeany=nanmean(temp2);
scatter((cuedsuccmeanx+cuedfailmeanx+uncuedsuccmeanx+uncuedfailmeanx)/4,(cuedsuccmeany+cuedfailmeany+uncuedsuccmeany+uncuedfailmeany)/4,[],'k','filled');
xlabel('Gp 1 average unit firing rate'); ylabel('Gp 2 average unit firing rate');

%% Trying a different mapping
figure(); nBoot=100; nUnits=40;
takeThese_gp1_cuedir1=nan(nBoot,nUnits); takeThese_gp1_cuedir2=nan(nBoot,nUnits); takeThese_gp2_cuedir1=nan(nBoot,nUnits); takeThese_gp2_cuedir2=nan(nBoot,nUnits);
for i=1:nBoot
    takeThese_gp1_cuedir1(i,:)=randsample(length(allgp1_cuedsuccFR_cuedir1),nUnits); 
    takeThese_gp1_cuedir2(i,:)=randsample(length(allgp1_cuedsuccFR_cuedir2),nUnits); 
    takeThese_gp2_cuedir1(i,:)=randsample(length(allgp2_cuedsuccFR_cuedir1),nUnits); 
    takeThese_gp2_cuedir2(i,:)=randsample(length(allgp2_cuedsuccFR_cuedir2),nUnits); 
end
temp1=nan(nBoot,nUnits); temp2=nan(nBoot,nUnits);
for i=1:nBoot
    temp1(i,:)=-nanmean(allgp1_cuedsuccFR_cuedir2(takeThese_gp1_cuedir2(i,:)))+nanmean(allgp1_cuedsuccFR_cuedir1(takeThese_gp1_cuedir1(i,:))); 
    temp2(i,:)=nanmean(allgp2_cuedsuccFR_cuedir2(takeThese_gp2_cuedir2(i,:)))-nanmean(allgp2_cuedsuccFR_cuedir1(takeThese_gp2_cuedir1(i,:))); 
end
scatter(temp1,temp2,[],'g'); hold on; scatter(nanmean(temp1),nanmean(temp2),[],'g','filled'); cuedsuccmeanx=nanmean(temp1); cuedsuccmeany=nanmean(temp2);
temp1=nan(nBoot,nUnits); temp2=nan(nBoot,nUnits);
for i=1:nBoot
    temp1(i,:)=-nanmean(allgp1_cuedfailFR_cuedir2(takeThese_gp1_cuedir2(i,:)))+nanmean(allgp1_cuedfailFR_cuedir1(takeThese_gp1_cuedir1(i,:))); 
    temp2(i,:)=nanmean(allgp2_cuedfailFR_cuedir2(takeThese_gp2_cuedir2(i,:)))-nanmean(allgp2_cuedfailFR_cuedir1(takeThese_gp2_cuedir1(i,:))); 
end
scatter(temp1,temp2,[],'r'); scatter(nanmean(temp1),nanmean(temp2),[],'r','filled'); cuedfailmeanx=nanmean(temp1); cuedfailmeany=nanmean(temp2);
temp1=nan(nBoot,nUnits); temp2=nan(nBoot,nUnits);
for i=1:nBoot
    temp1(i,:)=-nanmean(allgp1_uncuedsuccFR_cuedir2(takeThese_gp1_cuedir2(i,:)))+nanmean(allgp1_uncuedsuccFR_cuedir1(takeThese_gp1_cuedir1(i,:))); 
    temp2(i,:)=nanmean(allgp2_uncuedsuccFR_cuedir2(takeThese_gp2_cuedir2(i,:)))-nanmean(allgp2_uncuedsuccFR_cuedir1(takeThese_gp2_cuedir1(i,:)));
end
scatter(temp1,temp2,[],'b'); scatter(nanmean(temp1),nanmean(temp2),[],'b','filled'); uncuedsuccmeanx=nanmean(temp1); uncuedsuccmeany=nanmean(temp2);
temp1=nan(nBoot,nUnits); temp2=nan(nBoot,nUnits);
for i=1:nBoot
    temp1(i,:)=-nanmean(allgp1_uncuedfailFR_cuedir2(takeThese_gp1_cuedir2(i,:)))+nanmean(allgp1_uncuedfailFR_cuedir1(takeThese_gp1_cuedir1(i,:))); 
    temp2(i,:)=nanmean(allgp2_uncuedfailFR_cuedir2(takeThese_gp2_cuedir2(i,:)))-nanmean(allgp2_uncuedfailFR_cuedir1(takeThese_gp2_cuedir1(i,:)));
end
scatter(temp1,temp2,[],'y'); scatter(nanmean(temp1),nanmean(temp2),[],'y','filled'); uncuedfailmeanx=nanmean(temp1); uncuedfailmeany=nanmean(temp2);
scatter((cuedsuccmeanx+cuedfailmeanx+uncuedsuccmeanx+uncuedfailmeanx)/4,(cuedsuccmeany+cuedfailmeany+uncuedsuccmeany+uncuedfailmeany)/4,[],'k','filled');
xlabel('Gp 1 cue vs uncue'); ylabel('Gp 2 cue vs uncue');

%% CUED TUNING FROM GLM COEFFS
clear r
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\training\TCA\idx_groupLabelsFromTCA.mat');
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\test set\cued_success_Response.mat'); cued_success_Response.idx=idx; r{1}=cued_success_Response;
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\test set\cued_drop_Response.mat'); cued_drop_Response.idx=idx; r{2}=cued_drop_Response;
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\test set\uncued_success_Response.mat'); uncued_success_Response.idx=idx; r{3}=uncued_success_Response;
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\test set\uncued_drop_Response.mat'); uncued_drop_Response.idx=idx; r{4}=uncued_drop_Response;
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\all trials\cued_reach_Response.mat'); r{5}=cued_reach_Response;
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\all trials\uncued_reach_Response.mat'); r{6}=uncued_reach_Response;
r=matchAllUnits(r);
cued_success_Response=r{1};
cued_drop_Response=r{2};
uncued_success_Response=r{3};
uncued_drop_Response=r{4};
cued_reach_Response=r{5};
uncued_reach_Response=r{6};
groupLabelsFromTCA=cued_success_Response.idx;
% Exclude non-SPN units, i.e., firing rate > 4 Hz
nonSPNs=[762 797 1541]; trmv=zeros(length(cued_success_Response.excluded),1); trmv(nonSPNs)=1; trmv=logical(trmv);
cued_success_Response=removeUnitFromResponse(cued_success_Response,trmv);
cued_drop_Response=removeUnitFromResponse(cued_drop_Response,trmv);
uncued_success_Response=removeUnitFromResponse(uncued_success_Response,trmv);
uncued_drop_Response=removeUnitFromResponse(uncued_drop_Response,trmv);
cued_reach_Response=removeUnitFromResponse(cued_reach_Response,trmv);
uncued_reach_Response=removeUnitFromResponse(uncued_reach_Response,trmv);
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\matlab glm\unitnames_glm.mat');
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\matlab glm\all_glm_coef.mat');
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\matlab glm\fromWhichSess_glm.mat');
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\matlab glm\unitbyunit_names.mat');
[indexGLMcellsIntoUnitNames,indexUnitNamesIntoGLMcells,unitnames_glm_notInThisList]=getNamesIndexIntoNamesList(unitnames_glm,unitbyunit_names);

% line up Response and glm units
excluInds_unitbyunit_names=find(unitbyunit_names.excluded==0);
excluInds_cued_success_Response=find(cued_success_Response.excluded==0);
idxInExcludedFormat=nan(size(cued_success_Response.excluded)); 
idxInExcludedFormat(excluInds_cued_success_Response)=cued_success_Response.idx;
idxForUnitbyunitNames=idxInExcludedFormat(excluInds_unitbyunit_names);
idxForGLM=idxForUnitbyunitNames(indexUnitNamesIntoGLMcells(~isnan(indexUnitNamesIntoGLMcells(:,1)),1));
% and bring D1 and A2a tags
tagsInExcludedFormat=nan(size(cued_success_Response.excluded)); 
tagsInExcludedFormat(excluInds_cued_success_Response)=cued_success_Response.D1tag(cued_success_Response.excluded==0);
tagsForUnitbyunitNames=tagsInExcludedFormat(excluInds_unitbyunit_names);
D1tagForGLM=tagsForUnitbyunitNames(indexUnitNamesIntoGLMcells(~isnan(indexUnitNamesIntoGLMcells(:,1)),1));
% and A2a
tagsInExcludedFormat=nan(size(cued_success_Response.excluded)); 
tagsInExcludedFormat(excluInds_cued_success_Response)=cued_success_Response.A2atag(cued_success_Response.excluded==0);
tagsForUnitbyunitNames=tagsInExcludedFormat(excluInds_unitbyunit_names);
A2atagForGLM=tagsForUnitbyunitNames(indexUnitNamesIntoGLMcells(~isnan(indexUnitNamesIntoGLMcells(:,1)),1));

% isCuedFromGLM=any(all_glm_coef(:,9:11)>1.5,2); %isCuedFromGLM=any(all_glm_coef(:,9:11)>0.15,2);
% isCuedFromGLM=all(all_glm_coef(:,9:11)<0.15,2); 
% isCuedFromGLM=nanmean(all_glm_coef(:,8:11),2)-nanmean(all_glm_coef(:,1:7),2)>0.3;
% isCuedFromGLM=nanmean(all_glm_coef(:,8:11),2)-nanmean(all_glm_coef(:,1:7),2)<-0.1;
% isCuedFromGLM=nanmean(all_glm_coef(:,8:40),2)-nanmean(all_glm_coef(:,1:7),2)>0.1;
isCuedFromGLM=nanmean(all_glm_coef(:,8:40),2)-nanmean(all_glm_coef(:,1:7),2)>0.1;
isCuedFromGLM_in_cued_success_Response=excluInds_unitbyunit_names(indexGLMcellsIntoUnitNames(isCuedFromGLM));

% isolate cued units
leaveOnlyCued=zeros(length(cued_success_Response.excluded)); leaveOnlyCued(isCuedFromGLM_in_cued_success_Response)=1; leaveOnlyCued=logical(leaveOnlyCued);
cued_success_Response=removeUnitFromResponse(cued_success_Response,leaveOnlyCued);
cued_drop_Response=removeUnitFromResponse(cued_drop_Response,leaveOnlyCued);
uncued_success_Response=removeUnitFromResponse(uncued_success_Response,leaveOnlyCued);
uncued_drop_Response=removeUnitFromResponse(uncued_drop_Response,leaveOnlyCued);
plotUnitSummariesAfterTCAlabels(cued_success_Response.idx,cuez,cued_success_Response,cued_drop_Response,uncued_success_Response,uncued_drop_Response,[],'justAvs','justAvs');

% use glm coeffs for binning units before tuning plots
isSusFailFromGLM=nanmean(all_glm_coef(:,150:160+40),2);
isSusFailFromGLM_inExcluInds=nan(size(cued_success_Response.excluded)); 
isSusFailFromGLM_inExcluInds(excluInds_unitbyunit_names(indexGLMcellsIntoUnitNames))=isSusFailFromGLM;
isSusFailForCurrResponseUnits=isSusFailFromGLM_inExcluInds(cued_success_Response.excluded==0);
plotUnitSummariesAfterTCAlabels(cued_success_Response.idx,isSusFailForCurrResponseUnits,cued_success_Response,cued_drop_Response,uncued_success_Response,uncued_drop_Response,[],'cued','tuning');

isOutcomeDiffer=nanmean(all_glm_coef(:,150:160),2)-(nanmean(all_glm_coef(:,189:199),2)+nanmean(all_glm_coef(:,230:240),2))./2;
[ts,allco]=plotGLMcoef(all_glm_coef(idxForGLM==1,:),[],fnames,10*0.01,9,'mean'); title('mat glm');

%% TCA for just success vs failure
% response_to_plot='all_success'; plotUnitCriteria=[1 0 0 1 0]; getCriteriaForUnitsToPlot(plotUnitCriteria);
% dd_more=cell(1,length(dd)); 
% for i=1:length(dd)
%     dd_more{i}=[dd{i} sep response_to_plot];
% end
% whichUnitsToGrab='_'; success_Response=getAndSaveResponse(dd_more,whichUnitsToGrab,settingsForStriatumUnitPlots,[]); save('Z:\MICROSCOPE\Kim\20230205 all SU alignments\all trials averaged not downsampled\success_Response.mat','success_Response');
% response_to_plot='all_failure'; plotUnitCriteria=[1 0 0 1 0]; getCriteriaForUnitsToPlot(plotUnitCriteria);
% dd_more=cell(1,length(dd)); 
% for i=1:length(dd)
%     dd_more{i}=[dd{i} sep response_to_plot];
% end
% whichUnitsToGrab='_'; failure_Response=getAndSaveResponse(dd_more,whichUnitsToGrab,settingsForStriatumUnitPlots,[]); save('Z:\MICROSCOPE\Kim\20230205 all SU alignments\all trials averaged not downsampled\failure_Response.mat','failure_Response');
load('Z:\MICROSCOPE\Kim\20230205 all SU alignments\all trials averaged not downsampled\success_Response.mat');
load('Z:\MICROSCOPE\Kim\20230205 all SU alignments\all trials averaged not downsampled\failure_Response.mat');
trial_n_cutoff=0;
out=plotVariousSUResponsesAlignedToBeh('matchUnitsAcrossResponses',excludeTooFewTrials(success_Response,trial_n_cutoff,false),excludeTooFewTrials(failure_Response,trial_n_cutoff,false),[],[],[]);
success_Response=out.Response1; failure_Response=out.Response2; 
takePointsBeforeZero=30; %15;
takePointsAfterZero=450;
dataMatrix=setUpDataMatrix(success_Response,success_Response,failure_Response,success_Response,success_Response,takePointsBeforeZero,takePointsAfterZero);
dataMatrix(dataMatrix<0)=0; dataMatrix=dataMatrix(:,31:end,:); dataMatrix(:,1:end-73,3)=dataMatrix(:,74:end,3); dataMatrix(:,1:end-73,5)=dataMatrix(:,74:end,5);
clear newDataMatrix
for i=1:size(dataMatrix,3)
    temp=reshape(dataMatrix(:,:,i),size(dataMatrix(:,:,i),1),size(dataMatrix(:,:,i),2));
    for j=1:size(temp,1)
        temp(j,:)=smoothdata(temp(j,:),'gaussian',42);%10);
    end
    newDataMatrix(:,:,i)=temp; 
end
[groupLabelsFromTCA,cuez]=principaledCA(newDataMatrix,{'units','time','conditions'},6,1,'low');

%% D-prime within unit across trials over time
analyzeProbabilityOfOnAfterOutcome(dd,[],[],[],'cued_failure','uncued_failure','overTime');
% analyzeProbabilityOfOnAfterOutcome(dd,[0 2],[],[]);

%% Attempt trial by trial classification, using labels from training set
% attemptTrialByTrialClassification(dd,[],[],'cued_success','cued_failure',[0 5],[]);
% RAN THIS NEXT LINE TO GET 'Z:\MICROSCOPE\Kim\Final Figs\Fig5\Main
% figure\attempt trial by trial
% classification\cuedsuccess_vs_cuedfailure.mat' and
% 'uncuedsuccess_vs_uncuedfailure.mat'
% fortbytclass=attemptTrialByTrialClassification(dd,[],[],'cued_success','cued_failure',[2 5],[],[],[],[]);
% still need to divide fortbytclass firing rates by number of bins in time
% range, e.g., [2 5], because attemptTrialByTrialClassification.m takes the
% sum, not average
nbins=(5-2)/0.01;
% figure out how these trial by trial units map onto
% cued_success_Response.consensus_idx
% exactly 1063 unit from cued_success_Respone tbyt, so good, those map
% 1065 units from cued_failure_Response tbyt. Based on excluded from load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\all trials\discard trials where opto during cue\cued_failure_Response.mat')
% can find mapping. Units 1:1390 are same in cuedsuccess and cuedfailure. Then cuedfailure has an extra unit, then 1 unit same, then 1 extra unit
% in failure. So need to drop units 1391 and 1393 (indices into excluded) to get cuedfailure to match cuedsuccess. In indices 1:1065, 1:1391 is 
% unit indices 1:558. So drop unit 559 and 560 of cuedfailure to get it to match cuedsuccess. See adjustment below.
% Ok, now figure this out for uncuedsuccess and uncuedfailure wrt cuedsuccess units. uncued_success_Response already has 1063 units, and
% they are the same 1063 units as cued_success_Response. Uncued_failure_Response indexes up to unit 1065, like
% cued_failure_Response. Units 559 and 560 already missing from uncued_failure_Response, but need to reassign unitids, as I did for cued_failure_Response. 
% plot results
a=load('Z:\MICROSCOPE\Kim\Final Figs\Fig5\Main figure\attempt trial by trial classification\cuedsuccess_vs_cuedfailure.mat'); a.fortbytclass.unitfr_success=a.fortbytclass.unitfr_success./nbins; a.fortbytclass.unitfr_failure=a.fortbytclass.unitfr_failure./nbins;
% a.fortbytclass.unitfr_success(a.fortbytclass.unitfr_success<0.01)=0; a.fortbytclass.unitfr_failure(a.fortbytclass.unitfr_failure<0.01)=0; 
% a.fortbytclass.unitfr_success(a.fortbytclass.unitfr_success>0 & a.fortbytclass.unitfr_success<=1)=1; a.fortbytclass.unitfr_failure(a.fortbytclass.unitfr_failure>0 & a.fortbytclass.unitfr_failure<=1)=1;
% a.fortbytclass.unitfr_success(a.fortbytclass.unitfr_success>1 & a.fortbytclass.unitfr_success<=2)=2; a.fortbytclass.unitfr_failure(a.fortbytclass.unitfr_failure>1 & a.fortbytclass.unitfr_failure<=2)=2;
% a.fortbytclass.unitfr_success(a.fortbytclass.unitfr_success>2 & a.fortbytclass.unitfr_success<=3)=3; a.fortbytclass.unitfr_failure(a.fortbytclass.unitfr_failure>2 & a.fortbytclass.unitfr_failure<=3)=3;
% a.fortbytclass.unitfr_success(a.fortbytclass.unitfr_success>3 & a.fortbytclass.unitfr_success<=4)=4; a.fortbytclass.unitfr_failure(a.fortbytclass.unitfr_failure>3 & a.fortbytclass.unitfr_failure<=4)=4;
% a.fortbytclass.unitfr_success(a.fortbytclass.unitfr_success>4 & a.fortbytclass.unitfr_success<=5)=5; a.fortbytclass.unitfr_failure(a.fortbytclass.unitfr_failure>4 & a.fortbytclass.unitfr_failure<=5)=5;
% a.fortbytclass.unitfr_success(a.fortbytclass.unitfr_success>5 & a.fortbytclass.unitfr_success<=6)=6; a.fortbytclass.unitfr_failure(a.fortbytclass.unitfr_failure>5 & a.fortbytclass.unitfr_failure<=6)=6;
% a.fortbytclass.unitfr_success(a.fortbytclass.unitfr_success>6 & a.fortbytclass.unitfr_success<=7)=7; a.fortbytclass.unitfr_failure(a.fortbytclass.unitfr_failure>6 & a.fortbytclass.unitfr_failure<=7)=7;
% a.fortbytclass.unitfr_success(a.fortbytclass.unitfr_success>7 & a.fortbytclass.unitfr_success<=8)=8; a.fortbytclass.unitfr_failure(a.fortbytclass.unitfr_failure>7 & a.fortbytclass.unitfr_failure<=8)=8;
% a.fortbytclass.unitfr_success(a.fortbytclass.unitfr_success>3)=3; a.fortbytclass.unitfr_failure(a.fortbytclass.unitfr_failure>3)=3;
figure(); plot(unique(a.fortbytclass.fromWhichUnit_success)); title('cuedsuccess units'); disp('max of cuedsuccessunits'); disp(nanmax(a.fortbytclass.fromWhichUnit_success));
dropinds=ismember(a.fortbytclass.fromWhichUnit_failure,559:560); keepinds=~ismember(1:length(a.fortbytclass.fromWhichUnit_failure),dropinds);
a.fortbytclass.unitfr_failure=a.fortbytclass.unitfr_failure(keepinds);
a.fortbytclass.fromWhichUnit_failure=a.fortbytclass.fromWhichUnit_failure(keepinds);
a.fortbytclass.fromWhichTrial_failure=a.fortbytclass.fromWhichTrial_failure(keepinds);
a.fortbytclass.fromWhichSess_failure=a.fortbytclass.fromWhichSess_failure(keepinds);
a.fortbytclass.fromWhichTrialID_failure=a.fortbytclass.fromWhichTrialID_failure(keepinds);
% reassign unit ids
currunitids=[1:558 561:1065];
newunitids=1:1063;
for i=1:length(currunitids)
    a.fortbytclass.fromWhichUnit_failure(ismember(a.fortbytclass.fromWhichUnit_failure,currunitids(i)))=newunitids(i);
end
figure(); plot(unique(a.fortbytclass.fromWhichUnit_failure)); title('cuedfailure units'); disp('max of cuedfailureunits'); disp(nanmax(a.fortbytclass.fromWhichUnit_failure));
b=load('Z:\MICROSCOPE\Kim\Final Figs\Fig5\Main figure\attempt trial by trial classification\uncuedsuccess_vs_uncuedfailure.mat'); b.fortbytclass.unitfr_success=b.fortbytclass.unitfr_success./nbins; b.fortbytclass.unitfr_failure=b.fortbytclass.unitfr_failure./nbins;
% b.fortbytclass.unitfr_success(b.fortbytclass.unitfr_success<0.01)=0; b.fortbytclass.unitfr_failure(b.fortbytclass.unitfr_failure<0.01)=0; 
% b.fortbytclass.unitfr_success(b.fortbytclass.unitfr_success>0 & b.fortbytclass.unitfr_success<=1)=1; b.fortbytclass.unitfr_failure(b.fortbytclass.unitfr_failure>0 & b.fortbytclass.unitfr_failure<=1)=1;
% b.fortbytclass.unitfr_success(b.fortbytclass.unitfr_success>1 & b.fortbytclass.unitfr_success<=2)=2; b.fortbytclass.unitfr_failure(b.fortbytclass.unitfr_failure>1 & b.fortbytclass.unitfr_failure<=2)=2;
% b.fortbytclass.unitfr_success(b.fortbytclass.unitfr_success>2 & b.fortbytclass.unitfr_success<=3)=3; b.fortbytclass.unitfr_failure(b.fortbytclass.unitfr_failure>2 & b.fortbytclass.unitfr_failure<=3)=3;
% b.fortbytclass.unitfr_success(b.fortbytclass.unitfr_success>3 & b.fortbytclass.unitfr_success<=4)=4; b.fortbytclass.unitfr_failure(b.fortbytclass.unitfr_failure>3 & b.fortbytclass.unitfr_failure<=4)=4;
% b.fortbytclass.unitfr_success(b.fortbytclass.unitfr_success>4 & b.fortbytclass.unitfr_success<=5)=5; b.fortbytclass.unitfr_failure(b.fortbytclass.unitfr_failure>4 & b.fortbytclass.unitfr_failure<=5)=5;
% b.fortbytclass.unitfr_success(b.fortbytclass.unitfr_success>5 & b.fortbytclass.unitfr_success<=6)=6; b.fortbytclass.unitfr_failure(b.fortbytclass.unitfr_failure>5 & b.fortbytclass.unitfr_failure<=6)=6;
% b.fortbytclass.unitfr_success(b.fortbytclass.unitfr_success>6 & b.fortbytclass.unitfr_success<=7)=7; b.fortbytclass.unitfr_failure(b.fortbytclass.unitfr_failure>6 & b.fortbytclass.unitfr_failure<=7)=7;
% b.fortbytclass.unitfr_success(b.fortbytclass.unitfr_success>7 & b.fortbytclass.unitfr_success<=8)=8; b.fortbytclass.unitfr_failure(b.fortbytclass.unitfr_failure>7 & b.fortbytclass.unitfr_failure<=8)=8;
% b.fortbytclass.unitfr_success(b.fortbytclass.unitfr_success>3)=3; b.fortbytclass.unitfr_failure(b.fortbytclass.unitfr_failure>3)=3;
disp('max of uncuedsuccessunits'); disp(nanmax(b.fortbytclass.fromWhichUnit_success));
% reassign unit ids
currunitids=[1:558 561:1065];
newunitids=1:1063;
for i=1:length(currunitids)
    b.fortbytclass.fromWhichUnit_failure(ismember(b.fortbytclass.fromWhichUnit_failure,currunitids(i)))=newunitids(i);
end
disp('max of uncuedfailureunits'); disp(nanmax(b.fortbytclass.fromWhichUnit_failure));
decodeTrialByTrialType(a.fortbytclass,b.fortbytclass,cued_success_Response.consensus_idx,100,100,80,false,false,false,false,false); % nBoots,nUnits,nTrials,withReplacement,addThirdAxis,nanAllZeros,justBoostrapTrials,collapseWithinUnit

% gp1 bins: binsForTuning{1}=[-10 -0.0001 10]; % greater than bottom of bin, less than or equal to top of bin
% gp2 bins: binsForTuning{1}=[-10 0 10];
load('Z:\MICROSCOPE\Kim\Final Figs\Fig5\Main figure\cued_success_Response_w_py_metrics.mat'); backup_consensus_idx=cued_success_Response.consensus_idx;
cidx=backup_consensus_idx; cidx(cued_success_Response.cXfail_sus1to5sec-cued_success_Response.allfail_sus1to5sec<=-0.0001 & backup_consensus_idx==1)=nan; % throw out bin1 for gp1
cidx(cued_success_Response.cXsucc_sus1to5sec<=0 & backup_consensus_idx==2)=nan; % throw out bin1 for gp2
out_cuedir2=decodeTrialByTrialType(a.fortbytclass,b.fortbytclass,cidx,100,100,100,false,false,false,true,false); close all;
cidx=backup_consensus_idx; cidx(cued_success_Response.cXfail_sus1to5sec-cued_success_Response.allfail_sus1to5sec>-0.0001 & backup_consensus_idx==1)=nan; % throw out bin2 for gp1
cidx(cued_success_Response.cXsucc_sus1to5sec>0 & backup_consensus_idx==2)=nan; % throw out bin2 for gp2
out_cuedir1=decodeTrialByTrialType(a.fortbytclass,b.fo+rtbytclass,cidx,100,100,100,false,false,false,true,false); close all;
figure(); 
scatter(-out_cuedir2.cuedsucc_temp1+out_cuedir1.cuedsucc_temp1,out_cuedir2.cuedsucc_temp2-out_cuedir1.cuedsucc_temp2,[],'g'); hold on;
scatter(-out_cuedir2.cuedfail_temp1+out_cuedir1.cuedfail_temp1,out_cuedir2.cuedfail_temp2-out_cuedir1.cuedfail_temp2,[],'r');
scatter(-out_cuedir2.uncuedsucc_temp1+out_cuedir1.uncuedsucc_temp1,out_cuedir2.uncuedsucc_temp2-out_cuedir1.uncuedsucc_temp2,[],'b');
scatter(-out_cuedir2.uncuedfail_temp1+out_cuedir1.uncuedfail_temp1,out_cuedir2.uncuedfail_temp2-out_cuedir1.uncuedfail_temp2,[],'y');
xlabel('cue vs uncue gp1'); ylabel('cue vs uncue gp2');
% neuron type shuffle
% decodeTrialByTrialType(a.fortbytclass,b.fortbytclass,cued_success_Response.consensus_idx(randperm(length(cued_success_Response.consensus_idx))),100,100,70,false,false,false);

%% Behavior controls
% load some behavior data
plotPhotometryResult(beh2_tbt,beh2_tbt,[],'all_reachBatch','isChewing','cueZone_onVoff','first',[0 3],[]);

%% colormaps
% SUCCESS V FAILURE 
% or consider figure(); imagesc(1:100); colormap(brewermap(100,"PiYG")); to
% deconflict with opto
% [figure(); imagesc(1:100); colormap(brewermap(10,"RdYlGn"));] 
% OPTO red 'r'
% CUE blue [50 136 189]./255
% learning 'cool'
% CUE TUNING
% figure(); imagesc(1:100); colormap(brewermap(100,"Purples"));
% TRIAL IS CUED OR UNCUED
% cued jet(1), uncued jet(100,:)