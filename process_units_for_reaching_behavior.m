% function process_units_for_reaching_behavior()
% for processing unit alignments to behavior

clear all

onVPN=true; % if want to skip reading in spikes from raw data
goodUnitLabel=2; 
% for discarding units dead or moved away
dsinds=225;
percentThresh=5;
timeStretchThresh=60*10; % in seconds
plotInference=true;

%% 1. Get locations of data: aligned phys events and spikes, raw phys data, probe location and depth, etc.
% Cell array mapping
% Date, mouse, raw data physiology file name and directory, probe depth,
% recording location, location of trial by trial data, spikes location,
% location of SU aligned to behavior, raw data binary name
data_loc_array=cell(2,6);

i=1; data_loc_array{i,1}='20210311'; data_loc_array{i,2}='Oct_B'; 
data_loc_array{i,3}=[]; %Z:\MICROSCOPE\Kim\WHISPER recs\Oct_B\20210311\phys\OctB2__g0\ % raw data phys location, or empty
data_loc_array{i,4}=2300; % ventral depth in microns relative to bregma
data_loc_array{i,5}='pDMS-tail';
data_loc_array{i,6}='Z:\MICROSCOPE\Kim\WHISPER recs\Oct_B\20210311\tbt';
data_loc_array{i,7}='Z:\MICROSCOPE\Kim\WHISPER recs\Oct_B\20210311\spike output';
data_loc_array{i,8}='Z:\MICROSCOPE\Kim\WHISPER recs\Oct_B\20210311\SU aligned to behavior';
data_loc_array{i,9}='OctB2__g0_t0.nidq.bin';

i=2; data_loc_array{i,1}='20210311'; data_loc_array{i,2}='dLight1'; 
data_loc_array{i,3}=[]; %Z:\MICROSCOPE\Kim\WHISPER recs\dLight1\20210311\phys\dLight1__g0 % raw data phys location, or empty
data_loc_array{i,4}=2300; % ventral depth in microns relative to bregma
data_loc_array{i,5}='pDMS-tail';
data_loc_array{i,6}='/Users/kim/Desktop/Example data/20210302 dLight1/tbt';
data_loc_array{i,7}='/Users/kim/Desktop/Example data/20210302 dLight1/spike output dLight_3__g0';
data_loc_array{i,8}='/Users/kim/Desktop/Example data/20210302 dLight1/SU aligned to behavior';
data_loc_array{i,9}='dLight1.bin';

%% 2. Auto-populate SU QC and behavior alignments
% Discarding trials where unit dead or moved away
for i=1:size(data_loc_array,1)
    % load tbt (trial by trial) data
    [physiology_tbt,beh2_tbt,behavior_tbt,photometry_tbt]=loadReachingExptPhysData(data_loc_array(i,:));
    % get spikes
    dd=dir(data_loc_array{i,7});
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
            % make opto-tagged alignment
            [~,tbtspikes]=organizeSpikesToMatch_physiology_tbt(spikes,physiology_tbt);
            optoAligned_phys_tbt=alignToOpto(addGoodUnitsAsFields(physiology_tbt,tbtspikes,2,1,false,true));
            optoAligned_phys_tbt=checkWaveformsDuringOpto(optoAligned_phys_tbt,tbtspikes);
            % save opto alignment
            if ~exist([data_loc_array{i,8} sep 'opto_aligned'],'dir')
                mkdir([data_loc_array{i,8} sep 'opto_aligned']);
            end
            save([data_loc_array{i,8} sep 'opto_aligned' sep 'phys_tbt_for_' dd(j).name],'optoAligned_phys_tbt');
            clear tbtspikes
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
                trodeChsForSpikes=trodeChsForSpikes(si);
                forplot_trodeChs=trodeChsForSpikes;
                if length(forplot_trodeChs)<4
                    forplot_trodeChs=[forplot_trodeChs ones(1,4-length(forplot_trodeChs))*forplot_trodeChs(1)];
                end
                % check whether QC figure already exists for this unit
                [sue,qc_fname]=SU_QC_file_exists(data_loc_array{i,8}, currAssign, trodeChsForSpikes(end));
                if ~sue && onVPN==false
                    % make behavior alignments and single unit QC figs
                    unit_data=saveBehaviorAlignmentsSingleNeuron(physiology_tbt,spikes,currAssign,beh2_tbt,forplot_trodeChs,trodeChsForSpikes(end),data_loc_array{i,8},'',data_loc_array{i,3},data_loc_array{i,9});
                else
                    if sue==false
                        % make QC fig without reading in raw data
                        unit_data=saveBehaviorAlignmentsSingleNeuron(physiology_tbt,spikes,currAssign,beh2_tbt,forplot_trodeChs,trodeChsForSpikes(end),data_loc_array{i,8},'',data_loc_array{i,3},[]);
                    else
                        % just redo behavior alignments, skipping QC fig
                        spikes.skipQC=true;
                        saveBehaviorAlignmentsSingleNeuron(physiology_tbt,spikes,currAssign,beh2_tbt,forplot_trodeChs,trodeChsForSpikes(end),data_loc_array{i,8},'',data_loc_array{i,3},[]);
                        unit_data=[data_loc_array{i,8} sep qc_fname];
                    end
                end
                % discard trials where unit dead or moved away
                [physiology_tbt.(['unit' num2str(currAssign) '_dontUseTrials']),physiology_tbt.(['unit' num2str(currAssign) '_meanFR'])]=inferUnitStability(unit_data,physiology_tbt,dsinds,percentThresh,timeStretchThresh,plotInference);
                % get depth of unit

                % get waveform features
                [halfwidth, peakToTrough, amp, avWaveforms]=getWaveformFeatures(filtspikes_without_sweeps(spikes,0,'assigns',currAssign),spikes.params.Fs);
                physiology_tbt.(['unit' num2str(currAssign) '_dontUseTrials'])
                % save unit details for these spikes
                if ~exist([data_loc_array{i,8} sep 'unit_details'],'dir')
                    mkdir([data_loc_array{i,8} sep 'unit_details']);
                end
                save([data_loc_array{i,8} sep 'unit_details' sep 'phys_tbt_for_' dd(j).name],'physiology_tbt');
            end
        end
    end
end

%% 3. Sort units by depth, discard units not in pDMS-tail

%% 4. Sort units by type

%% 5. Identify opto-tagged units



