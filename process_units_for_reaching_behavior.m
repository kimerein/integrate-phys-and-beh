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

i=1; data_loc_array{i,1}='20210621'; data_loc_array{i,2}='Mar_6'; 
data_loc_array{i,3}=[]; % raw data phys location, or empty
data_loc_array{i,4}=2300; % ventral depth in microns relative to bregma
data_loc_array{i,5}='pDMS-tail';
data_loc_array{i,6}='/Users/kim/Desktop/Example data/20210621 Mar_6/tbt';
data_loc_array{i,7}='/Users/kim/Desktop/Example data/20210621 Mar_6/spike output Mar_6__g1';
data_loc_array{i,8}='/Users/kim/Desktop/Example data/20210621 Mar_6/SU aligned to behavior';
data_loc_array{i,9}='Mar_6.bin';

i=2; data_loc_array{i,1}='20210302'; data_loc_array{i,2}='dLight1'; 
data_loc_array{i,3}=[]; % raw data phys location, or empty
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
                otherwise 
                    error('do not recognize name of spikes mat file');
            end
            % find good units
            gu=find(spikes.labels(:,2)==goodUnitLabel);
            % unit by unit
            for k=1:length(gu)
                currU=gu(k);
                currAssign=spikes.labels(currU,1);
                % find ch where waveform biggest
                amp=nan(1,size(spikes.waveforms,3));
                for l=1:size(spikes.waveforms,3)
                    amp(l)=abs(min(reshape(mean(spikes.waveforms(spikes.assigns==3,:,l),1,'omitnan'),1,size(spikes.waveforms,2)),[],2,'omitnan'));
                end
                [~,si]=sort(amp);
                % sort trode chs for units
                trodeChsForSpikes=trodeChsForSpikes(si);
                % check whether QC figure already exists for this unit
                [sue,qc_fname]=SU_QC_file_exists(data_loc_array{i,8}, currAssign, trodeChsForSpikes(end));
                if ~sue && onVPN==false
                    % make behavior alignments and single unit QC figs
                    forplot_trodeChs=trodeChsForSpikes;
                    if length(forplot_trodeChs)<4
                        forplot_trodeChs=[forplot_trodeChs ones(1,4-length(forplot_trodeChs))*forplot_trodeChs(1)];
                    end
                    unit_data=saveBehaviorAlignmentsSingleNeuron(physiology_tbt,spikes,currAssign,beh2_tbt,forplot_trodeChs,trodeChsForSpikes(end),data_loc_array{i,8},'',data_loc_array{i,3},data_loc_array{i,9});
                else
                    if sue==false
                        % make QC fig without reading in raw data
                        unit_data=saveBehaviorAlignmentsSingleNeuron(physiology_tbt,spikes,currAssign,beh2_tbt,trodeChsForSpikes,trodeChsForSpikes(end),data_loc_array{i,8},'',data_loc_array{i,3},[]);
                    else
                        % just redo behavior alignments, skipping QC fig
                        spikes.skipQC=true;
                        saveBehaviorAlignmentsSingleNeuron(physiology_tbt,spikes,currAssign,beh2_tbt,trodeChsForSpikes,trodeChsForSpikes(end),data_loc_array{i,8},'',data_loc_array{i,3},[]);
                        unit_data=[data_loc_array{i,8} sep qc_fname];
                    end
                end
                % discard trials where unit dead or moved away
                dontUseTrials=inferUnitStability(unit_data,physiology_tbt,dsinds,percentThresh,timeStretchThresh,plotInference);
                
                
                % make opto-tagged alignment
                
                
            end
        end
    end
end

%% 3. Sort units by depth, discard units not in pDMS-tail

%% 4. Sort units by type

%% 5. Identify opto-tagged units

