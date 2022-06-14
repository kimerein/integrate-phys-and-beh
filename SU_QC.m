function SU_QC(spikes, unit_assign, unit_on_channel, raw_data_filename, raw_data_directory, behavior_timepoints, varAdditionalInputs)

varAdditionalInputs.furtherProcessData=@furtherProcessWHISPER;

% Creates summary figure for quality control (QC) of single units (SU)

% spikes structure
% spikes.params.Fs          sampling rate of signal
% spikes.waveforms          for each spike, waveform around spike minimum
%                           (peak of action potential), as a 3D array
%                           N X M X P, where N is total number of spikes,
%                           M is number of samples over time, and P is
%                           number of channels used to detect the spike 
%                           (e.g., a tetrode would have P = 4)
% spikes.params.refractory_period 
%                           refractory period in milliseconds
%
% varAdditionalInputs       structure containing any other necessary inputs
% varAdditionalInputs may contain
%       varAdditionalInputs.timeWindowToPlotAtSpike 
%                           time window around a spike for plotting
%                           high-passed and low-passed raw data

% Layout 
% Column 1, Row 1 - average spike waveform
% Column 1, Row 2 - spike amplitude distribution
% Column 2, Row 1 - LFP segment
% Column 2, Row 2 - HF segment
% Column 3, Row 1 - LFP closeup
% Column 3, Row 2 - HF closeup
% Next columns - PCA space
% Next columns - PSTH aligned to behavior events
n_rows=2;
show_N_views_of_PCA_space=4;
N_behavior_alignments=4;
PCA_columns=ceil(show_N_views_of_PCA_space/n_rows);
beh_columns=ceil(N_behavior_alignments/n_rows);
n_columns=3+PCA_columns+beh_columns;

% Format figure
[mainfig, ha]=formatFigure(n_rows, n_columns);

% # refractory period violations (RPVs)
rpv=count_RPV(spikes, unit_assign);

% Unit waveform
wvfms=unitWaveform(spikes, unit_assign);
ind_into_ha=1;
ind_into_ha=nextPlot(ha, ind_into_ha);
wvfmPlot(wvfms, spikes.params.Fs, rpv);

% Spike amplitude histogram
ind_into_ha=nextPlot(ha, ind_into_ha);
plot_detection_criterion(spikes, unit_assign);

% Get some raw data
s=getUnitSpiketimes(spikes, unit_assign);
varAdditionalInputs.firstNSpikes=20;
varAdditionalInputs.firstNmins=s(varAdditionalInputs.firstNSpikes)/60; % for WHISPER only
[HFValues,LFPValues,ADFreq]=rawDataFromChannel(raw_data_filename, raw_data_directory, unit_on_channel, varAdditionalInputs);
times=timepoints(ADFreq,HFValues);

% LFP
ind_into_ha=nextPlot(ha, ind_into_ha);
plotLFP(LFPValues, times, varAdditionalInputs, s);

% High pass-filtered
ind_into_ha=nextPlot(ha, ind_into_ha);
plotHF(HFValues, times, varAdditionalInputs, s);

% PCA space

% Behavior-triggered PSTHs

end

function times=timepoints(ADFreq, HFValues)

timestep=1/ADFreq;
times=0:timestep:(length(HFValues)-1)*timestep;

end

function plotHF(HFValues, times, varAdditionalInputs, spiketimes)

if isfield(varAdditionalInputs, 'timeWindowToPlotAtSpike')
    window=varAdditionalInputs.timeWindowToPlotAtSpike;
else
    window=5; % default time window in secs
end
if isfield(varAdditionalInputs, 'firstNSpikes')
    plotSpike=ceil(varAdditionalInputs.firstNSpikes/2);
else
    plotSpike=1; % default spike to plot
end
if spiketimes(plotSpike)-(window/2)<times(1)
    startAtTime=times(1);
else
    startAtTime=spiketimes(plotSpike)-(window/2);
end
if spiketimes(plotSpike)+(window/2)>times(end)
    endAtTime=times(end);
else
    endAtTime=spiketimes(plotSpike)+(window/2);
end
[~,startind]=min(abs(times-startAtTime),[],'all','omitnan');
[~,endind]=min(abs(times-endAtTime),[],'all','omitnan');

plot(times(startind:endind),HFValues(startind:endind),'Color','k');

% Label spike times
s=spiketimes(spiketimes>startAtTime & spiketimes<endAtTime);
hold on;
scatter(s,max(HFValues(startind:endind),[],'all','omitnan'),50,'r','Marker','*');
axis tight
xlabel('Time (s)');
ylabel('High-passed');

end

function plotLFP(LFPValues, times, varAdditionalInputs, spiketimes)

if isfield(varAdditionalInputs, 'timeWindowToPlotAtSpike')
    window=varAdditionalInputs.timeWindowToPlotAtSpike;
else
    window=5; % default time window in secs
end
if isfield(varAdditionalInputs, 'firstNSpikes')
    plotSpike=ceil(varAdditionalInputs.firstNSpikes/2);
else
    plotSpike=1; % default spike to plot
end
if spiketimes(plotSpike)-(window/2)<times(1)
    startAtTime=times(1);
else
    startAtTime=spiketimes(plotSpike)-(window/2);
end
if spiketimes(plotSpike)+(window/2)>times(end)
    endAtTime=times(end);
else
    endAtTime=spiketimes(plotSpike)+(window/2);
end
[~,startind]=min(abs(times-startAtTime),[],'all','omitnan');
[~,endind]=min(abs(times-endAtTime),[],'all','omitnan');

plot(times(startind:endind),LFPValues(startind:endind),'Color','k');

% Label spike times
s=spiketimes(spiketimes>startAtTime & spiketimes<endAtTime);
hold on;
scatter(s,max(LFPValues(startind:endind),[],'all','omitnan'),50,'r','Marker','*');
axis tight
xlabel('Time (s)');
ylabel('LFP');

end

function s=getUnitSpiketimes(spikes, unit_assign)

s=spikes.spiketimes(ismember(spikes.assigns, unit_assign));

end

function data=bandpass(data, Fs, lowCutoff, highCutoff)

% Low-pass-filter data
data=fftFilt_short(data',Fs,highCutoff,1);
data=data';
% High-pass-filter data
data=fftFilt_short(data',Fs,lowCutoff,2);
data=real(data');
% Align all traces to initial value = 0
for i=1:size(data,1)
    initVal=data(i,1);
    data(i,:)=data(i,:)-initVal;
end

end

function [HFout,LFPout,Fs]=rawDataFromChannel(raw_data_filename, raw_data_directory, unit_on_channel, varAdditionalInputs)

% Replace with your own function

% For WHISPER system
[HFout,LFPout,Fs]=readRawData(raw_data_filename, raw_data_directory, unit_on_channel, varAdditionalInputs.firstNmins, varAdditionalInputs);

end

function [Values,ADFreq]=getRawData(saveTo, savedAs)

% Formatted for WHISPER data 
% can replace with your own function
% read raw data and return
% Values = raw data values
% ADFreq = sampling rate
a=load([saveTo filesep savedAs]);
Values=a.data.Values;
ADFreq=a.data.ADFreq;

end

function [HFout,LFPout,Fs]=readRawData(raw_data_filename, raw_data_directory, readCh, firstNmins, varAdditionalInputs)

% Replace with your own code

% Kim uses this for WHISPER system
[data,Fs]=readWHISPER(raw_data_filename, raw_data_directory, varAdditionalInputs.trodeChs, firstNmins);
f=find(varAdditionalInputs.trodeChs==readCh);
LFPout=bandpass(data(f,:), Fs, 0.1, 300);
for i=1:size(data,1)
    data(i,:)=bandpass(data(i,:), Fs, 300, Fs);
end
if isfield(varAdditionalInputs,'furtherProcessData')
    if ~isempty(varAdditionalInputs.furtherProcessData) % this field should contain a function handle
        [data(1,:),data(2,:),data(3,:),data(4,:)]=varAdditionalInputs.furtherProcessData(data(1,:),data(2,:),data(3,:),data(4,:),Fs);
    end
end
HFout=data(f,:);

end

function rpv=count_RPV(spikes, unit_assign)

% counts refractory period violations
spiketimes =  sort( spikes.unwrapped_times(ismember(spikes.assigns, unit_assign)) );
rpv  = sum( diff(spiketimes)  <= (spikes.params.refractory_period * .001) );

end

function ind_into_ha=nextPlot(ha, ind_into_ha)

axes(ha(ind_into_ha));
ind_into_ha=ind_into_ha+1;

end

function wvfms=unitWaveform(spikes, unit_assign)

wvfms=nan(size(spikes.waveforms,2),size(spikes.waveforms,3));
for i=1:size(spikes.waveforms,3)
    wvfms(:,i)=reshape(mean(spikes.waveforms(spikes.assigns==unit_assign,:,i),1,'omitnan'),[size(spikes.waveforms,2) 1]);
end

end

function wvfmPlot(wvfms, Fs, rpv)

wvfms=wvfms(:,[2 1 3 4]);
wvfms=[wvfms; nan(1,size(wvfms,2))];
wvfms=wvfms(1:end);
plot(0:1/Fs:(size(wvfms,2)-1)*(1/Fs),wvfms,'Color','k','LineWidth',1);
xlabel('Time (s)');
ylabel(['# RPVs = ' num2str(rpv)]);

end

function setXlims(ha, xranges)

for i=1:length(ha)
    set(ha(i),'XLim',xranges{i});
end

end

function spawnIndividualFigures(ha, whichSubfigs)

for i=1:length(ha)
    if ~ismember(i, whichSubfigs)
        continue
    end
    % spawn individual figures
    f=figure();
    newax=copyobj(ha(i),f);
    set(newax,'Position',[0.1 0.1 0.85 0.85]);
end

end

function [mainfig, ha]=formatFigure(n_rows, n_columns)

% Set up figure layout
mainfig=figure();
Nh=n_rows; % number of rows
Nw=n_columns; % number of columns
gap=[.1 .03]; % between plots
marg_h=[.1 .01]; % margin
marg_w=[.1 .01];% marg_w=[.01 .01]; % margin
[ha,pos]=tight_subplot(Nh,Nw,gap,marg_h,marg_w);
% reorder ha so populates down rows within column first
reorder=[];
for i=0:Nw-1
    reorder=[reorder i+[1:Nw:Nw*Nh]];
end
ha=ha(reorder);

end