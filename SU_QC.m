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
% Column 2, Row 1 - autocorrelation
% Column 2, Row 2 - waveform stability and firing rate
% Column 3, Row 1 - LFP segment
% Column 3, Row 2 - HF segment
% Column 4, Row 1 - LFP closeup
% Column 4, Row 2 - HF closeup
% Next columns - PCA space
% Next columns - PSTH aligned to behavior events
n_rows=2;
PC_mode='unique'; % either 'unique' or 'all'
switch PC_mode
    case 'all'
        show_first_N_PCs=10; % number of PCs to show
        show_N_views_of_PCA_space=factorial(show_first_N_PCs)/(factorial(2)*factorial(show_first_N_PCs-2));
    case 'unique'
        show_first_N_PCs=4; % will show combos of this many PCs (top PCs that separate this unit from rest of spikes)
        show_N_views_of_PCA_space=factorial(show_first_N_PCs)/(factorial(2)*factorial(show_first_N_PCs-2));
end
N_behavior_alignments=4;
PCA_columns=ceil(show_N_views_of_PCA_space/n_rows);
beh_columns=ceil(N_behavior_alignments/n_rows);
n_columns=3+PCA_columns+beh_columns;

% Format figure
[mainfig, ha]=formatFigure(n_rows, n_columns);

% # refractory period violations (RPVs)
rpv=count_RPV(spikes, unit_assign);
total_spikes=countTotalSpikes(spikes, unit_assign);

% Unit waveform
wvfms=unitWaveform(spikes, unit_assign);
ind_into_ha=1;
ind_into_ha=nextPlot(ha, ind_into_ha);
wvfmPlot(wvfms, spikes.params.Fs, rpv, total_spikes);

% Spike amplitude histogram
ind_into_ha=nextPlot(ha, ind_into_ha);
plot_detection_criterion(spikes, unit_assign);

% Autocorrelation
ind_into_ha=nextPlot(ha, ind_into_ha);
plot_autocorr(spikes, unit_assign);

% Waveform stability and firing rate
ind_into_ha=nextPlot(ha, ind_into_ha);
plotWaveformAmplitudeStability(spikes, unit_assign); hold on;
plotFiringRate(spikes, unit_assign, 100); % last argument is binsize in ms for histogram of firing rate

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
switch PC_mode
    case 'all'
        for i=1:show_first_N_PCs
            for j=i+1:show_first_N_PCs
                ind_into_ha=nextPlot(ha, ind_into_ha);
                plotPCspace(spikes, unit_assign, i, j);
            end
        end
    case 'unique'
        distance=findUniquePCs(spikes, unit_assign);
        [~,si]=sort(distance,'descend');
        for i=1:show_first_N_PCs
            for j=i+1:show_first_N_PCs
                ind_into_ha=nextPlot(ha, ind_into_ha);
                plotPCspace(spikes, unit_assign, si(i), si(j));
            end
        end
end

% Behavior-triggered PSTHs

end

function plotFiringRate(spikes, unit_assign, binsize)

% binsize in ms
spiketimes=spikes.spiketimes(ismember(spikes.assigns, unit_assign));
[n,x]=histcounts(spiketimes,0:binsize/1000:max(spiketimes));
[n,x]=cityscape_hist(n,x);
yyaxis right
plot(x,n);
xlabel('Time (s)');
ylabel('Firing rate (1/s)');

end

function [new_n,new_x]=cityscape_hist(n,x)

new_x=nan(1,length(x)*2-1);
j=1;
x_step=mode(diff(x));
new_n=nan(1,length(x)*2-1);
for i=1:length(x)-1
    new_x(j)=x(i);
    new_n(j)=n(i);
    j=j+1;
    new_x(j)=x(i)+x_step;
    new_n(j)=n(i);
    j=j+1;
end

end

function plotWaveformAmplitudeStability(spikes, unit_assign)

which=ismember(spikes.assigns, unit_assign);
onChannel=mode(spikes.info.detect.event_channel(ismember(spikes.assigns, unit_assign)));
y=range(squeeze(spikes.waveforms(which,:,onChannel)),2);
x=spikes.spiketimes(ismember(spikes.assigns, unit_assign));
scatter(x,y,[],'b','filled','MarkerFaceAlpha',0.2);
xlabel('Time (s)');
ylabel('Waveform amp','Color','b');

end

function n_spikes=countTotalSpikes(spikes, unit_assign)

n_spikes=sum(ismember(spikes.assigns, unit_assign));

end

function distance=findUniquePCs(spikes, unit_assign)

% subsample other spikes
which_other=randsample(find(~ismember(spikes.assigns, unit_assign)),sum(ismember(spikes.assigns, unit_assign)));
% spikes in unit
which=ismember(spikes.assigns, unit_assign);

distance=nan(1,size(spikes.info.pca.v,2));
for i=1:size(spikes.info.pca.v,2)
    pc_other=spikes.waveforms(which_other,:) * spikes.info.pca.v(:,i);
    pc=spikes.waveforms(which,:) * spikes.info.pca.v(:,i);
    distance(i)=abs(mean(pc_other)-mean(pc))/sqrt(std(pc_other).^2 + std(pc).^2);
end

end

function plotPCspace(spikes, unit_assign, pc_dim1, pc_dim2)

% subsample other spikes
which=randsample(find(~ismember(spikes.assigns, unit_assign)),sum(ismember(spikes.assigns, unit_assign)));
plotUnitInPCspace(spikes, which, pc_dim1, pc_dim2, 'c'); hold on;

% plot spikes in unit
which=ismember(spikes.assigns, unit_assign);
plotUnitInPCspace(spikes, which, pc_dim1, pc_dim2, 'y');
xlabel(['PC' num2str(pc_dim1)]);
ylabel(['PC' num2str(pc_dim2)]);

end

function plotUnitInPCspace(spikes, which, pc_dim1, pc_dim2, c)

x = spikes.waveforms(which,:) * spikes.info.pca.v(:,pc_dim1);
y = spikes.waveforms(which,:) * spikes.info.pca.v(:,pc_dim2);
scatter(x,y,[],c,'filled','MarkerFaceAlpha',0.2);

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

function wvfmPlot(wvfms, Fs, rpv, total_spikes)

wvfms=wvfms(:,[2 1 3 4]);
wvfms=[wvfms; nan(1,size(wvfms,2))];
wvfms=wvfms(1:end);
plot(0:1/Fs:(size(wvfms,2)-1)*(1/Fs),wvfms,'Color','k','LineWidth',1);
xlabel('Time (s)');
ylabel([num2str(total_spikes) ' spikes (# RPVs = ' num2str(rpv) ')']);

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

function plot_autocorr(spikes, unit_assign)     

spiketimes = spikes.spiketimes(ismember(spikes.assigns, unit_assign));
shadow = spikes.params.shadow;
rp     = spikes.params.refractory_period;

maxlag = 0.1; % in seconds
corr_bin_size=2; % in ms

% calculate autocorrelation
if length(spiketimes) > 1
    [cross,lags] = pxcorr(spiketimes,spiketimes, round(1000/corr_bin_size), maxlag);
else
    cross = 0;  lags = 0;
end
cross(lags==0) = 0;

% plot autocorrelation histogram
bb = bar(lags*1000,cross,1.0); hold on;
set(bb,'FaceColor',[0 0 0 ],'EdgeColor',[0 0 0 ])
% place patches to represent shadow and refractory period
ymax = max(cross) + 1;
patch(shadow*[-1 1 1 -1], [0 0 ymax ymax], [.5 .5 .5],'EdgeColor', 'none'); hold on;
patch([shadow [rp rp] shadow ], [0 0 ymax ymax], [1 0 0],'EdgeColor','none');
patch(-[shadow [rp rp] shadow ], [0 0 ymax ymax], [1 0 0],'EdgeColor','none');

% set axes
set(gca, 'XLim', maxlag*1000*[-1 1]);
set(gca,'YLim',[0 ymax])
xlabel('Time lag (msec)')
ylabel('Autocorrelation (Hz)')
    
end