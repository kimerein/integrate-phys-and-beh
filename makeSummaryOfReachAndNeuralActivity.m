function makeSummaryOfReachAndNeuralActivity(physiology_tbt,phys_beh_tbt,photometry_tbt,photo_beh_tbt,spikesDir,spikes)

spikeName='spikes';
binsize=10; % in ms, for PSTHs
goodUnitLabel=2; 
bsmooth=1;
maxTrialLength=9; % in seconds
normalizeSU=true;

% figure 1
alignments={'cue','all_reachBatch'};
timewindows={[-1 16],[-1 16]};
withintimewindow={[],[]};
beh_fields={'all_reachBatch','fidgetData'};
photo_fields={'green_ch'};
phys_fields={'unit_by_unit','sum_over_singleunit'}; % each field must contain "unit"
xranges={[0 maxTrialLength],[0 maxTrialLength]};
physthenphoto_fields(1:length(phys_fields))=phys_fields;
physthenphoto_fields(length(phys_fields)+1:length(phys_fields)+length(photo_fields))=photo_fields;
isPhysField=[ones(size(1:length(phys_fields))) zeros(size(1:length(photo_fields)))];

% load spikes from directory
spike_d=dir(spikesDir);
all_wvfms=[];
all_halfWidths=[];
all_depths=[];
all_assigns=[];
unitnames={};
whichTrode=[];
firstloadedspikes=true;
if ~isempty(spikes)
    [~,spikes]=organizeSpikesToMatch_physiology_tbt(spikes,physiology_tbt);
    % get PSTH for each unit
    if firstloadedspikes==true
        [physiology_tbt,~,unit_fieldnames]=addGoodUnitsAsFields(physiology_tbt,spikes,goodUnitLabel,binsize,bsmooth,true);
        firstloadedspikes=false;
    else
        [physiology_tbt,~,unit_fieldnames]=addGoodUnitsAsFields(physiology_tbt,spikes,goodUnitLabel,binsize,bsmooth,false);
    end
    if ~isfield(physiology_tbt,'sum_over_singleunit')
        physiology_tbt.sum_over_singleunit=zeros(size(physiology_tbt.unitsum));
    end
    physiology_tbt.sum_over_singleunit=physiology_tbt.sum_over_singleunit+physiology_tbt.unitsum;
    unitnames(length(unitnames)+1:length(unitnames)+length(unit_fieldnames))=unit_fieldnames;
    % get info to classify units
    [unit_wvfms,unit_halfWidths,unit_depths,useAssigns]=getUnitWaveformAndDepth(spikes);
    all_wvfms=[all_wvfms; unit_wvfms(spikes.labels(:,2)==goodUnitLabel,:)];
    all_assigns=[all_assigns useAssigns(spikes.labels(:,2)==goodUnitLabel)];
    all_halfWidths=[all_halfWidths unit_halfWidths(spikes.labels(:,2)==goodUnitLabel)];
    all_depths=[all_depths unit_depths(spikes.labels(:,2)==goodUnitLabel)];
    whichTrode=[whichTrode ones(size(unit_depths(spikes.labels(:,2)==goodUnitLabel)))*1];
else
    for i=1:length(spike_d)
        spind=regexp(spike_d(i).name,spikeName);
        if ~isempty(spind)
            % load spikes
            [sta,en]=regexp(spike_d(i).name,'(?<name>\d+)');
            if isempty(sta)
                currtrode=1;
            else
                currtrode=str2double(spike_d(i).name(sta:en));
            end
            disp(['loading ' spike_d(i).name]);
            a=load([spikesDir '\' spike_d(i).name]);
            spikes=a.spikes;
        else
            continue
        end
        [~,spikes]=organizeSpikesToMatch_physiology_tbt(spikes,physiology_tbt);
        % get PSTH for each unit
        if firstloadedspikes==true
            [physiology_tbt,~,unit_fieldnames]=addGoodUnitsAsFields(physiology_tbt,spikes,goodUnitLabel,binsize,bsmooth,true);
            firstloadedspikes=false;
        else
            [physiology_tbt,~,unit_fieldnames]=addGoodUnitsAsFields(physiology_tbt,spikes,goodUnitLabel,binsize,bsmooth,false);
        end
        if ~isfield(physiology_tbt,'sum_over_singleunit')
            physiology_tbt.sum_over_singleunit=zeros(size(physiology_tbt.unitsum));
        end
        physiology_tbt.sum_over_singleunit=physiology_tbt.sum_over_singleunit+physiology_tbt.unitsum;
        unitnames(length(unitnames)+1:length(unitnames)+length(unit_fieldnames))=unit_fieldnames;
        % get info to classify units
        [unit_wvfms,unit_halfWidths,unit_depths,useAssigns]=getUnitWaveformAndDepth(spikes);
        all_wvfms=[all_wvfms; unit_wvfms(spikes.labels(:,2)==goodUnitLabel,:)];
        all_assigns=[all_assigns useAssigns(spikes.labels(:,2)==goodUnitLabel)];
        all_halfWidths=[all_halfWidths unit_halfWidths(spikes.labels(:,2)==goodUnitLabel)];
        all_depths=[all_depths unit_depths(spikes.labels(:,2)==goodUnitLabel)];
        whichTrode=[whichTrode ones(size(unit_depths(spikes.labels(:,2)==goodUnitLabel)))*currtrode];
    end
end

% make figures
% Set up figure layout 1
mainfig1=figure();
Nh=length(beh_fields)+length(photo_fields)*2+length(phys_fields); % number of rows
Nw=length(alignments); % number of columns
gap=[.01 .03]; % between plots
marg_h=[.1 .01]; % margin
marg_w=[.1 .01];% marg_w=[.01 .01]; % margin
[ha,pos]=tight_subplot(Nh,Nw,gap,marg_h,marg_w);
% reorder ha so populates down rows within column first
reorder=[];
for i=0:Nw-1
    reorder=[reorder i+[1:Nw:Nw*Nh]];
end
ha=ha(reorder);
indintoha=1;
for i=1:Nw
    % each alignment
    for j=1:length(physthenphoto_fields)
        % each type of data to plot
        for k=1:length(beh_fields) % get behavior fields aligned
            if j==1
                if isPhysField(j)==1
                    if k==1
                        [fout,dataout,n_events_in_av,alignmentCompanion,f_heatmap,plotBehFieldOut,su]=plotPhysiologyResult(physiology_tbt,phys_beh_tbt,[],alignments{i},physthenphoto_fields{j},beh_fields{k},withintimewindow{i},timewindows{i},{nan ha(indintoha)});
                        plottedIntoWhichAxes=indintoha;
                        if ~strcmp(physthenphoto_fields{j},'unit_by_unit')
                            strtoplot=[physthenphoto_fields{j} ' aligned to ' alignments{i}];
                            f_=regexp(strtoplot,'_');
                            strtoplot(f_)=' ';
                            text(nanmin(xranges{i}),nanmax(nanmean(dataout.y,1)),strtoplot);
                        end
                        indintoha=indintoha+1;
                        closeAllBut(mainfig1);
                    else
                        [fout,dataout,n_events_in_av,alignmentCompanion,f_heatmap,plotBehFieldOut,su]=plotPhysiologyResult(physiology_tbt,phys_beh_tbt,[],alignments{i},physthenphoto_fields{j},beh_fields{k},withintimewindow{i},timewindows{i},{nan nan});
                        closeAllBut(mainfig1);
                    end
                    axes(ha(indintoha));
                    indintoha=indintoha+1;
                    plotWStderr(plotBehFieldOut.y,plotBehFieldOut.x,'k',[],size(plotBehFieldOut.y,1));
                    hold on;
                    rescale=nanmax(nanmean(plotBehFieldOut.y,1))/nanmax(nanmean(alignmentCompanion.y,1));
                    plotWStderr(alignmentCompanion.y.*rescale,alignmentCompanion.x,'g',[],size(alignmentCompanion.y,1));
                    strtoplot=[beh_fields{k} ' aligned to ' alignments{i}];
                    f_=regexp(strtoplot,'_');
                    strtoplot(f_)=' ';
                    text(nanmin(xranges{i}),nanmax(nanmean(plotBehFieldOut.y,1)),strtoplot);
                else
                    if k==1
                        [fout,dataout,n_events_in_av,alignmentCompanion,f_heatmap,plotBehFieldOut]=plotPhotometryResult(photometry_tbt,photo_beh_tbt,[],alignments{i},physthenphoto_fields{j},beh_fields{k},withintimewindow{i},timewindows{i},{ha(indintoha) ha(indintoha+1)});
                        text(nanmin(xranges{i}),nanmax(dataout.y),[physthenphoto_fields{j} ' aligned to ' alignments{i}]);
                        indintoha=indintoha+2;
                        closeAllBut(mainfig1);
                    else
                        [fout,dataout,n_events_in_av,alignmentCompanion,f_heatmap,plotBehFieldOut]=plotPhotometryResult(photometry_tbt,photo_beh_tbt,[],alignments{i},physthenphoto_fields{j},beh_fields{k},withintimewindow{i},timewindows{i},{nan nan});
                        closeAllBut(mainfig1);
                    end
                    axes(ha(indintoha));
                    indintoha=indintoha+1;
                    plotWStderr(plotBehFieldOut.y,plotBehFieldOut.x,'k',[],size(plotBehFieldOut.y,1));
                    hold on;
                    rescale=nanmax(nanmean(plotBehFieldOut.y,1))/nanmax(nanmean(alignmentCompanion.y,1));
                    plotWStderr(alignmentCompanion.y.*rescale,alignmentCompanion.x,'g',[],size(alignmentCompanion.y,1));
                    strtoplot=[beh_fields{k} ' aligned to ' alignments{i}];
                    f_=regexp(strtoplot,'_');
                    strtoplot(f_)=' ';
                    text(nanmin(xranges{i}),nanmax(nanmean(plotBehFieldOut.y,1)),strtoplot);
                end
            else
                if isPhysField(j)==1
                    if k==1
                        [fout,dataout,n_events_in_av,alignmentCompanion,f_heatmap,plotBehFieldOut,su]=plotPhysiologyResult(physiology_tbt,phys_beh_tbt,[],alignments{i},physthenphoto_fields{j},[],withintimewindow{i},timewindows{i},{nan ha(indintoha)});
                        plottedIntoWhichAxes=indintoha;
                        if ~strcmp(physthenphoto_fields{j},'unit_by_unit')
                            strtoplot=[physthenphoto_fields{j} ' aligned to ' alignments{i}];
                            f_=regexp(strtoplot,'_');
                            strtoplot(f_)=' ';
                            text(nanmin(xranges{i}),nanmax(nanmean(dataout.y,1)),strtoplot);
                        end
                        indintoha=indintoha+1;
                        closeAllBut(mainfig1);
                    else
                        continue
                    end
                else
                    if k==1
                        [fout,dataout,n_events_in_av,alignmentCompanion,f_heatmap,plotBehFieldOut]=plotPhotometryResult(photometry_tbt,photo_beh_tbt,[],alignments{i},physthenphoto_fields{j},[],withintimewindow{i},timewindows{i},{ha(indintoha) ha(indintoha+1)});
                        strtoplot=[physthenphoto_fields{j} ' aligned to ' alignments{i}];
                        f_=regexp(strtoplot,'_');
                        strtoplot(f_)=' ';
                        text(nanmin(xranges{i}),nanmax(nanmean(dataout.y,1)),strtoplot);
                        indintoha=indintoha+2;
                        closeAllBut(mainfig1);
                    else
                        continue
                    end
                end
            end
            if strcmp(physthenphoto_fields{j},'unit_by_unit') && k==1
                % replace with unit by unit plot
                axes(ha(plottedIntoWhichAxes));
                cla(ha(plottedIntoWhichAxes));
                if normalizeSU==true
                    alignmentCompanion.y=alignmentCompanion.y-nanmin(nanmean(alignmentCompanion.y,1));
                    alignmentCompanion.y=length(su).*alignmentCompanion.y./nanmax(nanmean(alignmentCompanion.y,1));
                end
                plotWStderr(alignmentCompanion.y,alignmentCompanion.x,'g',[],size(alignmentCompanion.y,1));
                hold on;
                unitByUnitPlot(su,dataout.x,maxTrialLength,normalizeSU);
                strtoplot=[physthenphoto_fields{j} ' aligned to ' alignments{i}];
                f_=regexp(strtoplot,'_');
                strtoplot(f_)=' ';
                text(nanmin(xranges{i}),nanmax(nanmean(alignmentCompanion.y,1)),strtoplot);
            end
        end
    end
end

k=1;
for i=1:Nw
    % each alignment
    for j=1:Nh
        set(ha(k),'XLim',xranges{i});
        % spawn individual figures
        f=figure();
        newax=copyobj(ha(k),f);
        set(newax,'Position',[0.1 0.1 0.85 0.85]);
        k=k+1;
    end
end

return

% figure 2
alignments={'success_fromPerchOrWheel','success batch when pellet dislodged','drop_fromPerchOrWheel','misses_and_pelletMissing'};
timewindows={[-1 16],[-1 16],[-1 16],[-1 16]};
withintimewindow={'first','first','first'};
beh_fields={'all_reachBatch','fidgetData'};
photo_fields={'green_ch'};
phys_fields={'unit_by_unit','sum_across_SU','MU'};
isPhysField=[ones(size(1:length(phys_fields))) zeros(size(1:length(photo_fields)))];
% figure 3
alignments={'success_fromPerchOrWheel','drop_fromPerchOrWheel','misses_and_pelletMissing'};
timewindows={[0 5],[0 5],[0 5],[0 5]};
withintimewindow={'first','first','first'};
beh_fields={'all_reachBatch','fidgetData'};
photo_fields={'green_ch'};
phys_fields={'unit_by_unit','sum_across_SU','MU'};
isPhysField=[ones(size(1:length(phys_fields))) zeros(size(1:length(photo_fields)))];
% figure 4

end

function unitByUnitPlot(su,su_times,maxTrialLength,normalizeSU)

ds=5;
spaceBetween=0;
offset=0;
endPlotAtInd=find(downSampAv(su_times,ds)>maxTrialLength,1,'first')-1;
for i=1:length(su)
    temp=downSampMatrix(su(i).alignedData,ds);
    if normalizeSU==true
        % normalize to y range 0 to 1
        sumin=nanmin(nanmean(temp(:,1:endPlotAtInd),1));
        temp=temp-sumin;
        sumax=nanmax(nanmean(temp(:,1:endPlotAtInd),1));
        temp=temp./sumax;
    end
    thisismax=plotWStderr(temp,downSampAv(su_times,ds),'k',endPlotAtInd,size(temp,1),offset);
    hold on;
    offset=thisismax+spaceBetween;
end
ylim([0 offset]);

end

function thisismax=plotWStderr(varargin)

dataMatrix=varargin{1};
times=varargin{2};
c=varargin{3};
plotUntilInd=varargin{4};
nEvents=varargin{5};
if length(varargin)>5
    offset=varargin{6};
else
    offset=0;
end

showStdevInstead=false;

if isempty(plotUntilInd)
    plotUntilInd=size(dataMatrix,2);
end

plot(times(1:plotUntilInd),offset+nanmean(dataMatrix(:,1:plotUntilInd),1),'Color',c,'LineWidth',1);
hold on;
if showStdevInstead==true
    plot(times(1:plotUntilInd),offset+nanmean(dataMatrix(:,1:plotUntilInd),1)-nanstd(dataMatrix(:,1:plotUntilInd),[],1),'Color',c,'LineWidth',0.5);
    plot(times(1:plotUntilInd),offset+nanmean(dataMatrix(:,1:plotUntilInd),1)+nanstd(dataMatrix(:,1:plotUntilInd),[],1),'Color',c,'LineWidth',0.5);
    thisismax=offset+nanmax(nanmean(dataMatrix(:,1:plotUntilInd),1)+nanstd(dataMatrix(:,1:plotUntilInd),[],1));
else
    plot(times(1:plotUntilInd),offset+nanmean(dataMatrix(:,1:plotUntilInd),1)-nanstd(dataMatrix(:,1:plotUntilInd),[],1)./sqrt(nEvents),'Color',c,'LineWidth',0.5);
    plot(times(1:plotUntilInd),offset+nanmean(dataMatrix(:,1:plotUntilInd),1)+nanstd(dataMatrix(:,1:plotUntilInd),[],1)./sqrt(nEvents),'Color',c,'LineWidth',0.5);
    thisismax=offset+nanmax(nanmean(dataMatrix(:,1:plotUntilInd),1)+nanstd(dataMatrix(:,1:plotUntilInd),[],1)./sqrt(nEvents));
end
    
end

function closeAllBut(f)

figs2keep = f;
all_figs = findobj(0, 'type', 'figure');
delete(setdiff(all_figs, figs2keep));

end