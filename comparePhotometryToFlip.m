function comparePhotometryToFlip(response_to_plot,data_loc_array,whichSessionsToLoad,dapeakoffset)

% dapeakoffset is time delay in seconds between aligncomp peak and START of
% success-related DA peak
% negative means DA peak comes BEFORE aligncomp peak

% Is it a flip or a reshuffle? idk
tset=settingsForStriatumUnitPlots;
if tset.keepAllSingleTrials==false
    error('settingsForStriatumUnitPlots.keepAllSingleTrials must be true');
end

% Load units
dd=cell(1,length(whichSessionsToLoad));
for i=1:length(whichSessionsToLoad)
    % load locations of SU data aligned to behavior
    % e.g., 'Z:\MICROSCOPE\Kim\WHISPER recs\Mar_3\20210721\SU aligned to behavior';
    dd{i}=data_loc_array{whichSessionsToLoad(i),8};
end
% choose type of response to plot
% response_to_plot='cue'; % can be any of the directories created in saveBehaviorAlignmentsSingleNeuron.m
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

% Load photometry
dd_photo=cell(1,length(whichSessionsToLoad));
for i=length(whichSessionsToLoad)
    % load locations of photo data aligned to behavior
    % e.g., 'Z:\MICROSCOPE\Kim\WHISPER recs\Mar_3\20210721\photometry aligned to behavior';
    dd_photo{i}=data_loc_array{whichSessionsToLoad(i),15};
end
% choose type of response to plot
% response_to_plot='uncued_success'; % can be any of the directories created in saveBehaviorAlignmentsSingleNeuron.m
dd_photo_more=cell(1,length(dd_photo)); 
for i=1:length(dd_photo)
    dd_photo_more{i}=[dd_photo{i} sep response_to_plot];
end
Response_photo=getAndSaveResponse(dd_photo_more,'_',settingsForStriatumPhotometryPlots,[]);

if length(dd)>1
    error('Not yet implemented for multiple sessions');
end
a=load([regexprep(dd{1}, 'SU aligned to behavior', 'tbt', 'ignorecase') '\behavior_tbt.mat']);
b=load([regexprep(dd{1}, 'SU aligned to behavior', 'tbt', 'ignorecase') '\beh2_tbt.mat']);
thesePhotoTrials=1:size(a.behavior_tbt.cue,1);
referToThesePhysTrials=a.behavior_tbt.reference_into_beh2trialinds(:,1);

% Plot simple averages
ds=5; figure(); plot(nanmean(Response_photo.unitbyunit_x,1),nanmean(Response_photo.unitbyunit_y,1),'Color','r');
hold on; plot(downSampAv(nanmean(Response.unitbyunit_x,1),ds),downSampAv(nanmean(Response.unitbyunit_y,1),ds),'Color','b');
hold on; plot(nanmean(Response.aligncomp_x,1),nanmean(Response.aligncomp_y,1),'Color','k');
legend({'DA','Average SU FR across trials','aligncomp'});

f=find(Response_photo.isEventInThisTrial==1);
figure();
plot(nanmean(Response_photo.unitbyunit_x,1),Response_photo.unitbyunit_y(f(randperm(length(f),5)),:)');
hold on; plot(nanmean(Response.aligncomp_x,1),nanmean(Response.aligncomp_y,1),'Color','k');

dathresh = input('DA event thresh: ');
thresh=dathresh;

% Relative to DA after aligncomp peak, what is timing of first next spike?
uniunits=unique(Response.fromWhichUnit);
for uniu=1:length(uniunits)
    ui=uniunits(uniu);
    %thresh=2; % DA peak thresh
    spithresh=0.05;
    ds=5;
    f=find(Response_photo.isEventInThisTrial==1);
    dapeaks=nan(1,length(f));
    firstspiketimes=nan(1,length(f));
    temp=nanmean(Response.aligncomp_x,1);
    aligncomp_time=temp(find(nanmean(Response.aligncomp_y,1)>0.5,1,'first'));

    figure();
    subplot(2,1,1);
    plot(downSampAv(nanmean(Response.unitbyunit_x(Response.fromWhichUnit==ui,:),1),ds),downSampAv(nanmean(Response.unitbyunit_y(Response.fromWhichUnit==ui,:),1),ds),'Color','b');
    hold on; plot(nanmean(Response_photo.unitbyunit_x,1),nanmean(Response_photo.unitbyunit_y,1),'Color','r');
    hold on; plot(nanmean(Response.aligncomp_x,1),nanmean(Response.aligncomp_y,1),'Color','k');
    %figure(); plot(nanmean(Response_photo.unitbyunit_x,1),Response_photo.unitbyunit_y');
    for i=1:length(f)
        % Get which phys trial corresponds to this photometry trial
        tri=Response_photo.fromWhichTrial(f(i));
        if length(unique(Response_photo.fromWhichTrial))~=length(thesePhotoTrials)
            error('inconsistent number of photometry trials');
        end
        phys_tri=referToThesePhysTrials(f(i));

        DAtimes=nanmean(Response_photo.unitbyunit_x,1);
        spiketimes=downSampAv(nanmean(Response.unitbyunit_x,1),ds);
        currDA=smooth(Response_photo.unitbyunit_y(Response_photo.fromWhichTrial==tri,:),10);
        currspiking=downSampAv(Response.unitbyunit_y(Response.fromWhichTrial==phys_tri & Response.fromWhichUnit==ui,:),ds);
        if isempty(currspiking)
            continue
        end
        if isempty(currDA)
            continue
        end
        % find first DA deflection above thresh
        fda=find(currDA(DAtimes>aligncomp_time+dapeakoffset)>thresh,1,'first');
        if ~isempty(fda)
            DAtimes=DAtimes(DAtimes>aligncomp_time+dapeakoffset);
            dapeaks(i)=DAtimes(fda);
        end
        % find time of first spike
        fda=find(currspiking(spiketimes>aligncomp_time+dapeakoffset)>spithresh,1,'first');
        if ~isempty(fda)
            spiketimes=spiketimes(spiketimes>aligncomp_time+dapeakoffset);
            firstspiketimes(i)=spiketimes(fda);
        end
    end
    subplot(2,1,2);
    scatter(dapeaks,firstspiketimes);
    xlabel('da peak time (sec)');
    ylabel('first spike time (sec)');
    [r,p]=corrcoef(dapeaks(~isnan(dapeaks) & ~isnan(firstspiketimes)),firstspiketimes(~isnan(dapeaks) & ~isnan(firstspiketimes)));
    title(['unit ' num2str(ui) ' r ' num2str(r(1,2)) ' p ' num2str(p(1,2))]);
end

% Test SU pop vec flip first


end