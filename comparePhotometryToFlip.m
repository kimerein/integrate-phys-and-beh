function comparePhotometryToFlip(response_to_plot,data_loc_array,whichSessionsToLoad,shiftPhotoBack)

% Is it a flip or a reshuffle? idk

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

% Check photometry and unit trial alignment
if max(Response_photo.fromWhichTrial)~=max(Response.fromWhichTrial)
    disp('Trial lengths of photo and units do not match');
    pause;
end
[~,si]=unique(Response_photo.fromWhichTrial);
[~,siu]=unique(Response.fromWhichTrial);
figure(); plot(si,Response_photo.isEventInThisTrial(si),'Color','r'); hold on;
plot(siu,Response.isEventInThisTrial(siu),'Color','b');
title('is event in this trial');
if ~isempty(shiftPhotoBack)
    Response_photo.fromWhichTrial=Response_photo.fromWhichTrial+shiftPhotoBack;
    [rsi,si]=unique(Response_photo.fromWhichTrial);
    [rsiu,siu]=unique(Response.fromWhichTrial);
    figure(); plot(rsi,Response_photo.isEventInThisTrial(si),'Color','r'); hold on;
    plot(rsiu,Response.isEventInThisTrial(siu),'Color','b');
    title('is event in this trial AFTER PHOTO SHIFT BACK');
end

% Plot simple averages
ds=5; figure(); plot(nanmean(Response_photo.unitbyunit_x,1),nanmean(Response_photo.unitbyunit_y,1),'Color','r');
hold on; plot(downSampAv(nanmean(Response.unitbyunit_x,1),ds),downSampAv(nanmean(Response.unitbyunit_y,1),ds),'Color','b');
hold on; plot(nanmean(Response.aligncomp_x,1),nanmean(Response.aligncomp_y,1),'Color','k');
legend({'DA','Average SU FR across trials','aligncomp'});

% Relative to DA after aligncomp peak, what is timing of first next spike?
ui=4;
thresh=0.5; % DA peak thresh
spithresh=0.05; 
ds=5;
f=find(Response_photo.isEventInThisTrial==1);
dapeaks=nan(1,length(f));
firstspiketimes=nan(1,length(f));
temp=nanmean(Response.aligncomp_x,1);
aligncomp_time=temp(find(nanmean(Response.aligncomp_y,1)>0.5,1,'first'));
figure(); plot(downSampAv(nanmean(Response.unitbyunit_x(Response.fromWhichUnit==ui,:),1),ds),downSampAv(nanmean(Response.unitbyunit_y(Response.fromWhichUnit==ui,:),1),ds),'Color','b');
hold on; plot(nanmean(Response_photo.unitbyunit_x,1),nanmean(Response_photo.unitbyunit_y,1),'Color','r');
hold on; plot(nanmean(Response.aligncomp_x,1),nanmean(Response.aligncomp_y,1),'Color','k');
figure(); plot(nanmean(Response_photo.unitbyunit_x,1),Response_photo.unitbyunit_y');
for i=1:length(f)
    tri=Response_photo.fromWhichTrial(f(i));
    DAtimes=nanmean(Response_photo.unitbyunit_x,1);
    spiketimes=downSampAv(nanmean(Response.unitbyunit_x,1),ds);
    currDA=smooth(Response_photo.unitbyunit_y(Response_photo.fromWhichTrial==tri,:),10);
    currspiking=downSampAv(Response.unitbyunit_y(Response.fromWhichTrial==tri & Response.fromWhichUnit==ui,:),ds);
    if isempty(currspiking)
        continue
    end
    if isempty(currDA)
        continue
    end
    % find first DA deflection above thresh
    fda=find(currDA(DAtimes>aligncomp_time)>thresh,1,'first');
    if ~isempty(fda)
        DAtimes=DAtimes(DAtimes>aligncomp_time);
        dapeaks(i)=DAtimes(fda);
    end
    % find time of first spike
    fda=find(currspiking(spiketimes>aligncomp_time)>spithresh,1,'first');
    if ~isempty(fda)
        spiketimes=spiketimes(spiketimes>aligncomp_time);
        firstspiketimes(i)=spiketimes(fda);
    end
end
figure();
scatter(dapeaks,firstspiketimes);
xlabel('da peak time (sec)');
ylabel('first spike time (sec)');

% Test SU pop vec flip first


end
