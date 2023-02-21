function analyzeProbabilityOfOnAfterOutcome(dd, timeWindow, success_Response, failure_Response)

% timeWindow is in seconds wrt peak of aligncomp

%aboveBaseHz=baseFiringRate*(timeWindow(2)-timeWindow(1));

if isempty(success_Response)
    % choose type of response to plot
    response_to_plot='all_success'; % can be any of the directories created in saveBehaviorAlignmentsSingleNeuron.m
    % [inStructure isFS isTAN isSPN isLowFRThin]
    plotUnitCriteria=[1 0 0 1 0]; % -100 is a wildcard, else 0 (false) and 1 (true)
    getCriteriaForUnitsToPlot(plotUnitCriteria);
    % read in some units
    dd_more=cell(1,length(dd));
    for i=1:length(dd)
        dd_more{i}=[dd{i} sep response_to_plot];
    end
    whichUnitsToGrab='_'; % '_' for all units, or can be something like 'D1tagged'
    settings=settingsForStriatumUnitPlots;
    settings.maxUnitsPerSess=30;
    settings.keepAllSingleTrials=true;
    success_Response=getAndSaveResponse(dd_more,whichUnitsToGrab,settings,[]);
end
if isempty(failure_Response)
    % choose type of response to plot
    response_to_plot='all_drop'; % can be any of the directories created in saveBehaviorAlignmentsSingleNeuron.m
    % read in some units
    dd_more=cell(1,length(dd));
    for i=1:length(dd)
        dd_more{i}=[dd{i} sep response_to_plot];
    end
    failure_Response=getAndSaveResponse(dd_more,whichUnitsToGrab,settings,[]);
end

% get probability of response in time window for each unit
[~,alignPeakInd]=nanmax(nanmean(success_Response.aligncomp_y,1));
temp=nanmean(success_Response.aligncomp_x,1);
alignSuccessTime=temp(alignPeakInd);
[~,alignPeakInd]=nanmax(nanmean(failure_Response.aligncomp_y,1));
temp=nanmean(failure_Response.aligncomp_x,1);
alignFailureTime=temp(alignPeakInd);
[~,startAt]=nanmin(abs(nanmean(success_Response.unitbyunit_x,1)-(alignSuccessTime+timeWindow(1))));
[~,endAt]=nanmin(abs(nanmean(success_Response.unitbyunit_x,1)-(alignSuccessTime+timeWindow(2))));
successRange=[startAt endAt];
[~,startAt]=nanmin(abs(nanmean(failure_Response.unitbyunit_x,1)-(alignFailureTime+timeWindow(1))));
[~,endAt]=nanmin(abs(nanmean(failure_Response.unitbyunit_x,1)-(alignFailureTime+timeWindow(2))));
failureRange=[startAt endAt];

unitfr_success=sum(success_Response.unitbyunit_y(:,successRange(1):successRange(2)),2,'omitnan');
fromWhichUnit_success=success_Response.fromWhichUnit;
unitfr_success=unitfr_success(success_Response.isEventInThisTrial==1);
fromWhichUnit_success=fromWhichUnit_success(success_Response.isEventInThisTrial==1);
unitfr_failure=sum(failure_Response.unitbyunit_y(:,failureRange(1):failureRange(2)),2,'omitnan');
fromWhichUnit_failure=failure_Response.fromWhichUnit;
unitfr_failure=unitfr_failure(failure_Response.isEventInThisTrial==1);
fromWhichUnit_failure=fromWhichUnit_failure(failure_Response.isEventInThisTrial==1);

units=unique(success_Response.fromWhichUnit);
[p_success_unitbyunit,p_failure_unitbyunit]=getProbOfResponse(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure);

figure();
scatter(p_success_unitbyunit,p_failure_unitbyunit);
xlabel('p success');
ylabel('p failure');

p_bins=0:0.1:1.001;
p_bins(end)=1.001;
heatmap_p=nan(length(p_bins)-1,length(p_bins)-1);
for i=1:length(p_bins)-1
    for j=1:length(p_bins)-1
        successbin=[p_bins(i) p_bins(i+1)];
        failurebin=[p_bins(j) p_bins(j+1)];
        heatmap_p(i,j)=nansum(p_success_unitbyunit>=successbin(1) & p_success_unitbyunit<successbin(2) & p_failure_unitbyunit>=failurebin(1) & p_failure_unitbyunit<failurebin(2));
    end
end
figure();
imagesc(flipud(heatmap_p'));
ylabel('p failure');
xlabel('p success');

[unitfr_success,unitfr_failure]=shuffleTrialType(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure);
[p_success_unitbyunit,p_failure_unitbyunit]=getProbOfResponse(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure);

figure();
scatter(p_success_unitbyunit,p_failure_unitbyunit);
xlabel('p success SHUFFLE');
ylabel('p failure SHUFFLE');

p_bins=0:0.1:1.001;
p_bins(end)=1.001;
heatmap_p=nan(length(p_bins)-1,length(p_bins)-1);
for i=1:length(p_bins)-1
    for j=1:length(p_bins)-1
        successbin=[p_bins(i) p_bins(i+1)];
        failurebin=[p_bins(j) p_bins(j+1)];
        heatmap_p(i,j)=nansum(p_success_unitbyunit>=successbin(1) & p_success_unitbyunit<successbin(2) & p_failure_unitbyunit>=failurebin(1) & p_failure_unitbyunit<failurebin(2));
    end
end
figure();
imagesc(flipud(heatmap_p'));
ylabel('p failure SHUFFLE');
xlabel('p success SHUFFLE');

backup_success_Response=success_Response; backup_failure_Response=failure_Response;
out=plotVariousSUResponsesAlignedToBeh('matchUnitsAcrossResponses',excludeTooFewTrials(success_Response,6,false),excludeTooFewTrials(failure_Response,6,false),[],[],[]);
success_Response=out.Response1; failure_Response=out.Response2;
unitfr_success=sum(success_Response.unitbyunit_y(:,successRange(1):successRange(2)),2,'omitnan');
fromWhichUnit_success=success_Response.fromWhichUnit;
unitfr_success=unitfr_success(success_Response.isEventInThisTrial==1);
fromWhichUnit_success=fromWhichUnit_success(success_Response.isEventInThisTrial==1);
unitfr_failure=sum(failure_Response.unitbyunit_y(:,failureRange(1):failureRange(2)),2,'omitnan');
fromWhichUnit_failure=failure_Response.fromWhichUnit;
unitfr_failure=unitfr_failure(failure_Response.isEventInThisTrial==1);
fromWhichUnit_failure=fromWhichUnit_failure(failure_Response.isEventInThisTrial==1);
units=unique(success_Response.fromWhichUnit);
[p_success_unitbyunit,p_failure_unitbyunit]=getProbOfResponse(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure);
figure();
scatter(p_success_unitbyunit,p_failure_unitbyunit);
xlabel('p success');
ylabel('p failure');

end

function [p_success_unitbyunit,p_failure_unitbyunit]=getProbOfResponse(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure)

p_success_unitbyunit=nan(length(units),1);
p_failure_unitbyunit=nan(length(units),1);
disp([num2str(length(units)) ' units']);
for i=1:length(units)
    p_success_unitbyunit(i)=nansum(unitfr_success(fromWhichUnit_success==units(i))>0.5)./nansum(fromWhichUnit_success==units(i));
    p_failure_unitbyunit(i)=nansum(unitfr_failure(fromWhichUnit_failure==units(i))>0.5)./nansum(fromWhichUnit_failure==units(i));
end

end

function [unitfr_success,unitfr_failure]=shuffleTrialType(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure)

for i=1:length(units)
    currunit=units(i);
    resp=[unitfr_success(fromWhichUnit_success==currunit); unitfr_failure(fromWhichUnit_failure==currunit)];
    shuffle_resp=resp(randperm(length(resp)));
    howmanysuccesses=nansum(fromWhichUnit_success==currunit);
    unitfr_success(fromWhichUnit_success==currunit)=shuffle_resp(1:howmanysuccesses);
    unitfr_failure(fromWhichUnit_failure==currunit)=shuffle_resp(howmanysuccesses+1:end);
end

end



