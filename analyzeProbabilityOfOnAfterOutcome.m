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
forscatter.p_success_unitbyunit=p_success_unitbyunit;
forscatter.p_failure_unitbyunit=p_failure_unitbyunit;

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

[isSig,pval]=getSignificantUnits(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure,nBoots);
figure(); scatter(p_success_unitbyunit(~isSig),p_failure_unitbyunit(~isSig)); hold on; scatter(p_success_unitbyunit(isSig),p_failure_unitbyunit(isSig),[],'filled');
xlabel('p success'); ylabel('p failure');
binedges=-1-0.025:0.05:1; binedges(1)=-1.001; binedges(end)=1.001;
[n,x]=histcounts(p_success_unitbyunit(~isSig)-p_failure_unitbyunit(~isSig),binedges); [n,x]=cityscape_hist(n,x);
figure(); plot(x,n,'Color',[0.5 0.5 0.5]);
[n,x]=histcounts(p_success_unitbyunit(isSig)-p_failure_unitbyunit(isSig),binedges); [n,x]=cityscape_hist(n,x);
hold on; plot(x,n,'Color','k'); xlim([-1 1]);
% plot dprimes
maxdp=3;
binedges=-maxdp:0.1:maxdp; 
dp=norminv(p_success_unitbyunit)-norminv(p_failure_unitbyunit);
dp(dp<-maxdp)=-maxdp; dp(dp>maxdp)=maxdp;
[n,x]=histcounts(dp(~isSig),binedges); [n,x]=cityscape_hist(n,x);
figure(); plot(x,n,'Color',[0.5 0.5 0.5]);
[n,x]=histcounts(dp(isSig),binedges); [n,x]=cityscape_hist(n,x);
hold on; plot(x,n,'Color','k'); xlim([-maxdp maxdp]);
% sig from how many sess
disp(['sig units from ' num2str(length(unique(success_Response.fromWhichSess(isSig)))) ' out of ' num2str(length(unique(success_Response.fromWhichSess))) ' sessions']);


% get change in unit probability 
% sigs and plots
[isSig,pval]=getSignificantUnits_Differences(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure,nBoots);
figure(); scatter(p_success_unitbyunit(~isSig),p_failure_unitbyunit(~isSig)); hold on; scatter(p_success_unitbyunit(isSig),p_failure_unitbyunit(isSig),[],'filled');
xlabel('p success'); ylabel('p failure');
binedges=-1-0.025:0.05:1; binedges(1)=-1.001; binedges(end)=1.001;
[n,x]=histcounts(p_success_unitbyunit(~isSig)-p_failure_unitbyunit(~isSig),binedges); [n,x]=cityscape_hist(n,x);
figure(); plot(x,n,'Color',[0.5 0.5 0.5]);
[n,x]=histcounts(p_success_unitbyunit(isSig)-p_failure_unitbyunit(isSig),binedges); [n,x]=cityscape_hist(n,x);
hold on; plot(x,n,'Color','k'); xlim([-1 1]);
% plot dprimes
maxdp=3;
binedges=-maxdp:0.1:maxdp; 
dp=norminv(p_success_unitbyunit)-norminv(p_failure_unitbyunit);
dp(dp<-maxdp)=-maxdp; dp(dp>maxdp)=maxdp;
[n,x]=histcounts(dp(~isSig),binedges); [n,x]=cityscape_hist(n,x);
figure(); plot(x,n,'Color',[0.5 0.5 0.5]);
[n,x]=histcounts(dp(isSig),binedges); [n,x]=cityscape_hist(n,x);
hold on; plot(x,n,'Color','k'); xlim([-maxdp maxdp]);
% sig from how many sess
disp(['sig units from ' num2str(length(unique(success_Response.fromWhichSess(isSig)))) ' out of ' num2str(length(unique(success_Response.fromWhichSess))) ' sessions']);

% changes across time
[isSig,pval,real_success_per_unit,real_failure_per_unit]=getSignificantUnits_Differences(units,success_Response,failure_Response,100);
figure(); scatter(real_success_per_unit(~isSig),real_failure_per_unit(~isSig)); hold on; scatter(real_success_per_unit(isSig),real_failure_per_unit(isSig),[],'filled');
xlabel('change p success'); ylabel('change p failure');
binedges=-1-0.025:0.05:1; binedges(1)=-1.001; binedges(end)=1.001;
[n,x]=histcounts(real_success_per_unit(~isSig)-real_failure_per_unit(~isSig),binedges); [n,x]=cityscape_hist(n,x);
figure(); plot(x,n,'Color',[0.5 0.5 0.5]);
[n,x]=histcounts(real_success_per_unit(isSig)-real_failure_per_unit(isSig),binedges); [n,x]=cityscape_hist(n,x);
hold on; plot(x,n,'Color','k'); xlim([-1 1]);
% plot dprimes
maxdp=3;
binedges=-maxdp:0.1:maxdp; 
dp=norminv(real_success_per_unit)-norminv(real_failure_per_unit);
dp(dp<-maxdp)=-maxdp; dp(dp>maxdp)=maxdp;
[n,x]=histcounts(dp(~isSig),binedges); [n,x]=cityscape_hist(n,x);
figure(); plot(x,n,'Color',[0.5 0.5 0.5]);
[n,x]=histcounts(dp(isSig),binedges); [n,x]=cityscape_hist(n,x);
hold on; plot(x,n,'Color','k'); xlim([-maxdp maxdp]);
% sig from how many sess
disp(['sig units from ' num2str(length(unique(success_Response.fromWhichSess(isSig)))) ' out of ' num2str(length(unique(success_Response.fromWhichSess))) ' sessions']);
figure(); scatter(abs(real_success_per_unit(~isSig)),abs(real_failure_per_unit(~isSig))); hold on; scatter(abs(real_success_per_unit(isSig)),abs(real_failure_per_unit(isSig)),[],'filled');
xlabel('change p success'); ylabel('change p failure');

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

% get 95% CI for bootstrap (i.e., shuffle)
p_bins=0:0.01:1.001;
p_bins(end)=1.001;
perc2point5=nan(length(p_bins)-1,1);
perc97point5=nan(length(p_bins)-1,1);
binmids=nan(length(p_bins)-1,1);
for i=1:length(p_bins)-1
    successbin=[p_bins(i) p_bins(i+1)];
    binmids(i)=mean(successbin);
    P=prctile(p_failure_unitbyunit(p_success_unitbyunit>=successbin(1) & p_success_unitbyunit<successbin(2)),[2.5 97.5]);
    perc2point5(i)=P(1);
    perc97point5(i)=P(2);
end
figure();
scatter(p_success_unitbyunit,p_failure_unitbyunit); xlabel('p success SHUFFLE'); ylabel('p failure SHUFFLE'); hold on;
plot(binmids,perc2point5,'Color',[0.8 0.8 0.8]); plot(binmids,perc97point5,'Color',[0.8 0.8 0.8]);
L1=polyfit(binmids(~isnan(perc2point5)),perc2point5(~isnan(perc2point5)),1); yfit1=L1(1)*binmids+L1(2);
L2=polyfit(binmids(~isnan(perc97point5)),perc97point5(~isnan(perc97point5)),1); yfit2=L2(1)*binmids+L2(2);
plot(binmids,yfit1,'Color','k'); plot(binmids,yfit2,'Color','k'); xlim([0 1]); ylim([0 1]);
% [b,bint,r,rint,stats]=regress(p_failure_unitbyunit,[ones(size(p_success_unitbyunit)) p_success_unitbyunit],0.05);
% figure();
% scatter(p_success_unitbyunit,p_failure_unitbyunit); xlabel('p success SHUFFLE'); ylabel('p failure SHUFFLE'); hold on;
% plot(binmids,binmids*bint(2,1)+bint(1,1)); plot(binmids,binmids*bint(2,2)+bint(1,2));

figure();
scatter(forscatter.p_success_unitbyunit,forscatter.p_failure_unitbyunit);
xlabel('p success');
ylabel('p failure');
hold on; plot(binmids,yfit1,'Color','k'); plot(binmids,yfit2,'Color','k'); xlim([0 1]); ylim([0 1]);
% [out2point5,out97point5,outbinmids]=bootstrapPerctiles(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure,100);
% figure();
% scatter(forscatter.p_success_unitbyunit,forscatter.p_failure_unitbyunit);
% xlabel('p success');
% ylabel('p failure');
% hold on; plot(outbinmids,out2point5,'Color','k'); plot(outbinmids,out97point5,'Color','k'); xlim([0 1]); ylim([0 1]);
% count how many units are differentially modulated according to 95% CI
p_bins=0:0.01:1.001;
p_bins(end)=1.001;
nabove=nan(length(p_bins)-1,1);
nbelow=nan(length(p_bins)-1,1);
nwithin=nan(length(p_bins)-1,1);
for i=1:length(p_bins)-1
    successbin=[p_bins(i) p_bins(i+1)];
    temp=forscatter.p_failure_unitbyunit(forscatter.p_success_unitbyunit>=successbin(1) & forscatter.p_success_unitbyunit<successbin(2));
    nabove(i)=nansum(temp>=yfit2(i));
    nbelow(i)=nansum(temp<=yfit1(i));
    nwithin(i)=nansum(temp>yfit1(i) & temp<yfit2(i));
end
disp([num2str(sum(nabove)) ' units above 95% CI']);
disp([num2str(sum(nbelow)) ' units below 95% CI']);
disp([num2str(sum(nwithin)) ' units within 95% CI']);

backup_success_Response=success_Response; backup_failure_Response=failure_Response;
trialNCutoff=6;
out=plotVariousSUResponsesAlignedToBeh('matchUnitsAcrossResponses',excludeTooFewTrials(success_Response,trialNCutoff,false),excludeTooFewTrials(failure_Response,trialNCutoff,false),[],[],[]);
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

% figure(); scatter(p_success_unitbyunitBEFORE(success_Response.A2atag(success_Response.excluded==0)==1)-p_success_unitbyunitAFTER(success_Response.A2atag(success_Response.excluded==0)==1),p_failure_unitbyunitBEFORE(success_Response.A2atag(success_Response.excluded==0)==1)-p_failure_unitbyunitAFTER(success_Response.A2atag(success_Response.excluded==0)==1),[],'b');
% hold on; scatter(p_success_unitbyunitBEFORE(success_Response.D1tag(success_Response.excluded==0)==1)-p_success_unitbyunitAFTER(success_Response.D1tag(success_Response.excluded==0)==1),p_failure_unitbyunitBEFORE(success_Response.D1tag(success_Response.excluded==0)==1)-p_failure_unitbyunitAFTER(success_Response.D1tag(success_Response.excluded==0)==1),[],'r');

end

function [isSig,pval,real_success_per_unit,real_failure_per_unit]=getSignificantUnits_Differences(units,success_Response,failure_Response,nBoots)

boot_success_per_unit=nan(length(units),nBoots);
boot_failure_per_unit=nan(length(units),nBoots);
% get real prob values
timeWindow=[-2 0];
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
[p_success_unitbyunit_BEFORE,p_failure_unitbyunit_BEFORE]=getProbOfResponse(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure);
timeWindow=[1 3];
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
[p_success_unitbyunit_AFTER,p_failure_unitbyunit_AFTER]=getProbOfResponse(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure);
real_success_per_unit=p_success_unitbyunit_BEFORE-p_success_unitbyunit_AFTER;
real_failure_per_unit=p_failure_unitbyunit_BEFORE-p_failure_unitbyunit_AFTER;
for j=1:nBoots
    % shuffle and bootstrap
    timeWindow=[-2 0];
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
    % shuffle
    randie=cell(length(units),1);
    for i=1:length(units)
        currunit=units(i);
        resp=[unitfr_success(fromWhichUnit_success==currunit); unitfr_failure(fromWhichUnit_failure==currunit)];
        randie{i}=randperm(length(resp));
        shuffle_resp=resp(randie{i});
        howmanysuccesses=nansum(fromWhichUnit_success==currunit);
        unitfr_success(fromWhichUnit_success==currunit)=shuffle_resp(1:howmanysuccesses);
        unitfr_failure(fromWhichUnit_failure==currunit)=shuffle_resp(howmanysuccesses+1:end);
    end
    % get prob
    [p_success_unitbyunit_BEFORE,p_failure_unitbyunit_BEFORE]=getProbOfResponse(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure);
    % new time window
    timeWindow=[1 3];
    [~,startAt]=nanmin(abs(nanmean(success_Response.unitbyunit_x,1)-(alignSuccessTime+timeWindow(1))));
    [~,endAt]=nanmin(abs(nanmean(success_Response.unitbyunit_x,1)-(alignSuccessTime+timeWindow(2))));
    successRange=[startAt endAt];
    [~,startAt]=nanmin(abs(nanmean(failure_Response.unitbyunit_x,1)-(alignFailureTime+timeWindow(1))));
    [~,endAt]=nanmin(abs(nanmean(failure_Response.unitbyunit_x,1)-(alignFailureTime+timeWindow(2))));
    failureRange=[startAt endAt];
    unitfr_success=sum(success_Response.unitbyunit_y(:,successRange(1):successRange(2)),2,'omitnan');
    unitfr_success=unitfr_success(success_Response.isEventInThisTrial==1);
    unitfr_failure=sum(failure_Response.unitbyunit_y(:,failureRange(1):failureRange(2)),2,'omitnan');
    unitfr_failure=unitfr_failure(failure_Response.isEventInThisTrial==1);
    % shuffle
    for i=1:length(units)
        currunit=units(i);
        resp=[unitfr_success(fromWhichUnit_success==currunit); unitfr_failure(fromWhichUnit_failure==currunit)];
        shuffle_resp=resp(randie{i});
        howmanysuccesses=nansum(fromWhichUnit_success==currunit);
        unitfr_success(fromWhichUnit_success==currunit)=shuffle_resp(1:howmanysuccesses);
        unitfr_failure(fromWhichUnit_failure==currunit)=shuffle_resp(howmanysuccesses+1:end);
    end
    % get prob
    [p_success_unitbyunit_AFTER,p_failure_unitbyunit_AFTER]=getProbOfResponse(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure);
    boot_success_per_unit(:,j)=p_success_unitbyunit_BEFORE-p_success_unitbyunit_AFTER;
    boot_failure_per_unit(:,j)=p_failure_unitbyunit_BEFORE-p_failure_unitbyunit_AFTER;
end
% get significance for each unit
% is significant if real difference between success and failure is greater
% than 97.5th percentile of shuffled differences OR less than 2.5th
% percentile of shuffled differences
greaterthanwhichprctile=nan(length(units),1);
lessthanwhichprctile=nan(length(units),1);
pval=nan(length(units),1);
getptiles=1:100;
for i=1:length(units)
    ptiles=prctile(boot_success_per_unit(i,:)-boot_failure_per_unit(i,:),getptiles);
    temp=getptiles(find(ptiles<(real_success_per_unit(i)-real_failure_per_unit(i)),1,'last'));
    if ~isempty(temp)
        greaterthanwhichprctile(i)=temp;
    else
        greaterthanwhichprctile(i)=0;
    end
    temp=getptiles(find(ptiles>(real_success_per_unit(i)-real_failure_per_unit(i)),1,'first'));
    if ~isempty(temp)
        lessthanwhichprctile(i)=temp;
    else
        lessthanwhichprctile(i)=100;
    end
    if greaterthanwhichprctile(i)>50
        if greaterthanwhichprctile(i)==100
            pval(i)=1/100;
            continue
        end
        % above median
        pval(i)=((100-(greaterthanwhichprctile(i)+1))/100);
    else
        if lessthanwhichprctile(i)==0
            pval(i)=1/100;
            continue
        end
        % below median
        pval(i)=(((lessthanwhichprctile(i)-1)-0)/100);
    end
end
isSig=pval<0.025;

end

function [isSig,pval]=getSignificantUnits(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure,nBoots)

boot_success_per_unit=nan(length(units),nBoots);
boot_failure_per_unit=nan(length(units),nBoots);
% get real prob values
[real_success_per_unit,real_failure_per_unit]=getProbOfResponse(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure);
for j=1:nBoots
    % shuffle and bootstrap
    [unitfr_success,unitfr_failure]=shuffleTrialType(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure);
    [boot_success_per_unit(:,j),boot_failure_per_unit(:,j)]=getProbOfResponse(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure);
end
% get significance for each unit
% is significant if real difference between success and failure is greater
% than 97.5th percentile of shuffled differences OR less than 2.5th
% percentile of shuffled differences
greaterthanwhichprctile=nan(length(units),1);
lessthanwhichprctile=nan(length(units),1);
pval=nan(length(units),1);
getptiles=1:100;
for i=1:length(units)
    ptiles=prctile(boot_success_per_unit(i,:)-boot_failure_per_unit(i,:),getptiles);
    temp=getptiles(find(ptiles<(real_success_per_unit(i)-real_failure_per_unit(i)),1,'last'));
    if ~isempty(temp)
        greaterthanwhichprctile(i)=temp;
    else
        greaterthanwhichprctile(i)=0;
    end
    temp=getptiles(find(ptiles>(real_success_per_unit(i)-real_failure_per_unit(i)),1,'first'));
    if ~isempty(temp)
        lessthanwhichprctile(i)=temp;
    else
        lessthanwhichprctile(i)=100;
    end
    if greaterthanwhichprctile(i)>50
        if greaterthanwhichprctile(i)==100
            pval(i)=1/100;
            continue
        end
        % above median
        pval(i)=((100-(greaterthanwhichprctile(i)+1))/100);
    else
        if lessthanwhichprctile(i)==0
            pval(i)=1/100;
            continue
        end
        % below median
        pval(i)=(((lessthanwhichprctile(i)-1)-0)/100);
    end
end
isSig=pval<0.025;

end

function [out2point5,out97point5,outbinmids]=bootstrapPerctiles(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure,nBoots)

p_bins=0:0.01:1.001;
p_bins(end)=1.001;
perc2point5=nan(length(p_bins)-1,nBoots);
perc97point5=nan(length(p_bins)-1,nBoots);
binmids=nan(length(p_bins)-1,nBoots);
for j=1:nBoots
    [unitfr_success,unitfr_failure]=shuffleTrialType(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure);
    [p_success_unitbyunit,p_failure_unitbyunit]=getProbOfResponse(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure);
    % get 95% CI for bootstrap (i.e., shuffle)
    for i=1:length(p_bins)-1
        successbin=[p_bins(i) p_bins(i+1)];
        binmids(i,j)=mean(successbin);
        P=prctile(p_failure_unitbyunit(p_success_unitbyunit>=successbin(1) & p_success_unitbyunit<successbin(2)),[2.5 97.5]);
        perc2point5(i,j)=P(1);
        perc97point5(i,j)=P(2);
    end
end
out2point5=mean(perc2point5,2,'omitnan');
out97point5=mean(perc97point5,2,'omitnan');
outbinmids=mean(binmids,2,'omitnan');

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



