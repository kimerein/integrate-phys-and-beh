function analyzeProbabilityOfOnAfterOutcome(dd, timeWindow, success_Response, failure_Response, response_to_plot1, response_to_plot2, overTimeOrJustTimeWindow)

% timeWindow is in seconds wrt peak of aligncomp

nBoots=100;

%aboveBaseHz=baseFiringRate*(timeWindow(2)-timeWindow(1));

if isempty(success_Response)
    % choose type of response to plot
    response_to_plot=response_to_plot1; %'cued_success'; %'all_success'; % can be any of the directories created in saveBehaviorAlignmentsSingleNeuron.m
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
    response_to_plot=response_to_plot2; %'uncued_success'; %'all_drop'; % can be any of the directories created in saveBehaviorAlignmentsSingleNeuron.m
    % read in some units
    dd_more=cell(1,length(dd));
    for i=1:length(dd)
        dd_more{i}=[dd{i} sep response_to_plot];
    end
    failure_Response=getAndSaveResponse(dd_more,whichUnitsToGrab,settings,[]);
end

% make aligncomp peaks the same
ti=nanmean(success_Response.aligncomp_x,1);
aligncompy=nanmean(success_Response.aligncomp_y,1);
[~,apeak]=nanmax(aligncompy); peakAt=ti(apeak);
success_Response.aligncomp_x=success_Response.aligncomp_x-peakAt;
success_Response.unitbyunit_x=success_Response.unitbyunit_x-peakAt;
ti=nanmean(failure_Response.aligncomp_x,1);
aligncompy=nanmean(failure_Response.aligncomp_y,1);
[~,apeak]=nanmax(aligncompy); peakAt=ti(apeak);
failure_Response.aligncomp_x=failure_Response.aligncomp_x-peakAt;
failure_Response.unitbyunit_x=failure_Response.unitbyunit_x-peakAt;

switch overTimeOrJustTimeWindow
    case 'justTimeWindow'
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

        [fr_success_unitbyunit,fr_failure_unitbyunit]=getFROfResponse(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure);
        figure(); scatter(fr_success_unitbyunit(~isSig)./200,fr_failure_unitbyunit(~isSig)./200); hold on; scatter(fr_success_unitbyunit(isSig)./200,fr_failure_unitbyunit(isSig)./200,[],'filled');
        xlabel('p success'); ylabel('p failure');

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

    case 'overTime'
        timeBins=nan(length(-3:0.125:9.5),2); 
        timeBins(:,1)=-3:0.125:9.5;
        timeBins(:,2)=[-3:0.125:9.5]+2;
        [dp,dpfr,isSig,pval,p_success_ubyu,p_failure_ubyu,fr_success_ubyu,fr_failure_ubyu,timebinMeans]=dprimesOverTime(success_Response,failure_Response,timeBins);
%         dp(isinf(dp) & dp>0)=10;
%         dp(isinf(dp) & dp<0)=-10;
%         dp(dp==-10)=nan; dp(dp==10)=nan;
%         load('Z:\MICROSCOPE\Kim\20230205 all SU alignments\all trials averaged not downsampled\groupLabelsFromTCA.mat');
%         figure(); plot(timebinMeans,nanmean(dp(:,groupLabelsFromTCA==1),2)); hold on;
%         plot(timebinMeans,nanmean(dp(:,groupLabelsFromTCA==1),2)+nanstd(dp(:,groupLabelsFromTCA==1),[],2)./sqrt(size(dp(:,groupLabelsFromTCA==1),1)));
%         plot(timebinMeans,nanmean(dp(:,groupLabelsFromTCA==1),2)-nanstd(dp(:,groupLabelsFromTCA==1),[],2)./sqrt(size(dp(:,groupLabelsFromTCA==1),1)));
%         cuez=nanmean(fr_success_ubyu(timebinMeans<0,:),1)-nanmean(fr_failure_ubyu(timebinMeans<0,:),1);
%         ta=cuez>prctile(cuez,90);
%         gpLab=2; figure(); plot(timebinMeans,nanmean(dp(:,groupLabelsFromTCA==gpLab & ta'),2),'Color','b'); hold on;
%         plot(timebinMeans,nanmean(dp(:,groupLabelsFromTCA==gpLab & ta'),2)+nanstd(dp(:,groupLabelsFromTCA==gpLab & ta'),[],2)./sqrt(size(dp(:,groupLabelsFromTCA==gpLab & ta'),1)),'Color','b');
%         plot(timebinMeans,nanmean(dp(:,groupLabelsFromTCA==gpLab & ta'),2)-nanstd(dp(:,groupLabelsFromTCA==gpLab & ta'),[],2)./sqrt(size(dp(:,groupLabelsFromTCA==gpLab & ta'),1)),'Color','b');
%         hold on; gpLab=1; plot(timebinMeans,nanmean(dp(:,groupLabelsFromTCA==gpLab & ta'),2),'Color','r'); hold on;
%         plot(timebinMeans,nanmean(dp(:,groupLabelsFromTCA==gpLab & ta'),2)+nanstd(dp(:,groupLabelsFromTCA==gpLab & ta'),[],2)./sqrt(size(dp(:,groupLabelsFromTCA==gpLab & ta'),1)),'Color','r');
%         plot(timebinMeans,nanmean(dp(:,groupLabelsFromTCA==gpLab & ta'),2)-nanstd(dp(:,groupLabelsFromTCA==gpLab & ta'),[],2)./sqrt(size(dp(:,groupLabelsFromTCA==gpLab & ta'),1)),'Color','r');

%         dpfr(dpfr>10)=10; dpfr(dpfr<-10)=10; dpfr(dpfr==-10)=nan; dpfr(dpfr==10)=nan;
%         gpLab=2; figure(); plot(timebinMeans,nanmean(dpfr(:,groupLabelsFromTCA==gpLab & ta'),2),'Color','b'); hold on;
%         plot(timebinMeans,nanmean(dpfr(:,groupLabelsFromTCA==gpLab & ta'),2)+nanstd(dpfr(:,groupLabelsFromTCA==gpLab & ta'),[],2)./sqrt(size(dpfr(:,groupLabelsFromTCA==gpLab & ta'),1)),'Color','b');
%         plot(timebinMeans,nanmean(dpfr(:,groupLabelsFromTCA==gpLab & ta'),2)-nanstd(dpfr(:,groupLabelsFromTCA==gpLab & ta'),[],2)./sqrt(size(dpfr(:,groupLabelsFromTCA==gpLab & ta'),1)),'Color','b');
%         hold on; gpLab=1; plot(timebinMeans,nanmean(dpfr(:,groupLabelsFromTCA==gpLab & ta'),2),'Color','r'); hold on;
%         plot(timebinMeans,nanmean(dpfr(:,groupLabelsFromTCA==gpLab & ta'),2)+nanstd(dpfr(:,groupLabelsFromTCA==gpLab & ta'),[],2)./sqrt(size(dpfr(:,groupLabelsFromTCA==gpLab & ta'),1)),'Color','r');
%         plot(timebinMeans,nanmean(dpfr(:,groupLabelsFromTCA==gpLab & ta'),2)-nanstd(dpfr(:,groupLabelsFromTCA==gpLab & ta'),[],2)./sqrt(size(dpfr(:,groupLabelsFromTCA==gpLab & ta'),1)),'Color','r');
% 
%         ta=cuez<prctile(cuez,10); gpLab=2; figure(); plot(timebinMeans,nanmean(dp(:,groupLabelsFromTCA==gpLab & ta'),2),'Color','b'); hold on; 
%         plot(timebinMeans,nanmean(dp(:,groupLabelsFromTCA==gpLab & ta'),2)+nanstd(dp(:,groupLabelsFromTCA==gpLab & ta'),[],2)./sqrt(size(dp(:,groupLabelsFromTCA==gpLab & ta'),1)),'Color','b'); 
%         plot(timebinMeans,nanmean(dp(:,groupLabelsFromTCA==gpLab & ta'),2)-nanstd(dp(:,groupLabelsFromTCA==gpLab & ta'),[],2)./sqrt(size(dp(:,groupLabelsFromTCA==gpLab & ta'),1)),'Color','b');
%         hold on; gpLab=1; plot(timebinMeans,nanmean(dp(:,groupLabelsFromTCA==gpLab & ta'),2),'Color','r'); hold on; 
%         plot(timebinMeans,nanmean(dp(:,groupLabelsFromTCA==gpLab & ta'),2)+nanstd(dp(:,groupLabelsFromTCA==gpLab & ta'),[],2)./sqrt(size(dp(:,groupLabelsFromTCA==gpLab & ta'),1)),'Color','r'); 
%         plot(timebinMeans,nanmean(dp(:,groupLabelsFromTCA==gpLab & ta'),2)-nanstd(dp(:,groupLabelsFromTCA==gpLab & ta'),[],2)./sqrt(size(dp(:,groupLabelsFromTCA==gpLab & ta'),1)),'Color','r');
    
    load('Z:\MICROSCOPE\Kim\20230205 all SU alignments\all trials averaged not downsampled\groupLabelsFromTCA.mat');
    backup.success_Response=success_Response;
    backup.failure_Response=failure_Response;
    timeBins=nan(length(-3:0.125:9.5),2); 
    timeBins(:,1)=-3:0.125:9.5;
    timeBins(:,2)=[-3:0.125:9.5]+2;
    ta=cuez>prctile(cuez,90) & groupLabelsFromTCA'==2;
    uns=unique(success_Response.fromWhichUnit);
    succFWU=success_Response.fromWhichUnit;
    failFWU=failure_Response.fromWhichUnit;
    success_Response.fromWhichUnit(ismember(succFWU,uns(ta)))=1;
    success_Response.fromWhichUnit(~ismember(succFWU,uns(ta)))=2;
    failure_Response.fromWhichUnit(ismember(failFWU,uns(ta)))=1;
    failure_Response.fromWhichUnit(~ismember(failFWU,uns(ta)))=2;
    
    combineTrials=true;
    if combineTrials==true
        newSuccess_Response=success_Response;
        newFailure_Response=failure_Response;
        usess=unique(success_Response.fromWhichSess_forTrials);
        k=1;
        for l=1:length(usess)
            thetrials=unique(success_Response.fromWhichTrial(success_Response.fromWhichSess_forTrials==usess(l)));
            for i=1:length(thetrials)
                currtrial=thetrials(i);
                for j=1:2 % 2 types of units
                    takeTrials=ismember(success_Response.fromWhichTrial,currtrial) & success_Response.fromWhichUnit==j & success_Response.fromWhichSess_forTrials==usess(l);
                    if nansum(takeTrials)>0
                        newSuccess_Response.unitbyunit_x(k,:)=nanmean(success_Response.unitbyunit_x(takeTrials,:),1);
                        newSuccess_Response.unitbyunit_y(k,:)=nanmean(success_Response.unitbyunit_y(takeTrials,:),1);
                        newSuccess_Response.aligncomp_x(k,:)=nanmean(success_Response.aligncomp_x(takeTrials,:),1);
                        newSuccess_Response.aligncomp_y(k,:)=nanmean(success_Response.aligncomp_y(takeTrials,:),1);
                        newSuccess_Response.fromWhichUnit(k,:)=nanmean(success_Response.fromWhichUnit(takeTrials,:),1);
                        newSuccess_Response.fromWhichTrial(k,:)=nanmean(success_Response.fromWhichTrial(takeTrials,:),1);
                        newSuccess_Response.isEventInThisTrial(k,:)=nanmean(success_Response.isEventInThisTrial(takeTrials,:),1);
                        newSuccess_Response.fromWhichSess_forTrials(k,:)=nanmean(success_Response.fromWhichSess_forTrials(takeTrials,:),1);
                        k=k+1;
                        if mod(k,100)==0
                            disp(k);
                        end
                    end
                end
            end
        end
        newSuccess_Response.unitbyunit_x=newSuccess_Response.unitbyunit_x(1:k-1,:);
        newSuccess_Response.unitbyunit_y=newSuccess_Response.unitbyunit_y(1:k-1,:);
        newSuccess_Response.aligncomp_x=newSuccess_Response.aligncomp_x(1:k-1,:);
        newSuccess_Response.aligncomp_y=newSuccess_Response.aligncomp_y(1:k-1,:);
        newSuccess_Response.fromWhichUnit=newSuccess_Response.fromWhichUnit(1:k-1,:);
        newSuccess_Response.fromWhichTrial=newSuccess_Response.fromWhichTrial(1:k-1,:);
        newSuccess_Response.isEventInThisTrial=newSuccess_Response.isEventInThisTrial(1:k-1,:);
        newSuccess_Response.fromWhichSess_forTrials=newSuccess_Response.fromWhichSess_forTrials(1:k-1,:);

        usess=unique(failure_Response.fromWhichSess_forTrials);
        k=1;
        for l=1:length(usess)
            thetrials=unique(failure_Response.fromWhichTrial(failure_Response.fromWhichSess_forTrials==usess(l)));
            for i=1:length(thetrials)
                currtrial=thetrials(i);
                for j=1:2 % 2 types of units
                    takeTrials=ismember(failure_Response.fromWhichTrial,currtrial) & failure_Response.fromWhichUnit==j & failure_Response.fromWhichSess_forTrials==usess(l);
                    if nansum(takeTrials)>0
                        newFailure_Response.unitbyunit_x(k,:)=nanmean(failure_Response.unitbyunit_x(takeTrials,:),1);
                        newFailure_Response.unitbyunit_y(k,:)=nanmean(failure_Response.unitbyunit_y(takeTrials,:),1);
                        newFailure_Response.aligncomp_x(k,:)=nanmean(failure_Response.aligncomp_x(takeTrials,:),1);
                        newFailure_Response.aligncomp_y(k,:)=nanmean(failure_Response.aligncomp_y(takeTrials,:),1);
                        newFailure_Response.fromWhichUnit(k,:)=nanmean(failure_Response.fromWhichUnit(takeTrials,:),1);
                        newFailure_Response.fromWhichTrial(k,:)=nanmean(failure_Response.fromWhichTrial(takeTrials,:),1);
                        newFailure_Response.isEventInThisTrial(k,:)=nanmean(failure_Response.isEventInThisTrial(takeTrials,:),1);
                        newFailure_Response.fromWhichSess_forTrials(k,:)=nanmean(failure_Response.fromWhichSess_forTrials(takeTrials,:),1);
                        k=k+1;
                    end
                end
            end
        end

        newFailure_Response.unitbyunit_x=newFailure_Response.unitbyunit_x(1:k-1,:);
        newFailure_Response.unitbyunit_y=newFailure_Response.unitbyunit_y(1:k-1,:);
        newFailure_Response.aligncomp_x=newFailure_Response.aligncomp_x(1:k-1,:);
        newFailure_Response.aligncomp_y=newFailure_Response.aligncomp_y(1:k-1,:);
        newFailure_Response.fromWhichUnit=newFailure_Response.fromWhichUnit(1:k-1,:);
        newFailure_Response.fromWhichTrial=newFailure_Response.fromWhichTrial(1:k-1,:);
        newFailure_Response.isEventInThisTrial=newFailure_Response.isEventInThisTrial(1:k-1,:);
        newFailure_Response.fromWhichSess_forTrials=newFailure_Response.fromWhichSess_forTrials(1:k-1,:);
        [dp,dpfr,isSig,pval,p_success_ubyu,p_failure_ubyu,fr_success_ubyu,fr_failure_ubyu,timebinMeans,fr_succ_sd,fr_fail_sd]=withinSessDprime(newSuccess_Response,newFailure_Response,timeBins);
    else
        [dp,dpfr,isSig,pval,p_success_ubyu,p_failure_ubyu,fr_success_ubyu,fr_failure_ubyu,timebinMeans,fr_succ_sd,fr_fail_sd]=withinSessDprime(success_Response,failure_Response,timeBins);
    end
    % axes of dp are now 1st dim (timebins) x 2nd dim (unit type) x 3rd dim
    % (n sessions)
    temp=p_success_ubyu(:,1,:); temp(isinf(temp))=nan;
%     temp=fr_success_ubyu(:,1,:); temp(isinf(temp))=nan;
    figure(); plot(timebinMeans,reshape(nanmean(temp(:,1,:),3),size(dp,1),1),'Color','k'); hold on;
    plot(timebinMeans,reshape(nanmean(temp(:,1,:),3),size(dp,1),1)+reshape(nanstd(temp(:,1,:),[],3)./sqrt(size(dp,3)),size(dp,1),1),'Color','k');
    plot(timebinMeans,reshape(nanmean(temp(:,1,:),3),size(dp,1),1)-reshape(nanstd(temp(:,1,:),[],3)./sqrt(size(dp,3)),size(dp,1),1),'Color','k');
    temp=p_success_ubyu(:,2,:); temp(isinf(temp))=nan;
%     temp=fr_success_ubyu(:,2,:); temp(isinf(temp))=nan;
    plot(timebinMeans,reshape(nanmean(temp(:,1,:),3),size(dp,1),1),'Color','m'); 
    plot(timebinMeans,reshape(nanmean(temp(:,1,:),3),size(dp,1),1)+reshape(nanstd(temp(:,1,:),[],3)./sqrt(size(dp,3)),size(dp,1),1),'Color','m');
    plot(timebinMeans,reshape(nanmean(temp(:,1,:),3),size(dp,1),1)-reshape(nanstd(temp(:,1,:),[],3)./sqrt(size(dp,3)),size(dp,1),1),'Color','m');
    legend({'ta true p active FIRST TRIAL TYPE','ta false p active FIRST TRIAL TYPE'});

%     temp=p_failure_ubyu(:,1,:); temp(isinf(temp))=nan;
    temp=fr_failure_ubyu(:,1,:); temp(isinf(temp))=nan;
    figure(); plot(timebinMeans,reshape(nanmean(temp(:,1,:),3),size(dp,1),1),'Color','k'); hold on;
    plot(timebinMeans,reshape(nanmean(temp(:,1,:),3),size(dp,1),1)+reshape(nanstd(temp(:,1,:),[],3)./sqrt(size(dp,3)),size(dp,1),1),'Color','k');
    plot(timebinMeans,reshape(nanmean(temp(:,1,:),3),size(dp,1),1)-reshape(nanstd(temp(:,1,:),[],3)./sqrt(size(dp,3)),size(dp,1),1),'Color','k');
%     temp=p_failure_ubyu(:,2,:); temp(isinf(temp))=nan;
    temp=fr_failure_ubyu(:,2,:); temp(isinf(temp))=nan;
    plot(timebinMeans,reshape(nanmean(temp(:,1,:),3),size(dp,1),1),'Color','m'); 
    plot(timebinMeans,reshape(nanmean(temp(:,1,:),3),size(dp,1),1)+reshape(nanstd(temp(:,1,:),[],3)./sqrt(size(dp,3)),size(dp,1),1),'Color','m');
    plot(timebinMeans,reshape(nanmean(temp(:,1,:),3),size(dp,1),1)-reshape(nanstd(temp(:,1,:),[],3)./sqrt(size(dp,3)),size(dp,1),1),'Color','m');
    legend({'ta true p active SECOND TRIAL TYPE','ta false p active SECOND TRIAL TYPE'});

    dpwithinsess=norminv(p_success_ubyu(:,1,:))-norminv(p_success_ubyu(:,2,:));
    dpwithinsess(isinf(dpwithinsess))=nan;
    temp=dpwithinsess;
    figure(); plot(timebinMeans,reshape(nanmean(temp(:,1,:),3),size(dp,1),1),'Color','k'); hold on;
    plot(timebinMeans,reshape(nanmean(temp(:,1,:),3),size(dp,1),1)+reshape(nanstd(temp(:,1,:),[],3)./sqrt(size(dp,3)),size(dp,1),1),'Color','k');
    plot(timebinMeans,reshape(nanmean(temp(:,1,:),3),size(dp,1),1)-reshape(nanstd(temp(:,1,:),[],3)./sqrt(size(dp,3)),size(dp,1),1),'Color','k');
    legend({'dprime prob ta true vs ta false FIRST TRIAL TYPE'});

    backup.fr_succ_sd=fr_succ_sd; backup.fr_fail_sd=fr_fail_sd;
    fr_succ_sd(fr_succ_sd<0.001)=0; fr_fail_sd(fr_fail_sd<0.001)=0;

    dpwithinsess=(fr_success_ubyu(:,1,:)-fr_success_ubyu(:,2,:))./sqrt(fr_succ_sd(:,1,:).^2+fr_succ_sd(:,2,:).^2);
    dpwithinsess(isinf(dpwithinsess))=nan;
    temp=dpwithinsess;
    figure(); plot(timebinMeans,reshape(nanmean(temp(:,1,:),3),size(dp,1),1),'Color','k'); hold on;
    plot(timebinMeans,reshape(nanmean(temp(:,1,:),3),size(dp,1),1)+reshape(nanstd(temp(:,1,:),[],3)./sqrt(size(dp,3)),size(dp,1),1),'Color','k');
    plot(timebinMeans,reshape(nanmean(temp(:,1,:),3),size(dp,1),1)-reshape(nanstd(temp(:,1,:),[],3)./sqrt(size(dp,3)),size(dp,1),1),'Color','k');
    legend({'dprime from fr ta true vs ta false FIRST TRIAL TYPE'});

    nneuronsinthissess=nan(length(usess),1);
    for i=1:length(usess)
    currsess=usess(i);
    nneuronsinthissess(i)=length(unique(backup.success_Response.fromWhichUnit(success_Response.fromWhichUnit==1 & success_Response.fromWhichSess_forTrials==currsess)));
    end
    figure(); histogram(nneuronsinthissess,200);
    unneur=unique(nneuronsinthissess); unneur=unneur(unneur~=0);
    figure();
    cs={'m','k','g','c','y','r','b'};
    dps=cell(1,length(unneur));
    for i=1:length(unneur)
        hold on;
        dpwithinsess=(fr_success_ubyu(:,1,nneuronsinthissess==unneur(i))-fr_success_ubyu(:,2,nneuronsinthissess==unneur(i)))./sqrt(fr_succ_sd(:,1,nneuronsinthissess==unneur(i)).^2+fr_succ_sd(:,2,nneuronsinthissess==unneur(i)).^2);
        dpwithinsess(isinf(dpwithinsess))=nan;
        temp=dpwithinsess;
        plot(timebinMeans,reshape(nanmean(temp(:,1,:),3),size(dp,1),1),'Color',cs{i}); hold on;
        plot(timebinMeans,reshape(nanmean(temp(:,1,:),3),size(dp,1),1)+reshape(nanstd(temp(:,1,:),[],3)./sqrt(size(dp,3)),size(dp,1),1),'Color',cs{i});
        plot(timebinMeans,reshape(nanmean(temp(:,1,:),3),size(dp,1),1)-reshape(nanstd(temp(:,1,:),[],3)./sqrt(size(dp,3)),size(dp,1),1),'Color',cs{i});
        dps{i}=nanmean(temp(timebinMeans>2 & timebinMeans<=4,1,:),1);
    end
    figure(); 
    for i=1:length(unneur)
        temp=dps{i}; temp=reshape(temp,length(temp),1); dps{i}=temp;
        scatter(ones(size(temp)).*i,temp); hold on;
    end
    figure();
    violin(dps); %,'bw',0.4);
end

end

function [dp_all,dpfr_all,isSig_all,pval_all,p_success_unitbyunit_all,p_failure_unitbyunit_all,fr_success_unitbyunit_all,fr_failure_unitbyunit_all,timebinMeans,frsd_success_unitbyunit_all,frsd_failure_unitbyunit_all]=withinSessDprime(success_Response,failure_Response,timeBins)

for i=1:size(timeBins,1)
    [dp,dpfr,isSig,pval,p_success_unitbyunit,p_failure_unitbyunit,fr_success_unitbyunit,fr_failure_unitbyunit,frsd_success_unitbyunit,frsd_failure_unitbyunit]=dprimeInTimeSessBySess(success_Response,failure_Response,timeBins(i,:));
    if i==1
        dp_all=nan(size(timeBins,1),size(dp,1),size(dp,2));
        dpfr_all=nan(size(timeBins,1),size(dp,1),size(dp,2));
        isSig_all=nan(size(timeBins,1),size(dp,1),size(dp,2));
        pval_all=nan(size(timeBins,1),size(dp,1),size(dp,2));
        p_success_unitbyunit_all=nan(size(timeBins,1),size(dp,1),size(dp,2));
        p_failure_unitbyunit_all=nan(size(timeBins,1),size(dp,1),size(dp,2));
        fr_success_unitbyunit_all=nan(size(timeBins,1),size(dp,1),size(dp,2));
        fr_failure_unitbyunit_all=nan(size(timeBins,1),size(dp,1),size(dp,2));
        frsd_success_unitbyunit_all=nan(size(timeBins,1),size(dp,1),size(dp,2));
        frsd_failure_unitbyunit_all=nan(size(timeBins,1),size(dp,1),size(dp,2));
        timebinMeans=nan(size(timeBins,1),1);
    end
    dp_all(i,:,:)=dp;
    dpfr_all(i,:,:)=dpfr;
    if ~isempty(isSig)
        isSig_all(i,:,:)=isSig;
        pval_all(i,:,:)=pval;
    end
    p_success_unitbyunit_all(i,:,:)=p_success_unitbyunit;
    p_failure_unitbyunit_all(i,:,:)=p_failure_unitbyunit;
    fr_success_unitbyunit_all(i,:,:)=fr_success_unitbyunit;
    fr_failure_unitbyunit_all(i,:,:)=fr_failure_unitbyunit;
    frsd_success_unitbyunit_all(i,:,:)=frsd_success_unitbyunit;
    frsd_failure_unitbyunit_all(i,:,:)=frsd_failure_unitbyunit;
    timebinMeans(i)=mean(timeBins(i,:));
end

end

function [dp_all,dpfr_all,isSig_all,pval_all,p_success_unitbyunit_all,p_failure_unitbyunit_all,fr_success_unitbyunit_all,fr_failure_unitbyunit_all,timebinMeans]=dprimesOverTime(success_Response,failure_Response,timeBins)

for i=1:size(timeBins,1)
    [dp,dpfr,isSig,pval,p_success_unitbyunit,p_failure_unitbyunit,fr_success_unitbyunit,fr_failure_unitbyunit]=dprimeInTime(success_Response,failure_Response,timeBins(i,:));
    if i==1
        dp_all=nan(size(timeBins,1),length(dp));
        dpfr_all=nan(size(timeBins,1),length(dp));
        isSig_all=nan(size(timeBins,1),length(dp));
        pval_all=nan(size(timeBins,1),length(dp));
        p_success_unitbyunit_all=nan(size(timeBins,1),length(dp));
        p_failure_unitbyunit_all=nan(size(timeBins,1),length(dp));
        fr_success_unitbyunit_all=nan(size(timeBins,1),length(dp));
        fr_failure_unitbyunit_all=nan(size(timeBins,1),length(dp));
        timebinMeans=nan(size(timeBins,1),1);
    end
    dp_all(i,:)=dp;
    dpfr_all(i,:)=dpfr;
    if ~isempty(isSig)
        isSig_all(i,:)=isSig;
        pval_all(i,:)=pval;
    end
    p_success_unitbyunit_all(i,:)=p_success_unitbyunit;
    p_failure_unitbyunit_all(i,:)=p_failure_unitbyunit;
    fr_success_unitbyunit_all(i,:)=fr_success_unitbyunit;
    fr_failure_unitbyunit_all(i,:)=fr_failure_unitbyunit;
    timebinMeans(i)=mean(timeBins(i,:));
end

end

function [dp,dpfr,isSig,pval,p_success_unitbyunit,p_failure_unitbyunit,fr_success_unitbyunit,fr_failure_unitbyunit,frsd_success_unitbyunit,frsd_failure_unitbyunit]=dprimeInTimeSessBySess(success_Response,failure_Response,timeWindow)

isSig=[];
pval=[];
nBoots=100;

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
fromWhichSess_success=success_Response.fromWhichSess_forTrials;
unitfr_success=unitfr_success(success_Response.isEventInThisTrial==1);
fromWhichUnit_success=fromWhichUnit_success(success_Response.isEventInThisTrial==1);
fromWhichSess_success=fromWhichSess_success(success_Response.isEventInThisTrial==1);
unitfr_failure=sum(failure_Response.unitbyunit_y(:,failureRange(1):failureRange(2)),2,'omitnan');
fromWhichUnit_failure=failure_Response.fromWhichUnit;
fromWhichSess_failure=failure_Response.fromWhichSess_forTrials;
unitfr_failure=unitfr_failure(failure_Response.isEventInThisTrial==1);
fromWhichUnit_failure=fromWhichUnit_failure(failure_Response.isEventInThisTrial==1);
fromWhichSess_failure=fromWhichSess_failure(failure_Response.isEventInThisTrial==1);

units=unique(success_Response.fromWhichUnit);
% [p_success_unitbyunit,p_failure_unitbyunit]=getProbOfResponse(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure);
[p_success_unitbyunit,p_failure_unitbyunit]=getProbOfResponseSessBySess(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure,fromWhichSess_success,fromWhichSess_failure);
% [isSig,pval]=getSignificantUnits(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure,nBoots);
% maxdp=3;
dp=norminv(p_success_unitbyunit)-norminv(p_failure_unitbyunit);
% dp(dp<-maxdp)=-maxdp; dp(dp>maxdp)=maxdp;

[fr_success_unitbyunit,fr_failure_unitbyunit,frsd_success_unitbyunit,frsd_failure_unitbyunit]=getFROfResponseSessBySess(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure,fromWhichSess_success,fromWhichSess_failure);
dpfr=dprime_from_FR(fr_success_unitbyunit,fr_failure_unitbyunit,frsd_success_unitbyunit,frsd_failure_unitbyunit);

end

function [dp,dpfr,isSig,pval,p_success_unitbyunit,p_failure_unitbyunit,fr_success_unitbyunit,fr_failure_unitbyunit]=dprimeInTime(success_Response,failure_Response,timeWindow)

isSig=[];
pval=[];
nBoots=100;

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

% [isSig,pval]=getSignificantUnits(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure,nBoots);
% maxdp=3;
dp=norminv(p_success_unitbyunit)-norminv(p_failure_unitbyunit);
% dp(dp<-maxdp)=-maxdp; dp(dp>maxdp)=maxdp;

[fr_success_unitbyunit,fr_failure_unitbyunit,frsd_success_unitbyunit,frsd_failure_unitbyunit]=getFROfResponse(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure);
dpfr=dprime_from_FR(fr_success_unitbyunit,fr_failure_unitbyunit,frsd_success_unitbyunit,frsd_failure_unitbyunit);

end

function dp=dprime_from_FR(data1_me,data2_me,data1_sd,data2_sd)

% RMS sd discriminability index
% suboptimal but simple
dp=(data1_me-data2_me)./sqrt(data1_sd.^2+data2_sd.^2);

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

function [fr_success_unitbyunit,fr_failure_unitbyunit,frsd_success_unitbyunit,frsd_failure_unitbyunit]=getFROfResponse(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure)

fr_success_unitbyunit=nan(length(units),1);
fr_failure_unitbyunit=nan(length(units),1);
frsd_success_unitbyunit=nan(length(units),1);
frsd_failure_unitbyunit=nan(length(units),1);
disp([num2str(length(units)) ' units']);
for i=1:length(units)
    fr_success_unitbyunit(i)=nanmean(unitfr_success(fromWhichUnit_success==units(i)));
    fr_failure_unitbyunit(i)=nanmean(unitfr_failure(fromWhichUnit_failure==units(i)));
    frsd_success_unitbyunit(i)=nanstd(unitfr_success(fromWhichUnit_success==units(i)));
    frsd_failure_unitbyunit(i)=nanstd(unitfr_failure(fromWhichUnit_failure==units(i)));
end

end

function [fr_success_unitbyunit,fr_failure_unitbyunit,frsd_success_unitbyunit,frsd_failure_unitbyunit]=getFROfResponseSessBySess(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure,fromWhichSess_success,fromWhichSess_failure)

uSess=unique([fromWhichSess_success; fromWhichSess_failure]);
fr_success_unitbyunit=nan(length(units),1);
fr_failure_unitbyunit=nan(length(units),1);
frsd_success_unitbyunit=nan(length(units),1);
frsd_failure_unitbyunit=nan(length(units),1);
for j=1:length(uSess)
    currSess=uSess(j);
    for i=1:length(units)
        fr_success_unitbyunit(i,j)=nanmean(unitfr_success(fromWhichUnit_success==units(i) & fromWhichSess_success==currSess));
        fr_failure_unitbyunit(i,j)=nanmean(unitfr_failure(fromWhichUnit_failure==units(i) & fromWhichSess_failure==currSess));
        frsd_success_unitbyunit(i,j)=nanstd(unitfr_success(fromWhichUnit_success==units(i) & fromWhichSess_success==currSess));
        frsd_failure_unitbyunit(i,j)=nanstd(unitfr_failure(fromWhichUnit_failure==units(i) & fromWhichSess_failure==currSess));
    end
end

end

function [p_success_unitbyunit,p_failure_unitbyunit]=getProbOfResponseSessBySess(units,unitfr_success,unitfr_failure,fromWhichUnit_success,fromWhichUnit_failure,fromWhichSess_success,fromWhichSess_failure)

uSess=unique([fromWhichSess_success; fromWhichSess_failure]);
p_success_unitbyunit=nan(length(units),length(uSess));
p_failure_unitbyunit=nan(length(units),length(uSess));
for j=1:length(uSess)
    currSess=uSess(j);
    for i=1:length(units)
        p_success_unitbyunit(i,j)=nansum(unitfr_success(fromWhichUnit_success==units(i) & fromWhichSess_success==currSess)>0.5)./nansum(fromWhichUnit_success==units(i) & fromWhichSess_success==currSess);
        p_failure_unitbyunit(i,j)=nansum(unitfr_failure(fromWhichUnit_failure==units(i) & fromWhichSess_failure==currSess)>0.5)./nansum(fromWhichUnit_failure==units(i) & fromWhichSess_failure==currSess);
    end
end

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



