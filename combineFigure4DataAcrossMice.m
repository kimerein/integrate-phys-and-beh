function combineFigure4DataAcrossMice(dd,useOptoForThisGp,whichtoplot,timeStep,cueind,doCDF,doLEDcdf)

rr_noLED_tri1_uncued=[];
rr_noLED_tri1_cued=[];
rr_noLED_trinext_uncued=[];
rr_noLED_trinext_cued=[];

rr_LED_tri1_uncued=[];
rr_LED_tri1_cued=[];
rr_LED_trinext_uncued=[];
rr_LED_trinext_cued=[];

returnThis_control=cell(length(dd));
returnThis_LED=cell(length(dd));

forcdf_rawreach_trial1=[];
forcdf_rawreach_trial2=[];
for i=1:length(dd)
    currdir=dd{i};
    if doCDF==false
        a=load(fullfile(currdir,whichtoplot,'pdfcdfcontrol','returnThis.mat'));
        returnThis_control{i}=a.returnThis;
        a=load(fullfile(currdir,whichtoplot,'pdfcdfcontrol','returnThis.mat'));
        returnThis_control{i}=a.returnThis;

        a=load(fullfile(currdir,whichtoplot,'pdfcdfpDMStinh','returnThis.mat'));
        returnThis_LED{i}=a.returnThis;
        a=load(fullfile(currdir,whichtoplot,'pdfcdfpDMStinh','returnThis.mat'));
        returnThis_LED{i}=a.returnThis;

        a=load(fullfile(currdir,whichtoplot,'returnout.mat'));
        rout=a.returnout;
        if isempty(rout.reachrates_noLED)
            % nothing to add
            continue
        end
        temp=rout.reachrates_noLED.trial1_alltrials_uncued; temp=temp(:).';
        rr_noLED_tri1_uncued=[rr_noLED_tri1_uncued temp(1:end)];
        temp=rout.reachrates_noLED.trial1_alltrials_cued; temp=temp(:).';
        rr_noLED_tri1_cued=[rr_noLED_tri1_cued temp(1:end)];
        temp=rout.reachrates_noLED.alltrials_uncued; temp=temp(:).';
        rr_noLED_trinext_uncued=[rr_noLED_trinext_uncued temp(1:end)];
        temp=rout.reachrates_noLED.alltrials_cued; temp=temp(:).';
        rr_noLED_trinext_cued=[rr_noLED_trinext_cued temp(1:end)];
        % dump nan, which are just padding
        ina=isnan(rr_noLED_tri1_uncued);
        rr_noLED_tri1_uncued=rr_noLED_tri1_uncued(~ina);
        rr_noLED_tri1_cued=rr_noLED_tri1_cued(~ina);
        rr_noLED_trinext_uncued=rr_noLED_trinext_uncued(~ina);
        rr_noLED_trinext_cued=rr_noLED_trinext_cued(~ina);

        if isempty(rout.reachrates_LED) || useOptoForThisGp(i)==0
            % nothing to add
            continue
        end
        temp=rout.reachrates_LED.trial1_alltrials_uncued; temp=temp(:).';
        rr_LED_tri1_uncued=[rr_LED_tri1_uncued temp(1:end)];
        temp=rout.reachrates_LED.trial1_alltrials_cued; temp=temp(:).';
        rr_LED_tri1_cued=[rr_LED_tri1_cued temp(1:end)];
        temp=rout.reachrates_LED.alltrials_uncued; temp=temp(:).';
        rr_LED_trinext_uncued=[rr_LED_trinext_uncued temp(1:end)];
        temp=rout.reachrates_LED.alltrials_cued; temp=temp(:).';
        rr_LED_trinext_cued=[rr_LED_trinext_cued temp(1:end)];
        % dump nan, which are just padding
        ina=isnan(rr_LED_tri1_uncued);
        rr_LED_tri1_uncued=rr_LED_tri1_uncued(~ina);
        rr_LED_tri1_cued=rr_LED_tri1_cued(~ina);
        rr_LED_trinext_uncued=rr_LED_trinext_uncued(~ina);
        rr_LED_trinext_cued=rr_LED_trinext_cued(~ina);
    else
        % do CDF only
        if doLEDcdf==false
            a=load(fullfile(currdir,whichtoplot,'pdfcdfpDMStinh','returnThisCDF.mat'));
            rthisCDF=a.returnThisCDF;
            forcdf_rawreach_trial1=[forcdf_rawreach_trial1; rthisCDF.trial1_rawReachMatrix];
            forcdf_rawreach_trial2=[forcdf_rawreach_trial2; rthisCDF.trial2_rawReachMatrix];
        else
            a=load(fullfile(currdir,whichtoplot,'pdfcdfcontrol','returnThisCDF.mat'));
            rthisCDF=a.returnThisCDF;
            forcdf_rawreach_trial1=[forcdf_rawreach_trial1; rthisCDF.trial1_rawReachMatrix];
            forcdf_rawreach_trial2=[forcdf_rawreach_trial2; rthisCDF.trial2_rawReachMatrix];
        end

        returnThisCDF.trial1_rawReachMatrix=forcdf_rawreach_trial1;
        returnThisCDF.trial2_rawReachMatrix=forcdf_rawreach_trial2;
        getMeanAndBootstrapForCDF(timeStep,returnThisCDF,cueind);
    end
end

% do bootstrapped scatter
figure();
[uncued_mean_out,cued_mean_out,bootMeans]=bootstrap(rr_noLED_tri1_uncued,rr_noLED_tri1_cued,'k','k',false);
plotMeAndSe(rr_noLED_tri1_uncued,rr_noLED_tri1_cued,'k',2,false);
[uncued_mean_out,cued_mean_out,bootMeans]=bootstrap(rr_noLED_trinext_uncued,rr_noLED_trinext_cued,[0.5 0.5 0.5],[0.5 0.5 0.5],false);
plotMeAndSe(rr_noLED_trinext_uncued,rr_noLED_trinext_cued,[0.5 0.5 0.5],2,false);

figure();
[uncued_mean_out,cued_mean_out,bootMeans]=bootstrap(rr_LED_tri1_uncued,rr_LED_tri1_cued,'k','r',false);
plotMeAndSe(rr_LED_tri1_uncued,rr_LED_tri1_cued,'k',2,false);
[uncued_mean_out,cued_mean_out,bootMeans]=bootstrap(rr_LED_trinext_uncued,rr_LED_trinext_cued,'r','r',false);
plotMeAndSe(rr_LED_trinext_uncued,rr_LED_trinext_cued,'r',2,false);


end

function plotMeAndSe(data1,data2,c,linewidth,suppressOutput)
% make inputs vectors if they are not
data1=data1(1:end);
data2=data2(1:end);
% average within each session
% data1=nanmean(data1,2); data1=data1';
% data2=nanmean(data2,2); data2=data2';
if suppressOutput==false
    line([nanmean(data1)-nanstd(data1,[],2)./sqrt(nansum(~isnan(data1))) nanmean(data1)+nanstd(data1,[],2)./sqrt(nansum(~isnan(data1)))],[nanmean(data2) nanmean(data2)],'Color',c,'LineWidth',linewidth);
    hold on;
    line([nanmean(data1) nanmean(data1)],[nanmean(data2)-nanstd(data2,[],2)./sqrt(nansum(~isnan(data2))) nanmean(data2)+nanstd(data2,[],2)./sqrt(nansum(~isnan(data2)))],'Color',c,'LineWidth',linewidth);
end

end

function getMeanAndBootstrapForCDF(timeStep,returnThisCDF,f)

startAtCue=false;
subtractPreCue=false;
startAtPrecue=true;
preCueWindow=[-2 -1];
cutCDFat=9; % cut cdf at this time

maxTrile=cutCDFat; % make this empty if want whole trial length

% timeStep=mode(diff(nanmean(alltbt.times,1)));
timeBinsForReaching=0:timeStep:(size(returnThisCDF.trial1_rawReachMatrix,2)-1)*timeStep;
% find cue ind
% [~,f]=nanmax(nanmean(alltbt.cueZone_onVoff,1));
if startAtPrecue==true
    cueTime=timeBinsForReaching(f);
    precuecutoff=cueTime+preCueWindow(1);

    temp=returnThisCDF.trial1_rawReachMatrix;
    temp(:,timeBinsForReaching<precuecutoff)=0;
    returnThisCDF.trial1_rawReachMatrix=temp;
    
    temp=returnThisCDF.trial2_rawReachMatrix;
    temp(:,timeBinsForReaching<precuecutoff)=0;
    returnThisCDF.trial2_rawReachMatrix=temp;
end
if startAtCue==true
    cueTime=timeBinsForReaching(f);
    preCueWindow=[cueTime+preCueWindow(1) cueTime+preCueWindow(2)];
else
    % make preCueWindow wrt cueTime=0
    cueTime=timeBinsForReaching(f);
    preCueWindow=[cueTime+preCueWindow(1) cueTime+preCueWindow(2)];
    cueTime=0;
end
[d1,d2]=plotCDF_rawReaches(returnThisCDF.trial1_rawReachMatrix,returnThisCDF.trial2_rawReachMatrix,timeBinsForReaching,cueTime,['CDF Raw Reaches: trial 1 (black) vs later (red)'],subtractPreCue,preCueWindow,maxTrile);

% Bootstrap CDFs
dat2=returnThisCDF.trial1_rawReachMatrix;
takeFracForBootstrap=0.66;
takeIndsForBootstrap=ceil(takeFracForBootstrap*size(dat2,1));
nRuns=100;
bootCDFs=nan(nRuns,length(timeBinsForReaching));
for i=1:nRuns
    takeTheseForBoot=randi(size(dat2,1),1,takeIndsForBootstrap); % with replacement
    sub_dat2=dat2(takeTheseForBoot,:);
    sub_dat2=sum(sub_dat2,1,'omitnan');
    if isempty(maxTrile)
    else
        sub_dat2(timeBinsForReaching>maxTrile)=0;
    end
    cond2_cdf=accumulateDistribution(sub_dat2);
    cond2_cdf=cond2_cdf./nanmax(cond2_cdf);
    bootCDFs(i,:)=cond2_cdf;
end
% Show bootstrapped 95% CI
sorted_bootCDFs=nan(size(bootCDFs));
fifthPerc=nan(1,size(bootCDFs,2));
ninetyfifthPerc=nan(1,size(bootCDFs,2));
for i=1:size(bootCDFs,2)
    sorted_bootCDFs(:,i)=sort(bootCDFs(:,i));
    fifthPerc(i)=prctile(sorted_bootCDFs(:,i),5);
    ninetyfifthPerc(i)=prctile(sorted_bootCDFs(:,i),95);
end
plot(timeBinsForReaching,fifthPerc,'Color','k'); hold on;
plot(timeBinsForReaching,ninetyfifthPerc,'Color','k');

dat2=returnThisCDF.trial2_rawReachMatrix;
takeFracForBootstrap=0.66;
takeIndsForBootstrap=ceil(takeFracForBootstrap*size(dat2,1));
nRuns=100;
bootCDFs=nan(nRuns,length(timeBinsForReaching));
for i=1:nRuns
    takeTheseForBoot=randi(size(dat2,1),1,takeIndsForBootstrap); % with replacement
    sub_dat2=dat2(takeTheseForBoot,:);
    sub_dat2=sum(sub_dat2,1,'omitnan');
    if isempty(maxTrile)
    else
        sub_dat2(timeBinsForReaching>maxTrile)=0;
    end
    cond2_cdf=accumulateDistribution(sub_dat2);
    cond2_cdf=cond2_cdf./nanmax(cond2_cdf);
    bootCDFs(i,:)=cond2_cdf;
end
% Show bootstrapped 95% CI
sorted_bootCDFs=nan(size(bootCDFs));
fifthPerc=nan(1,size(bootCDFs,2));
ninetyfifthPerc=nan(1,size(bootCDFs,2));
for i=1:size(bootCDFs,2)
    sorted_bootCDFs(:,i)=sort(bootCDFs(:,i));
    fifthPerc(i)=prctile(sorted_bootCDFs(:,i),5);
    ninetyfifthPerc(i)=prctile(sorted_bootCDFs(:,i),95);
end
plot(timeBinsForReaching,fifthPerc,'Color','r'); hold on;
plot(timeBinsForReaching,ninetyfifthPerc,'Color','r');

end


function [uncued_mean_out,cued_mean_out,bootMeans]=bootstrap(approach2_alltrials_uncued,approach2_alltrials_cued,colorForBootstrapPoints,scatterPointsEdgeColor,suppressOutput)

stillsuppressbootstrap=false;

% altogether_prob_cued=nanmean(approach2_alltrials_cued,2);
% altogether_prob_uncued=nanmean(approach2_alltrials_uncued,2);
altogether_prob_cued=approach2_alltrials_cued(1:end); % better to bootstrap across trials, not sessions, because in some sessions, mouse drops a lot
altogether_prob_uncued=approach2_alltrials_uncued(1:end);
takeTrials=~isnan(altogether_prob_cued) & ~isnan(altogether_prob_uncued);
% nan trials are from cases where a session had fewer trials, just nanned
% to fill in matrix
% disp(['dropping this many trials because of nan ' num2str(nansum(~takeTrials))]);
altogether_prob_cued=altogether_prob_cued(takeTrials==1);
altogether_prob_uncued=altogether_prob_uncued(takeTrials==1);
% Show bootstrapped 95% CI
takeFracForBootstrap=0.66;
takeIndsForBootstrap=ceil(takeFracForBootstrap*length(altogether_prob_cued));
nRuns=100;
bootMeans=nan(2,nRuns);
for i=1:nRuns
    takeTheseForBoot=randi(length(altogether_prob_cued),1,takeIndsForBootstrap); % with replacement
    sub_prob_cued=altogether_prob_cued(takeTheseForBoot);
    sub_prob_uncued=altogether_prob_uncued(takeTheseForBoot);
    bootMeans(1,i)=nanmean(sub_prob_uncued);
    bootMeans(2,i)=nanmean(sub_prob_cued);
end
if suppressOutput==false && stillsuppressbootstrap==false
    s=scatter(bootMeans(1,:),bootMeans(2,:),20,'filled','MarkerEdgeColor',scatterPointsEdgeColor,'MarkerFaceColor',colorForBootstrapPoints); hold on;
    s.AlphaData = 0.5*ones(1,size(bootMeans,2));
    s.MarkerFaceAlpha = 'flat';
    scatter(nanmean(altogether_prob_uncued),nanmean(altogether_prob_cued),50,'filled','MarkerEdgeColor',scatterPointsEdgeColor,'MarkerFaceColor',colorForBootstrapPoints);
end
uncued_mean_out=nanmean(altogether_prob_uncued);
cued_mean_out=nanmean(altogether_prob_cued);

end

function [out1,out2]=plotCDF_rawReaches(data1,data2,timesteps,cueTime,tit,subtractPreCue,preCueWindow,maxTrialLength)

% make raw reaching data a timeseries (i.e., a pdf)
% then simply accumulate distribution, selecting only time points after the
% cue

% maxTrialLength=[]; %9.5; % in sec, or empty if want to include all reaches

data1=sum(data1,1,'omitnan');
data2=sum(data2,1,'omitnan');

if subtractPreCue==true
    data1=data1-mean(data1(timesteps>=preCueWindow(1) & timesteps<=preCueWindow(2)),2,'omitnan');
    data2=data2-mean(data2(timesteps>=preCueWindow(1) & timesteps<=preCueWindow(2)),2,'omitnan');
    data1(data1<0)=0;
    data2(data2<0)=0;
end
if ~isempty(maxTrialLength)
    data1(timesteps>maxTrialLength)=0;
    data2(timesteps>maxTrialLength)=0;
end

suppPlots=whetherToSuppressPlots();

if suppPlots==false
    cond1_cdf=accumulateDistribution(data1(timesteps>cueTime));
    figure();
    plot(timesteps(timesteps>cueTime),cond1_cdf./nanmax(cond1_cdf),'Color','k');
    xlabel('CDF');
    ylabel('Count');
    title(tit);
    out1=cond1_cdf./nanmax(cond1_cdf);
    
    hold on;
    cond2_cdf=accumulateDistribution(data2(timesteps>cueTime));
    plot(timesteps(timesteps>cueTime),cond2_cdf./nanmax(cond2_cdf),'Color','r');
    out2=cond2_cdf./nanmax(cond2_cdf);
end

end

function cdf=accumulateDistribution(data)

cdf=nan(size(data));
for i=1:length(data)
    cdf(i)=nansum(data(1:i));
end

end