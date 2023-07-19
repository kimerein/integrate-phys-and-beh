function combineFigure4DataAcrossMice(dd,whichtoplot,timeStep,cueind,doCDF,doLEDcdf)

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
        rout=load(fullfile(currdir,whichtoplot,'returnout.mat'));
        temp=rout.reachrates_noLED.trial1_alltrials_uncued;
        rr_noLED_tri1_uncued=[rr_noLED_tri1_uncued temp(1:end)];
        temp=rout.reachrates_noLED.trial1_alltrials_cued;
        rr_noLED_tri1_cued=[rr_noLED_tri1_cued temp(1:end)];
        temp=rout.reachrates_noLED.alltrials_uncued;
        rr_noLED_trinext_uncued=[rr_noLED_trinext_uncued temp(1:end)];
        temp=rout.reachrates_noLED.alltrials_cued;
        rr_noLED_trinext_cued=[rr_noLED_trinext_cued temp(1:end)];
        % dump nan, which are just padding
        ina=isnan(rr_noLED_tri1_uncued) | isnan(rr_noLED_tri1_cued) | isnan(rr_noLED_trinext_uncued) | isnan(rr_noLED_trinext_cued);
        rr_noLED_tri1_uncued=rr_noLED_tri1_uncued(~ina);
        rr_noLED_tri1_cued=rr_noLED_tri1_cued(~ina);
        rr_noLED_trinext_uncued=rr_noLED_trinext_uncued(~ina);
        rr_noLED_trinext_cued=rr_noLED_trinext_cued(~ina);

        temp=rout.reachrates_LED.trial1_alltrials_uncued;
        rr_LED_tri1_uncued=[rr_LED_tri1_uncued temp(1:end)];
        temp=rout.reachrates_LED.trial1_alltrials_cued;
        rr_LED_tri1_cued=[rr_LED_tri1_cued temp(1:end)];
        temp=rout.reachrates_LED.alltrials_uncued;
        rr_LED_trinext_uncued=[rr_LED_trinext_uncued temp(1:end)];
        temp=rout.reachrates_LED.alltrials_cued;
        rr_LED_trinext_cued=[rr_LED_trinext_cued temp(1:end)];
        % dump nan, which are just padding
        ina=isnan(rr_LED_tri1_uncued) | isnan(rr_LED_tri1_cued) | isnan(rr_LED_trinext_uncued) | isnan(rr_LED_trinext_cued);
        rr_LED_tri1_uncued=rr_LED_tri1_uncued(~ina);
        rr_LED_tri1_cued=rr_LED_tri1_cued(~ina);
        rr_LED_trinext_uncued=rr_LED_trinext_uncued(~ina);
        rr_LED_trinext_cued=rr_LED_trinext_cued(~ina);

        rthis=load(fullfile(currdir,whichtoplot,'pdfcdfcontrol','returnThis.mat'));
        returnThis_control{i}=rthis;
        rthis=load(fullfile(currdir,whichtoplot,'pdfcdfcontrol','returnThis.mat'));
        returnThis_control{i}=rthis;

        rthis=load(fullfile(currdir,whichtoplot,'pdfcdfpDMStinh','returnThis.mat'));
        returnThis_LED{i}=rthis;
        rthis=load(fullfile(currdir,whichtoplot,'pdfcdfpDMStinh','returnThis.mat'));
        returnThis_LED{i}=rthis;
    else
        % do CDF only
        if doLEDcdf==false
            rthisCDF=load(fullfile(currdir,whichtoplot,'pdfcdfpDMStinh','returnThisCDF.mat'));
            forcdf_rawreach_trial1=
            forcdf_rawreach_trial2
        else
            rthisCDF=load(fullfile(currdir,whichtoplot,'pdfcdfcontrol','returnThisCDF.mat'));
        end
    end
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
dat2=dataset.rawReaching_event_trial1InSeq{1};
returnThis.trial1_rawReachMatrix=dat2;
% dat2=dataset.rawReaching_event_trialiInSeq{1};
takeFracForBootstrap=0.66;
takeIndsForBootstrap=ceil(takeFracForBootstrap*size(dat2,1));
nRuns=100;
bootCDFs=nan(nRuns,length(timeBinsForReaching));
if size(dat2,1)==0
    returnThis.trial2_rawReachMatrix=dataset.rawReaching_event_trialiInSeq{1};
    return
end
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
% dat1=dataset.rawReaching_event_trial1InSeq{1};
dat2=dataset.rawReaching_event_trialiInSeq{1};
returnThis.trial2_rawReachMatrix=dat2;
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

if startAtCue==false
    line([timeBinsForReaching(f) timeBinsForReaching(f)],[0 1],'Color','b');
    c=polyfit(timeBinsForReaching(timeBinsForReaching>=preCueWindow(1) & timeBinsForReaching<=preCueWindow(2)),d2(timeBinsForReaching>=preCueWindow(1) & timeBinsForReaching<=preCueWindow(2)),1);
    yest=polyval(c,timeBinsForReaching);
    plot(timeBinsForReaching,yest,'Color','g');
    ylim([0 1]);
    returnThis.uncued_fit=c;
end
returnThis.timeBins=timeBinsForReaching;
returnThis.cueTime=timeBinsForReaching(f);
returnThis.cdf_trial1=d1;
returnThis.cdf_trial2=d2;

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