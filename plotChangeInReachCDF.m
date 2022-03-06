function plotChangeInReachCDF(dataset,alltbt)

preCueWindow=[-2 -1];
maxTrialLength=9.5;
maxTrialLength_for_CDFs=[];

timeStep=mode(diff(nanmean(alltbt.times,1)));
timeBinsForReaching=0:timeStep:(size(dataset.rawReaching_allTrialsSequence_trial1InSeq{1},2)-1)*timeStep;
% find cue ind
[~,f]=nanmax(nanmean(alltbt.cueZone_onVoff,1));
% make preCueWindow wrt cueTime=0
cueTime=timeBinsForReaching(f);
preCueWindow=[cueTime+preCueWindow(1) cueTime+preCueWindow(2)];
cueTime=0;

for i=1:length(dataset.allTrialsSequence_RT_trial1InSeq)
    [d1,d2]=plotCDF_rawReaches(dataset.rawReaching_allTrialsSequence_trial1InSeq{i},dataset.rawReaching_allTrialsSequence_trialiInSeq{i},timeBinsForReaching,cueTime,['CDF Raw Reaches all trials reference: trial 1 (black) vs ' num2str(dataset.nInSequence(i)-1) ' later (red)'],false,preCueWindow,maxTrialLength_for_CDFs);
end

line([timeBinsForReaching(f) timeBinsForReaching(f)],[0 1],'Color','b');
c=polyfit(timeBinsForReaching(timeBinsForReaching>=preCueWindow(1) & timeBinsForReaching<=preCueWindow(2)),d1(timeBinsForReaching>=preCueWindow(1) & timeBinsForReaching<=preCueWindow(2)),1);
yest=polyval(c,timeBinsForReaching);
plot(timeBinsForReaching,yest,'Color','k');
ylim([0 1]);
xlim([0 maxTrialLength]);

temp=dataset.event_name;
temp(regexp(temp,'_'))=' ';
for i=1:length(dataset.event_RT_trial1InSeq)
    [d1,d2]=plotCDF_rawReaches(dataset.rawReaching_event_trial1InSeq{i},dataset.rawReaching_event_trialiInSeq{i},timeBinsForReaching,cueTime,['CDF Raw Reaches fx of ' temp ': trial 1 (black) vs ' num2str(dataset.nInSequence(i)-1) ' later (red)'],false,preCueWindow,maxTrialLength_for_CDFs);
end
% Bootstrap CDFs
% dat1=dataset.rawReaching_event_trial1InSeq{1};
dat2=dataset.rawReaching_event_trialiInSeq{1};
takeFracForBootstrap=0.66;
takeIndsForBootstrap=ceil(takeFracForBootstrap*size(dat2,1));
nRuns=100;
bootCDFs=nan(nRuns,length(timeBinsForReaching));
for i=1:nRuns
    takeTheseForBoot=randi(size(dat2,1),1,takeIndsForBootstrap); % with replacement
    sub_dat2=dat2(takeTheseForBoot,:);
    sub_dat2=sum(sub_dat2,1,'omitnan');
    if isempty(maxTrialLength_for_CDFs)
    else
        sub_dat2(timeBinsForReaching>maxTrialLength)=0;
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
plot(timeBinsForReaching,fifthPerc,'Color','r');
plot(timeBinsForReaching,ninetyfifthPerc,'Color','r');

line([timeBinsForReaching(f) timeBinsForReaching(f)],[0 1],'Color','b');
c2=polyfit(timeBinsForReaching(timeBinsForReaching>=preCueWindow(1) & timeBinsForReaching<=preCueWindow(2)),d2(timeBinsForReaching>=preCueWindow(1) & timeBinsForReaching<=preCueWindow(2)),1);
yest=polyval(c2,timeBinsForReaching);
plot(timeBinsForReaching,yest,'Color','r');

c1=polyfit(timeBinsForReaching(timeBinsForReaching>=preCueWindow(1) & timeBinsForReaching<=preCueWindow(2)),d1(timeBinsForReaching>=preCueWindow(1) & timeBinsForReaching<=preCueWindow(2)),1);
yest=polyval(c1,timeBinsForReaching);
plot(timeBinsForReaching,yest,'Color','k');
ylim([0 1]);
xlim([0 maxTrialLength]);

end

function out=whetherToSuppressPlots()

out=false;

end

function [out1,out2]=plotCDF_rawReaches(data1,data2,timesteps,cueTime,tit,subtractPreCue,preCueWindow,maxTrialLength)

% make raw reaching data a timeseries (i.e., a pdf)
% then simply accumulate distribution, selecting only time points after the
% cue

% maxTrialLength=9.5; % in sec, or empty if want to include all reaches

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