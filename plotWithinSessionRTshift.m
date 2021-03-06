function plotWithinSessionRTshift(dataset,alltbt,trialTypes,metadata,settings,considerWindowAfterCue)

n_trials_away=settings.n_trials_away;
cueName=settings.cueName;
useReachType=settings.useReachType;
bins=settings.bins;
n_sems=settings.n_sems;
baselineWindow=settings.baselineWindow;

backToBase=getReturnToBaseline(dataset,alltbt,cueName,n_sems,baselineWindow,n_trials_away,false);
backToBase=4.55;

% Get RT shift within session
% Also consider rate change within this session
% If rate change is too dominant, note or throw out this data

sessids=unique(metadata.sessid);

usedPairs=dataset.event_isSeq{n_trials_away};
isTrial1=dataset.templateSequence2_cond;
isTrial1_inds=find(isTrial1);
isTrial1_inds=isTrial1_inds(logical(dataset.realrtpair_seq2{n_trials_away}'));  % check that this is a real pair (i.e., doesn't span two videos or sessions)
isTrial1=zeros(size(isTrial1));
isTrial1(isTrial1_inds)=1;
isTrial2_inds=isTrial1_inds+1;
isTrial2=zeros(size(isTrial1));
isTrial2(isTrial2_inds)=1;
allrts1=dataset.event_RT_trial1InSeq{n_trials_away};
allrts2=dataset.event_RT_trialiInSeq{n_trials_away};

if size(alltbt.all_reachBatch,1)~=length(usedPairs)
    error('alltbt and dataset do not match');    
end

cueInd=find(nanmean(alltbt.(cueName),1)>0,1,'first');
timeStep=mode(diff(nanmean(alltbt.times,1)));
timeBinsForReaching=0:timeStep:(size(alltbt.times,2)-1)*timeStep;
cueTime=timeBinsForReaching(cueInd);
backToBaseline=backToBase-cueTime; % wrt cue onset
considerWindowAfterCueInInds(1)=cueInd+floor(considerWindowAfterCue(1)/timeStep);
considerWindowAfterCueInInds(2)=cueInd+ceil(considerWindowAfterCue(2)/timeStep);

rates_trial1=nan(1,length(sessids));
rates_trial2=nan(1,length(sessids));
rateChanges=nan(1,length(sessids));
earth_movers=nan(1,length(sessids));
early_earth_movers=nan(1,length(sessids));
late_earth_movers=nan(1,length(sessids));
dprimes=nan(1,length(sessids));
rt1_cdfs=nan(length(sessids),length(bins));
rt2_cdfs=nan(length(sessids),length(bins));
rt_cdf_mids=nan(length(sessids),length(bins));
for i=1:length(sessids)
    currsess=sessids(i);
    isTrialinSess=metadata.sessid==currsess;
    indsInSess=find(isTrialinSess);
    dprimes(i)=nanmean(metadata.dprimes(isTrialinSess));
    % Consider rate change within this session
    rates_trial1(i)=getRate(alltbt,useReachType,isTrial1 & isTrialinSess,considerWindowAfterCueInInds,timeStep);
    rates_trial2(i)=getRate(alltbt,useReachType,isTrial2 & isTrialinSess,considerWindowAfterCueInInds,timeStep);
    rateChanges(i)=rates_trial1(i)-rates_trial2(i);
    
    rts1=allrts1(ismember(isTrial1_inds,indsInSess));
    rts2=allrts2(ismember(isTrial2_inds,indsInSess));
    [earth_movers(i),rt1_cdfs(i,:),rt2_cdfs(i,:),rt_cdf_mids(i,:),early_earth_movers(i),late_earth_movers(i),rt1_cdfs_late(i,:),rt2_cdfs_late(i,:),rt_cdf_mids_late(i,:)]=getRTshift(rts1,rts2,[bins(1)-(bins(2)-bins(1)) bins],backToBaseline);
end

% Plot results
plotTimeseries(rt1_cdfs,[],'k',rt2_cdfs,[],'r',nanmean(rt_cdf_mids,1),[],[],true);
title('RT CDF');
plotTimeseries(rt1_cdfs_late,[],'k',rt2_cdfs_late,[],'r',nanmean(rt_cdf_mids_late,1),[],[],true);
title('RT CDF only late RTs after cued reaching');

% figure();
% scatter(dprimes,earth_movers);
% title('dprimes vs earth movers all RTs');
% plotHist(earth_movers,earth_movers,-1-0.025:0.05:1,'Histogram of earth movers','Earth mover');

% figure();
% scatter(dprimes,early_earth_movers-late_earth_movers);
% title('dprimes vs earth movers early minus earth movers late');
% plotHist(early_earth_movers-late_earth_movers,early_earth_movers-late_earth_movers,-0.1-0.0005:0.001:0.1,'Histogram of earth movers early minus late','Earth mover');
% disp(['P-val of signrank test of earth movers early minus late is ' num2str(signrank(early_earth_movers-late_earth_movers))]);

figure();
% cutOff=0.1;
cutOff=100;
scatter(dprimes(rateChanges<=cutOff),early_earth_movers(rateChanges<=cutOff)-late_earth_movers(rateChanges<=cutOff));
disp('Using n sessions');
disp(nansum(rateChanges<=cutOff));
title('Only sessions with little rate change');
figure(); scatter(rateChanges,early_earth_movers-late_earth_movers);
disp(['P-val of signrank test of earth movers early minus late no rate change only is ' num2str(signrank(early_earth_movers(rateChanges<=cutOff)-late_earth_movers(rateChanges<=cutOff)))]);
plotHist(early_earth_movers(rateChanges<=cutOff),early_earth_movers(rateChanges<=cutOff),-0.1-0.0001:0.001:0.1,'Histogram of earth movers early','Earth mover');
% plotHist(early_earth_movers(rateChanges<=cutOff)-late_earth_movers(rateChanges<=cutOff),early_earth_movers(rateChanges<=cutOff)-late_earth_movers(rateChanges<=cutOff),-0.1-0.0005:0.001:0.1,'Histogram of earth movers early minus late','Earth mover');
plotTimeseries(rt1_cdfs(rateChanges<=cutOff,:),[],'k',rt2_cdfs(rateChanges<=cutOff,:),[],'r',nanmean(rt_cdf_mids,1),[],[],true);
earth_movers_early=nan(1,length(sessids));
earlyRTs1=[];
earlyRTs2=[];
for i=1:length(sessids)
    currsess=sessids(i);
    isTrialinSess=metadata.sessid==currsess;
    indsInSess=find(isTrialinSess);
    rts1=allrts1(ismember(isTrial1_inds,indsInSess));
    rts2=allrts2(ismember(isTrial2_inds,indsInSess));
    temp=[bins(1)-(bins(2)-bins(1)):(bins(2)-bins(1))/100:bins(end)];
%     temp=[bins(1)-(bins(2)-bins(1)):(bins(2)-bins(1))/1:bins(end)];
    [earth_movers_early(i),rt1_cdfs_early(i,:),rt2_cdfs_early(i,:),rt_cdf_mids_early(i,:),~,~,~,~,~,rtpdfs1(i,:),rtpdfs2(i,:)]=getRTshift(rts1,rts2,temp(temp<backToBaseline),backToBaseline);
    if rateChanges(i)<=cutOff
        earlyRTs1=[earlyRTs1 rts1(rts1<backToBaseline)];
        earlyRTs2=[earlyRTs2 rts2(rts2<backToBaseline)];
    end
end
plotTimeseries(rt1_cdfs_early(rateChanges<=cutOff,:),[],'k',rt2_cdfs_early(rateChanges<=cutOff,:),[],'r',nanmean(rt_cdf_mids_early,1),[],[],false);
[~,p]=kstest2(earlyRTs1,earlyRTs2);
disp(['P-val of 2 sample KS test comparing early RTs for trial 1 vs trial 2 is ' num2str(p)]);

plotTimeseries(rtpdfs1(rateChanges<=cutOff,:),[],'k',rtpdfs2(rateChanges<=cutOff,:),[],'r',nanmean(rt_cdf_mids_early,1),[],[],true);

end

function [x_backup,returnThis,f]=plotHist(data1,data2,bins,tit,xlab)

useLogForY=false;

[n,x]=histcounts(data1,bins);
x_backup=x;
[n,x]=cityscape_hist(n,x);
f=figure();
if useLogForY==true
    semilogy(x,n./nansum(n),'Color','k');
else
    plot(x,n./nansum(n),'Color','k');
end
xlabel(xlab);
ylabel('Count');
title(tit);

[n,x]=histcounts(data2,x_backup);
[n,x]=cityscape_hist(n,x);
hold on;
if useLogForY==true
    semilogy(x,n./nansum(n),'Color','r');
else
    plot(x,n./nansum(n),'Color','r');
end
leg={'data1','data2'};
legend(leg);

returnThis.x=x;
returnThis.y=n./nansum(n);

end

function [earth_mover,x_cdf,x_cdf_cond2,x_mids,earth_mover_early,earth_mover_late,x_cdf_late,x_cdf_cond2_late,x_mids_late,rtpdfs1,rtpdfs2]=getRTshift(rts1,rts2,bins,timeSplit)

[n,x]=histcounts(rts1,bins);
rtpdfs1=n./nansum(n);
x_mids=nanmean([x(1:end-1); x(2:end)],1);
x_cdf=accumulateDistribution(n);
x_cdf=x_cdf./nanmax(x_cdf);
x_backup=x;
[n,x]=histcounts(rts2,bins);
rtpdfs2=n./nansum(n);
x_cdf_cond2=accumulateDistribution(n);
x_cdf_cond2=x_cdf_cond2./nanmax(x_cdf_cond2);

differ=x_cdf_cond2-x_cdf; % so that positive is a speed up in RT in rts2
earth_mover=nansum(differ)*(bins(2)-bins(1));

[~,binForTimeSplit]=nanmin(abs(x_mids-timeSplit));
differ=x_cdf_cond2(1:binForTimeSplit)-x_cdf(1:binForTimeSplit);
earth_mover_early=nansum(differ)*(bins(2)-bins(1))/length(differ);

[n]=histcounts(rts1,bins);
x_cdf_late=accumulateDistribution(n(x_mids>=timeSplit));
x_cdf_late=x_cdf_late./nanmax(x_cdf_late);
[n]=histcounts(rts2,bins);
x_cdf_cond2_late=accumulateDistribution(n(x_mids>=timeSplit));
x_cdf_cond2_late=x_cdf_cond2_late./nanmax(x_cdf_cond2_late);
x_mids_late=x_mids(x_mids>timeSplit);

differ=x_cdf_cond2_late-x_cdf_late;
earth_mover_late=nansum(differ)*(bins(2)-bins(1))/length(differ);

end

function rate=getRate(alltbt,useReachType,theseTrials,considerWindowAfterCueInInds,timeStep)

temp=alltbt.(useReachType);
temp=temp(theseTrials,:);

rate=nanmean(nanmean(temp(:,considerWindowAfterCueInInds(1):considerWindowAfterCueInInds(2)),1)./timeStep,2);

end

function timeInd=getReturnToBaseline(dataset,alltbt,cueName,n_sems,baselineWindow,n_trials_away,suppressOutput)

timeStep=mode(diff(nanmean(alltbt.times,1)));
timeBinsForReaching=0:timeStep:(size(dataset.rawReaching_event_trial1InSeq{n_trials_away},2)-1)*timeStep;

% find cue
cueInd=find(nanmean(alltbt.(cueName),1)>0,1,'first');
cueInd=cueInd+floor(0.25/timeStep); % wait at least 250 ms after cue
baselineInds=timeBinsForReaching>baselineWindow(1) & timeBinsForReaching<baselineWindow(2);

% find when reach distribution returns to baseline after cue
% within n_sems of baseline
reachMeans=dataset.rawReaching_event_trial1InSeq{n_trials_away};
reachSEMs=dataset.se_rawReaching_event_trial1InSeq{n_trials_away};
baseline=nanmean(nanmean(reachMeans(:,baselineInds),1),2);
semAtBaseline=nanmean(sqrt(nansum(reachSEMs(:,baselineInds).^2,1)));
reachesAfterCue=nanmean(reachMeans(:,cueInd:end),1);
firstAtBaseAfterCue=find(reachesAfterCue<=baseline+n_sems*semAtBaseline,1,'first');
timeInd=timeBinsForReaching(cueInd+firstAtBaseAfterCue-1);

if suppressOutput==false
    i=n_trials_away;
    plotTimeseries(dataset.rawReaching_event_trial1InSeq{i},dataset.se_rawReaching_event_trial1InSeq{i},'k',dataset.rawReaching_event_trialiInSeq{i},dataset.se_rawReaching_event_trialiInSeq{i},'m',timeBinsForReaching,baseline+n_sems*semAtBaseline,[timeInd baseline+n_sems*semAtBaseline],true);
    temp=dataset.event_name;
    temp(regexp(temp,'_'))=' ';
    title(['Fx of ' temp ' first trial (black) vs trial ' num2str(dataset.nInSequence(i)-1) ' later (magenta)']);
    %legend({'me+-se','first trial','','','me+-se',['trial ' num2str(dataset.nInSequence(i)-1) ' later'],'',''});
    legend({'first trial','','',['trial ' num2str(dataset.nInSequence(i)-1) ' later'],'','','baseline','returnToBase'});
end
      
end

function plotTimeseries(data1,data1_se,color1,data2,data2_se,color2,timeBins,lineVal,scatterVal,plotAsCityscape)

% plotAsCityscape=true;

data1_mean=nanmean(data1,1);
data2_mean=nanmean(data2,1);
if isempty(data1_se)
    data1_se=nanstd(data1,[],1)./sqrt(size(data1,1));
else
    data1_se=sqrt(nansum(data1_se.^2,1));
end
if isempty(data2_se)
    data2_se=nanstd(data2,[],1)./sqrt(size(data2,1));
else
    data2_se=sqrt(nansum(data2_se.^2,1));
end

figure();
%fill([timeBins fliplr(timeBins)],[data1_mean+data1_se fliplr(data1_mean-data1_se)],[0.5 0.5 0.5]);
%hold on;
if plotAsCityscape==true
    [n,x]=cityscape_hist(data1_mean,timeBins);
    plot(x,n,'Color',color1,'LineWidth',1); hold on;
    [n,x]=cityscape_hist(data1_mean+data1_se,timeBins);
    plot(x,n,'Color',color1,'LineWidth',0.1);
    [n,x]=cityscape_hist(data1_mean-data1_se,timeBins);
    plot(x,n,'Color',color1,'LineWidth',0.1);
else
    plot(timeBins,data1_mean,'Color',color1); hold on;
%     plot(timeBins,data1_mean+data1_se,'Color',color1);
%     plot(timeBins,data1_mean-data1_se,'Color',color1);
end

%fill([timeBins fliplr(timeBins)],[data2_mean+data2_se fliplr(data2_mean-data2_se)],[0.1 0.7 0.5]);
%hold on;
if plotAsCityscape==true
    [n,x]=cityscape_hist(data2_mean,timeBins);
    plot(x,n,'Color',color2,'LineWidth',1); hold on;
    [n,x]=cityscape_hist(data2_mean+data2_se,timeBins);
    plot(x,n,'Color',color2,'LineWidth',0.1);
    [n,x]=cityscape_hist(data2_mean-data2_se,timeBins);
    plot(x,n,'Color',color2,'LineWidth',0.1);
else
    plot(timeBins,data2_mean,'Color',color2); hold on;
%     plot(timeBins,data2_mean+data2_se,'Color',color2);
%     plot(timeBins,data2_mean-data2_se,'Color',color2);
end

if ~isempty(lineVal)
    line([timeBins(1) timeBins(end)],[lineVal lineVal],'Color','c');
end
if ~isempty(scatterVal)
    scatter(scatterVal(1),scatterVal(2));
end

end
