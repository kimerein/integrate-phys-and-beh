function plotBehaviorEventFx(dataset,alltbt)

% dataset contains ...
% Get raw reaching (all reaches)
% Get reaction times
% Get change in reaction times (non-corrected input distributions)
% Get change in reaction times (corrected input distributions)

% Plot raw reaching data
timeStep=mode(diff(nanmean(alltbt.times,1)));
timeBinsForReaching=0:timeStep:(size(dataset.rawReaching_allTrialsSequence_trial1InSeq{1},2)-1)*timeStep;
for i=1:length(dataset.rawReaching_allTrialsSequence_trial1InSeq)
    plotTimeseries(dataset.rawReaching_allTrialsSequence_trial1InSeq{i},dataset.se_rawReaching_allTrialsSequence_trial1InSeq{i},'k',dataset.rawReaching_allTrialsSequence_trialiInSeq{i},dataset.se_rawReaching_allTrialsSequence_trialiInSeq{i},'m',timeBinsForReaching);
    title(['Reference data (all trials) first trial (black) vs trial ' num2str(dataset.nInSequence(i)-1) ' later (magenta)']);
    %legend({'me+-se','first trial','','','me+-se',['trial ' num2str(dataset.nInSequence(i)-1) ' later'],'',''});
    legend({'first trial','','',['trial ' num2str(dataset.nInSequence(i)-1) ' later'],'',''});
end
for i=1:length(dataset.rawReaching_event_trial1InSeq)
    plotTimeseries(dataset.rawReaching_event_trial1InSeq{i},dataset.se_rawReaching_event_trial1InSeq{i},'k',dataset.rawReaching_event_trialiInSeq{i},dataset.se_rawReaching_event_trialiInSeq{i},'m',timeBinsForReaching);
    temp=dataset.event_name;
    temp(regexp(temp,'_'))=' ';
    title(['Fx of ' temp ' first trial (black) vs trial ' num2str(dataset.nInSequence(i)-1) ' later (magenta)']);
    %legend({'me+-se','first trial','','','me+-se',['trial ' num2str(dataset.nInSequence(i)-1) ' later'],'',''});
    legend({'first trial','','',['trial ' num2str(dataset.nInSequence(i)-1) ' later'],'',''});
end
%         
% i=1;
% temp1=forHists_cond1{i};
% temp2=forHists_cond2{i};
% plotHist(temp1,temp2,200,'Histo','y-vals');
% plotCDF(temp1,temp2,200,'CDF');
% range_learning=2; % in seconds
% plotHist(temp1(temp1>-range_learning & temp1<range_learning),temp2(temp2>-range_learning & temp2<range_learning),200,'Close-up Histo','y-vals');
% plotCDF(temp1(temp1>-range_learning & temp1<range_learning),temp2(temp2>-range_learning & temp2<range_learning),200,'Close-up CDF');
% disp('CLOSE-UP P-VAL');
% testRanksum(temp1(temp1>-range_learning & temp1<range_learning),temp2(temp2>-range_learning & temp2<range_learning),1);
% 
% 
% for i=1:length(nInSequence)
%     disp([num2str(nanmean(ns.cond1(i,:),2)) ' trials for cond1 of nInSequence ' num2str(nInSequence(i))]);
%     disp([num2str(nanmean(ns.cond2(i,:),2)) ' trials for cond2 of nInSequence ' num2str(nInSequence(i))]);
% end

end

function plotTimeseries(data1_mean,data1_se,color1,data2_mean,data2_se,color2,timeBins)

data1_mean=nanmean(data1_mean,1);
data2_mean=nanmean(data2_mean,1);
data1_se=sqrt(nansum(data1_se.^2,1));
data2_se=sqrt(nansum(data2_se.^2,1));

figure();
%fill([timeBins fliplr(timeBins)],[data1_mean+data1_se fliplr(data1_mean-data1_se)],[0.5 0.5 0.5]);
%hold on;
plot(timeBins,data1_mean,'Color',color1); hold on;
plot(timeBins,data1_mean+data1_se,'Color',color1);
plot(timeBins,data1_mean-data1_se,'Color',color1);

%fill([timeBins fliplr(timeBins)],[data2_mean+data2_se fliplr(data2_mean-data2_se)],[0.1 0.7 0.5]);
%hold on;
plot(timeBins,data2_mean,'Color',color2); hold on;
plot(timeBins,data2_mean+data2_se,'Color',color2);
plot(timeBins,data2_mean-data2_se,'Color',color2);   

end

function plotScatter(RT_pairs1_x,RT_pairs1_y,RT_pairs2_x,RT_pairs2_y,jitter1,jitter2,tit)

figure();
if length(jitter1)==1
    useX=RT_pairs1_x+rand(size(RT_pairs1_x)).*jitter;
    useY=RT_pairs1_y+rand(size(RT_pairs1_y)).*jitter;
else
    useX=RT_pairs1_x+jitter1;
    useY=RT_pairs1_y+jitter1;
end
    
s=scatter(useX,useY,100,'k','filled');
set(gcf,'position',[10,10,1000,1000]);
pause;
m=get(s,'MarkerHandle');
alpha=0.3;
m.FaceColorData=uint8(255*[0;0;0;alpha]);
xlabel('real cue');
ylabel('fake cue');
title(tit);
% xlim([0 9.5]);
% ylim([0 9.5]);

hold on;
if length(jitter2)==1
    useX=RT_pairs2_x+rand(size(RT_pairs2_x)).*jitter;
    useY=RT_pairs2_y+rand(size(RT_pairs2_y)).*jitter;
else
    useX=RT_pairs2_x+jitter2;
    useY=RT_pairs2_y+jitter2;
end
s=scatter(useX,useY,100,'r','filled');
pause;
m=get(s,'MarkerHandle');
alpha=0.3;
m.FaceColorData=uint8(255*[1;0;0;alpha]);

end

function [p,med1,med2]=testRanksum(data1,data2,dispStuff)

med1=nanmedian(data1,2);
med2=nanmedian(data2,2);

if dispStuff==1
    disp('median of data1');
    disp(med1);
    disp('median of data2');
    disp(med2);
end

if all(isnan(data1) | all(isnan(data2)))
    p=nan;
    return
end
[p,h]=ranksum(data1,data2);
if dispStuff==1
    disp('p-value');
    disp(p);
end

end

function plotCDF(data1,data2,bins,tit)

[n,x]=histcounts(data1,bins);
x_backup=x;
x_mids=nanmean([x(1:end-1); x(2:end)],1);
cond1_cdf=accumulateDistribution(n);
figure();
plot(x_mids,cond1_cdf./nanmax(cond1_cdf),'Color','k');
xlabel('CDF');
ylabel('Count');
title(tit);

[n,x]=histcounts(data2,x_backup);
hold on;
cond2_cdf=accumulateDistribution(n);
plot(x_mids,cond2_cdf./nanmax(cond2_cdf),'Color','r');

end

function cdf=accumulateDistribution(data)

cdf=nan(size(data));
for i=1:length(data)
    cdf(i)=nansum(data(1:i));
end

end

function plotHist(data1,data2,bins,tit,xlab)

[n,x]=histcounts(data1,bins);
x_backup=x;
[n,x]=cityscape_hist(n,x);
figure();
plot(x,n./nansum(n),'Color','k');
xlabel(xlab);
ylabel('Count');
title(tit);

[n,x]=histcounts(data2,x_backup);
[n,x]=cityscape_hist(n,x);
hold on;
plot(x,n./nansum(n),'Color','r');
leg={'data1','data2'};
legend(leg);

end