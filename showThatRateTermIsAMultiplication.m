function showThatRateTermIsAMultiplication()

bin_centers = @(bins) nanmean([bins(1:end-1); bins(2:end)],1);

y=-exp(-bin_centers(bins)./3).*1;
randsample(bin_centers(bins),5000,true,y



return

reachesPerTrial=5;
nTrials=5000;
bins=0:0.01:9;

reachTimes_fromReachRateOnly=0+(0+9)*rand(nTrials,reachesPerTrial);

rate=reachesPerTrial;
bin_centers = @(bins) nanmean([bins(1:end-1); bins(2:end)],1);
for i=1:nTrials
    reachTimes(i,:)=randsample(bin_centers(bins),1,true,gampdf(bin_centers(bins),1,1/rate));
end

[~,reachPlot]=plotHist(reachTimes(1:end),reachTimes(1:end),bins,'PDF reaching','Time (sec)');

RTs=getRTs(reachTimes);

[~,RTplot]=plotHist(RTs,RTs,bins,'PDF RTs','Time (sec)');

figure();
plot(reachPlot.x,RTplot.y./reachPlot.y);
title('PDF RT / PDF reaching');


[~,reachPlot]=plotHist(reachTimes_fromReachRateOnly(1:end),reachTimes_fromReachRateOnly(1:end),bins,'PDF reaching from reach rate only','Time (sec)');

RTs=getRTs(reachTimes_fromReachRateOnly);

[~,RTplot]=plotHist(RTs,RTs,bins,'PDF RTs from reach rate only','Time (sec)');

figure();
plot(reachPlot.x(1:2:end-1),gampdf(bin_centers(bins),1,1/rate).*RTplot.y(1:2:end-1));
title('PDF reaching * PDF wait time');

end

function RTs=getRTs(reachTimes)

RTs=nan(1,size(reachTimes,1));
for i=1:size(reachTimes,1)
    RTs(i)=nanmin(reachTimes(i,:));
end

end

function [x_backup,returnThis]=plotHist(data1,data2,bins,tit,xlab)

useLogForY=false;

[n,x]=histcounts(data1,bins);
x_backup=x;
[n,x]=cityscape_hist(n,x);
figure();
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