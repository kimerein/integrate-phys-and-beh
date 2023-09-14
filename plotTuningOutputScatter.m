function plotTuningOutputScatter(tuningOutput,whichField1,whichField2,whichCuezOfUnits,timeWindowToCompare)

time1=tuningOutput.(whichField1).time{whichCuezOfUnits};
time2=tuningOutput.(whichField2).time{whichCuezOfUnits};

unitsfr1=tuningOutput.(whichField1).allunits{whichCuezOfUnits};
whichinds=time1>=timeWindowToCompare(1) & time1<=timeWindowToCompare(2);
fr1=nanmean(unitsfr1(:,whichinds),2);

unitsfr2=tuningOutput.(whichField2).allunits{whichCuezOfUnits};
whichinds=time2>=timeWindowToCompare(1) & time2<=timeWindowToCompare(2);
fr2=nanmean(unitsfr2(:,whichinds),2);

figure(); 
plot(time1,nanmean(unitsfr1),'Color','r'); hold on;
plot(time2,nanmean(unitsfr2),'Color','b'); 

figure();
scatter(ones(size(fr1)),fr1,[],'r');
hold on;
scatter(2*ones(size(fr2)),fr2,[],'b');
for i=1:length(fr1)
    line([1 2],[fr1(i) fr2(i)],'Color','k');
end
disp('pval from signrank');
p=signrank(fr1,fr2);
disp(p);

fr1(fr1<0.01)=0.01;
figure();
scatter(ones(size(fr1)),fr1./fr1,[],'r');
hold on;
scatter(2*ones(size(fr2)),fr2./fr1,[],'b');
for i=1:length(fr1)
    line([1 2],[fr1(i)/fr1(i) fr2(i)/fr1(i)],'Color','k');
end
title('normalized to fr1');

figure();
[n,x]=histcounts(fr2./fr1,0-0.25:0.5:15);
[n,x]=cityscape_hist(n,x);
figure(); plot(x,n,'Color','k'); hold on;
line([1 1],[0 15]);
title('Histogram of fr2 / fr1');

figure();
[n,x]=histcounts(fr2-fr1,-15-0.1:0.2:15);
[n,x]=cityscape_hist(n,x);
figure(); plot(x,n,'Color','k'); hold on;
line([0 0],[0 15]);
title('Histogram of fr2 - fr1');

end