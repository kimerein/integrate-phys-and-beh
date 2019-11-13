function fitRateTermToLateReaches(dataset)

rtBins=[6 9];
maxTrialLength=9.5; % in seconds
histo_nbins=[-4*12.4245:0.2510:4*12.4245];
i=1;
j=1;
temp1=dataset.allTrialsSequence_RT_trial1InSeq{i};
temp1_seq2=dataset.event_RT_trial1InSeq{i};
dimall_all=dataset.alldim_rtchanges_allTrialsSequence{i};
dimall_event=dataset.alldim_rtchanges_event{i};
a=temp1(dataset.realrtpair_seq1{i}==1);
b=temp1_seq2(dataset.realrtpair_seq2{i}==1);
temp_RTdiffs1=dimall_all(a>=rtBins(j,1) & a<rtBins(j,2));
temp_RTdiffsevent=dimall_event(b>=rtBins(j,1) & b<rtBins(j,2));
[histo_nbins,RTchangeHist]=plotHist(temp_RTdiffs1,temp_RTdiffsevent,histo_nbins,['CDF RT change for RTs ' num2str(dataset.nInSequence(i)-1) ' trials later, RTs less than ' num2str(rtBins(j,2))],'Change in RT');

temp=dataset.event_name;
temp(regexp(temp,'_'))=' ';
a=temp1(dataset.realrtpair_seq1{i}==1);
b=temp1_seq2(dataset.realrtpair_seq2{i}==1);
[histo_nbins,RTpdf]=plotHist(a(a>=rtBins(j,1) & a<=rtBins(j,2)),b(b>=rtBins(j,1) & b<=rtBins(j,2)),histo_nbins,['Reaction Times fx of ' temp ': trial 1 (black) vs ' num2str(dataset.nInSequence(i)-1) ' later (magenta)'],'RT (sec)');

bestTau=fitSingleTau(temp_RTdiffs1,a(a>=rtBins(j,1) & a<=rtBins(j,2)),[rtBins(2)-maxTrialLength rtBins(2)],histo_nbins);
disp(bestTau);

bestTau=fitSingleTau(temp_RTdiffsevent,b(b>=rtBins(j,1) & b<=rtBins(j,2)),[rtBins(2)-maxTrialLength rtBins(2)],histo_nbins);
disp(bestTau);

rtBins=[6 9];
[slope_early,slope_late]=fitTwoTaus(dimall_all(a>=rtBins(j,1) & a<rtBins(j,2)),histo_nbins);
disp(slope_early);
disp(slope_late);
x=0:0.001:20;
y=randsample(x,length(dimall_all(a>=rtBins(j,1) & a<rtBins(j,2))),true,-(exp(-abs(slope_late).*x)-exp(-abs(slope_early).*x)));
plotHist(a(a>=rtBins(j,1) & a<=rtBins(j,2))-dimall_all(a>=rtBins(j,1) & a<rtBins(j,2)),y,histo_nbins,'Fitting tau','RT (sec)');

[slope_early,slope_late]=fitTwoTaus(dimall_event(b>=rtBins(j,1) & b<rtBins(j,2)),histo_nbins);
disp(slope_early);
disp(slope_late);
x=0:0.001:20;
y=randsample(x,length(dimall_event(b>=rtBins(j,1) & b<rtBins(j,2))),true,-(exp(-abs(slope_late).*x)-exp(-abs(slope_early).*x)));
plotHist(b(b>=rtBins(j,1) & b<=rtBins(j,2))-dimall_event(b>=rtBins(j,1) & b<rtBins(j,2)),y,histo_nbins,'Fitting tau','RT (sec)');

end

function [slope_early,slope_late]=fitTwoTaus(dataToFit,histo_nbins)

[n_data,xhist]=histcounts(dataToFit,histo_nbins);
x_log=nanmean([xhist(1:end-1); xhist(2:end)],1);
n_data=n_data./nansum(n_data);
temp=n_data;
n_data=n_data(find(temp>0,1,'first'):end);
x_log=x_log(find(temp>0,1,'first'):end);
temp=temp(find(temp>0,1,'first'):end);
n_data=n_data(1:find(temp>0,1,'last'));
x_log=x_log(1:find(temp>0,1,'last'));
temp=temp(1:find(temp>0,1,'last'));
log_n_data=log(n_data);

% Fit early and late sections
figure();
plot(log_n_data);
earlyIndEnd=input('Last index for early exponent ');
lateIndBegin=input('First index for late exponent ');
% deriv=diff(smooth(log_n_data,5));
% deriv(isinf(deriv))=nan;
% earlyIndEnd=find(deriv<0,1,'first');
% lateIndBegin=find(deriv>0,1,'last');

% Fit late tau
log_n_data(isinf(log_n_data))=nan;
slope_early=nanmean(diff(log_n_data(1:earlyIndEnd))./diff(x_log(1:earlyIndEnd)));
slope_late=nanmean(diff(log_n_data(lateIndBegin:end))./diff(x_log(lateIndBegin:end)));

figure();
plot(x_log,log_n_data,'Color','k');
hold on;
plot(x_log(~isnan(log_n_data)),slope_early.*x_log(~isnan(log_n_data))+nanmean(log_n_data(1:earlyIndEnd)),'Color','r');
plot(x_log(~isnan(log_n_data)),slope_late.*x_log(~isnan(log_n_data))+log_n_data(end)-slope_late.*x_log(end),'Color','g');

end

function bestTau=fitSingleTau(dataToFit,firstRTs,fitRange,histo_nbins)

tauRange=0.0001:0.01:10;

x=0:0.001:20;
[n_data,xhist]=histcounts(dataToFit,histo_nbins);
n_data=n_data./nansum(n_data);
differs=nan(1,length(tauRange));
for i=1:length(tauRange)
    currTau=tauRange(i);
    y=randsample(x,length(firstRTs),true,exp(-currTau.*x));
    [n,xhist]=histcounts(firstRTs-y,histo_nbins);
    n=n./nansum(n);
    differ=nansum(abs(n_data(xhist>=fitRange(1) & xhist<=fitRange(2))-n(xhist>=fitRange(1) & xhist<=fitRange(2))));
    differs(i)=differ;
end
[~,mi]=nanmin(differs);
bestTau=tauRange(mi);
y=randsample(x,length(firstRTs),true,exp(-bestTau.*x));
plotHist(dataToFit,firstRTs-y,histo_nbins,'Fitting tau','RT (sec)');

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