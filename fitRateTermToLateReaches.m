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

bestTau=fitTau(temp_RTdiffs1,a(a>=rtBins(j,1) & a<=rtBins(j,2)),[rtBins(2)-maxTrialLength rtBins(2)],histo_nbins);
disp(bestTau);

bestTau=fitTau(temp_RTdiffsevent,b(b>=rtBins(j,1) & b<=rtBins(j,2)),[rtBins(2)-maxTrialLength rtBins(2)],histo_nbins);
disp(bestTau);

end

function bestTau=fitTau(dataToFit,firstRTs,fitRange,histo_nbins)

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