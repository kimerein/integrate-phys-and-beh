function compare_conditions_delta_rt_asFunc_rt(dataset1,dataset2)

% nBinsFor2Dhist=200;
nBinsFor2Dhist={[0:0.25:14]+0.125, [-14:0.25:14]+0.125};
% nBinsFor2Dhist={[0:0.125:14]+0.0625, [-14:0.125:14]+0.0625};
RPE_slice_at=[0 0.25]; % bin of delta_rts where expect to see RPE
% nonRPE_slice_at=[-0.75 -0.5]; % bin of delta_rts where do not expect to see RPE
nonRPE_slice_at=[-14 -0.5]; % bin of delta_rts where do not expect to see RPE

i=1;
dataset=dataset1;
temp1=dataset.allTrialsSequence_RT_trial1InSeq{i};
temp1=temp1(dataset.realrtpair_seq1{i}==1);
temp1_seq2=dataset.event_RT_trial1InSeq{i};
temp1_seq2=temp1_seq2(dataset.realrtpair_seq2{i}==1);
temp2_seq2=dataset.event_RT_trialiInSeq{i};
temp2_seq2=temp2_seq2(dataset.realrtpair_seq2{i}==1);
rtchanges_seq2=temp1_seq2-temp2_seq2;
cond1.temp1=temp1;
cond1.temp1_seq2=temp1_seq2;
cond1.rtchanges_seq2=rtchanges_seq2;

dataset=dataset2;
temp1=dataset.allTrialsSequence_RT_trial1InSeq{i};
temp1=temp1(dataset.realrtpair_seq1{i}==1);
temp1_seq2=dataset.event_RT_trial1InSeq{i};
temp1_seq2=temp1_seq2(dataset.realrtpair_seq2{i}==1);
temp2_seq2=dataset.event_RT_trialiInSeq{i};
temp2_seq2=temp2_seq2(dataset.realrtpair_seq2{i}==1);
rtchanges_seq2=temp1_seq2-temp2_seq2;
cond2.temp1=temp1;
cond2.temp1_seq2=temp1_seq2;
cond2.rtchanges_seq2=rtchanges_seq2;

[cond1_n,cond2_n,diff_n,x,y]=compareWithHeatmaps(cond1.temp1_seq2,cond1.rtchanges_seq2,cond2.temp1_seq2,cond2.rtchanges_seq2,nBinsFor2Dhist,'Comparing cond1 and cond2');

% Plot pdf slices for looking at RPE component
sliceDiff1=compareSlices(x,y,cond1_n,RPE_slice_at,nonRPE_slice_at,'r','b');
sliceDiff2=compareSlices(x,y,cond2_n,RPE_slice_at,nonRPE_slice_at,'r','b');

[~,histRT1,histRT2]=plotHist(cond1.temp1_seq2,cond2.temp1_seq2,nBinsFor2Dhist{1},'RTs of first reach in sequence','RT (sec)');

% Fit trial 1 RT distribution to slice difference in order to estimate RPE component
removeMeanRegression(cond1.temp1_seq2,cond1.temp1,cond1.temp1_seq2,cond1.rtchanges_seq2,0.05,0.5,'Event Cond 1',nBinsFor2Dhist); % get bootstrapped distribution from event then all trials pdf
removeMeanRegression(cond2.temp1_seq2,cond2.temp1,cond2.temp1_seq2,cond2.rtchanges_seq2,0.05,0.5,'Event Cond 2',nBinsFor2Dhist); % get bootstrapped distribution from event then all trials pdf

end

function [diff2Dhist,x,y]=removeMeanRegression(rts1,rts2,actual_first_rts,actual_rt_changes,scatterJitter,alpha,tit,nBinsFor2Dhist)

% generate a distribution from reaction time pdf
% bootstrap to sample change in rts
% note that we are not taking actual paired reaction times, take instead
% fake reaction time pairs that would result from a stationary rt pdf over
% the course of the experiment
% this will give the change resulting from regression to the mean
n_bootstrap=10000; % how many times to run bootstrap, will take 1 pair randomly per iteration of bootstrap
delta_rts=nan(1,n_bootstrap);
first_rt=nan(1,n_bootstrap);
second_rt=nan(1,n_bootstrap);
for i=1:n_bootstrap
    curr_rts=[rts1(randsample(length(rts1),1,true)) rts2(randsample(length(rts2),1,true))]; % with replacement
    delta_rts(i)=diff(curr_rts(end:-1:1));
    first_rt(i)=curr_rts(1);
    second_rt(i)=curr_rts(2);
end

% compare actual change in reaction times distribution to this bootstrapped
% distribution
% Scatter plot comparison
% plotScatter(first_rt,delta_rts,actual_first_rts,actual_rt_changes,scatterJitter,scatterJitter,[tit ' comparing bootstrapped vs real rt pairs'],'RT trial 1','Change in RT',[alpha/10 alpha]);
% Heatmap comparison
[diff2Dhist,x,y]=compareWithHeatmaps(first_rt,delta_rts,actual_first_rts,actual_rt_changes,nBinsFor2Dhist,tit);

end

function [x_backup,histOut1,histOut2]=plotHist(data1,data2,bins,tit,xlab)

[n,x]=histcounts(data1,bins);
x_backup=x;
[n,x]=cityscape_hist(n,x);
figure();
plot(x,n./nansum(n),'Color','k');
histOut1=n./nansum(n);
xlabel(xlab);
ylabel('Count');
title(tit);

[n,x]=histcounts(data2,x_backup);
[n,x]=cityscape_hist(n,x);
hold on;
plot(x,n./nansum(n),'Color','r');
histOut2=n./nansum(n);
leg={'data1','data2'};
legend(leg);

end

function sliceDiff=compareSlices(x,y,pdf2D,y_bin1,y_bin2,c1,c2)

[hax,slice1]=plotPDFat(x,y,pdf2D,y_bin1,c1,[]);
[hax,slice2]=plotPDFat(x,y,pdf2D,y_bin2,c2,hax);

sliceDiff=slice1-slice2;

figure();
plot(x,sliceDiff,'Color','k');
xlabel('Reaction time trial 1');
ylabel('Slice 1 minus slice 2');

end

function [hax,slice]=plotPDFat(x,y,pdf2D,y_bin,c,hax)

slice=nanmean(pdf2D(:,y>=y_bin(1) & y<=y_bin(2)),2);
if isempty(hax)
    figure();
    hax=axes();
end
plot(hax,x,slice,'Color',c);
hold on;
xlabel('Reaction time trial 1');
ylabel(['PDF of change in reaction times']);
title('Slice of PDF');

end

function [n,n2,diff_n,x,y]=compareWithHeatmaps(x1,y1,x2,y2,nBinsPerDim,tit)

% make 2D histogram
if size(x1,2)>1
    % make column vector
    x1=x1';
    y1=y1';
    x2=x2';
    y2=y2';
end
if length(nBinsPerDim)>1
    [n,bin_c]=hist3([x1 y1],'Ctrs',nBinsPerDim);
else
    [n,bin_c]=hist3([x1 y1],[nBinsPerDim nBinsPerDim]);
end
[n2]=hist3([x2 y2],'Ctrs',bin_c);

% Normalize each
n2=n2./nansum(nansum(n2,1),2);
n=n./nansum(nansum(n,1),2);

diff_n=n2-n;
figure();
imagesc(bin_c{1},bin_c{2},log(n'));
set(gca,'YDir','normal');
title([tit ' Cond 1 histogram']);
xlabel('Reaction time trial 1');
ylabel('Change in reaction times');

figure();
imagesc(bin_c{1},bin_c{2},log(n2'));
set(gca,'YDir','normal');
title([tit ' Cond 2 histogram']);
xlabel('Reaction time trial 1');
ylabel('Change in reaction times');

figure();
imagesc(bin_c{1},bin_c{2},log(diff_n'-nanmin(nanmin(diff_n'))));
set(gca,'YDir','normal');
title([tit ' Difference histogram']);
xlabel('Reaction time trial 1');
ylabel('Change in reaction times');

x=bin_c{1};
y=bin_c{2};

nSmooth=10;
K=ones(nSmooth);
smoothMat=conv2(diff_n'-nanmin(nanmin(diff_n'))',K,'same');
figure();
imagesc(x(floor(nSmooth/2):end-floor(nSmooth/2)),y(floor(nSmooth/2):end-floor(nSmooth/2)),log(smoothMat(floor(nSmooth/2):end-floor(nSmooth/2),floor(nSmooth/2):end-floor(nSmooth/2))));
set(gca,'YDir','normal');
title(['Smoothed']);
xlabel('Reaction time trial 1');
ylabel('Change in reaction times');

end