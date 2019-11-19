function fitLocalCorrelationStructureOfReaches(dataset)

% RTM = regression to the mean

histo_nbins=[-4*12.4245:0.2510:4*12.4245];
checkRange_prevRT=[6 8];
checkRange_deltaRT=[-20 9];

i=1;
temp1_seq2=dataset.event_RT_trial1InSeq{i};
% temp2_seq2=dataset.event_RT_trialiInSeq{i};
% rtchanges_seq2=temp1_seq2(dataset.realrtpair_seq2{i}==1)-temp2_seq2(dataset.realrtpair_seq2{i}==1);
x2=temp1_seq2(dataset.realrtpair_seq2{i}==1);
y2=dataset.alldim_rtchanges_event{i};
nBinsPerDim={histo_nbins; histo_nbins};
% make 2D histogram
if size(x2,2)>1
    % make column vector
    x2=x2';
    y2=y2';
end
if length(nBinsPerDim)>1
    [n,bin_c]=hist3([x2 y2],'Ctrs',nBinsPerDim);
else
    [n,bin_c]=hist3([x2 y2],[nBinsPerDim nBinsPerDim]);
end
n=n./nansum(nansum(n,1),2);
figure();
imagesc(bin_c{1},bin_c{2},log(n)');
set(gca,'YDir','normal');

% RTM -- shuffled
i=1;
temp1_seq2=dataset.event_RT_trial1InSeq{i};
temp2_seq2=dataset.event_RT_trialiInSeq{i};
% shuffle the second reaction time
temp2_seq2=temp2_seq2(randperm(length(temp2_seq2)));
x2=temp1_seq2;
y2=temp1_seq2-temp2_seq2;
nBinsPerDim={histo_nbins; histo_nbins};
% make 2D histogram
if size(x2,2)>1
    % make column vector
    x2=x2';
    y2=y2';
end
if length(nBinsPerDim)>1
    [n_data,bin_c]=hist3([x2 y2],'Ctrs',nBinsPerDim);
else
    [n_data,bin_c]=hist3([x2 y2],[nBinsPerDim nBinsPerDim]);
end
n_data=n_data./nansum(nansum(n_data,1),2);
figure();
imagesc(bin_c{1},bin_c{2},log(n_data)');
set(gca,'YDir','normal');

realShiftRPE=nanmean(n(bin_c{1}>=checkRange_prevRT(1) & bin_c{1}<checkRange_prevRT(2),bin_c{2}>=checkRange_deltaRT(1) & bin_c{2}<checkRange_deltaRT(2)),1);
disp(nanmean(realShiftRPE));
shuffleShiftRPE=nanmean(n_data(bin_c{1}>=checkRange_prevRT(1) & bin_c{1}<checkRange_prevRT(2),bin_c{2}>=checkRange_deltaRT(1) & bin_c{2}<checkRange_deltaRT(2)),1);
disp(nanmean(shuffleShiftRPE));

figure();
temp=bin_c{1};
temp=temp(bin_c{2}>=checkRange_deltaRT(1) & bin_c{2}<checkRange_deltaRT(2));
plot(temp,realShiftRPE./nansum(realShiftRPE),'Color','k');
hold on;
plot(temp,shuffleShiftRPE./nansum(shuffleShiftRPE),'Color','r');
legend({'real', 'shuffled'});

end