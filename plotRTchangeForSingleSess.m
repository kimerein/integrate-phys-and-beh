function plotRTchangeForSingleSess(dataset,alltbt,metadata,trialTypes)

whichsess=unique(metadata.sessid(trialTypes.dprime<1.6));
%whichsess=whichsess(3);
plot_rt_pdf=false;
plot_rt_cdf=false;
plot_delta_rt_asFunc_rt=false;
plot_delta_rt_asFunc_rt_removeMeanRegression=true;
plot_dim2_delta_asFunc_rt=false;
scatterJitter=0.05;
alpha=0.5;
histo_nbins=100; % number of bins for reaction time histogram
backup_histo_nbins=histo_nbins;
nBinsFor2Dhist=200;
smoothSize=15;

if isfield(trialTypes,'dprime') && length(whichsess)==1
    disp(['dprime for this session is ' num2str(unique(trialTypes.dprime(ismember(metadata.sessid,whichsess))))]);
end

% Plot reaction times distribution
if plot_rt_pdf==true
    histo_nbins=backup_histo_nbins;
    for i=1:length(dataset.allTrialsSequence_RT_trial1InSeq)
        temp1=dataset.allTrialsSequence_RT_trial1InSeq{i};
        temp1=temp1(ismember(metadata.sessid(dataset.allTrialsSequence_isSeq{i}==1),whichsess)');
        temp2=dataset.allTrialsSequence_RT_trialiInSeq{i};
        temp2=temp2(ismember(metadata.sessid(dataset.allTrialsSequence_isSeq{i}==1),whichsess)');
        histo_nbins=plotHist(temp1,temp2,histo_nbins,['Reaction Times all trials reference: trial 1 (black) vs ' num2str(dataset.nInSequence(i)-1) ' later (magenta)'],'RT (sec)');
    end
    temp=dataset.event_name;
    temp(regexp(temp,'_'))=' ';
    for i=1:length(dataset.event_RT_trial1InSeq)
        temp1=dataset.event_RT_trial1InSeq{i};
        temp1=temp1(ismember(metadata.sessid(dataset.event_isSeq{i}==1),whichsess)');
        temp2=dataset.event_RT_trialiInSeq{i};
        temp2=temp2(ismember(metadata.sessid(dataset.event_isSeq{i}==1),whichsess)');
        histo_nbins=plotHist(temp1,temp2,histo_nbins,['Reaction Times fx of ' temp ': trial 1 (black) vs ' num2str(dataset.nInSequence(i)-1) ' later (magenta)'],'RT (sec)');
    end
end

% Plot reaction times CDF
if plot_rt_cdf==true
    histo_nbins=backup_histo_nbins;
    for i=1:length(dataset.allTrialsSequence_RT_trial1InSeq)
        temp1=dataset.allTrialsSequence_RT_trial1InSeq{i};
        temp1=temp1(ismember(metadata.sessid(dataset.allTrialsSequence_isSeq{i}==1),whichsess)');
        temp2=dataset.allTrialsSequence_RT_trialiInSeq{i};
        temp2=temp2(ismember(metadata.sessid(dataset.allTrialsSequence_isSeq{i}==1),whichsess)');
        histo_nbins=plotCDF(temp1,temp2,histo_nbins,['CDF Reaction Times all trials reference: trial 1 (black) vs ' num2str(dataset.nInSequence(i)-1) ' later (magenta)']);
    end
    temp=dataset.event_name;
    temp(regexp(temp,'_'))=' ';
    for i=1:length(dataset.event_RT_trial1InSeq)
        temp1=dataset.event_RT_trial1InSeq{i};
        temp1=temp1(ismember(metadata.sessid(dataset.event_isSeq{i}==1),whichsess)');
        temp2=dataset.event_RT_trialiInSeq{i};
        temp2=temp2(ismember(metadata.sessid(dataset.event_isSeq{i}==1),whichsess)');
        histo_nbins=plotCDF(temp1,temp2,histo_nbins,['CDF Reaction Times fx of ' temp ': trial 1 (black) vs ' num2str(dataset.nInSequence(i)-1) ' later (magenta)']);
    end
end

% Plot change in RT as a function of RT
if plot_delta_rt_asFunc_rt==true
    for i=1:length(dataset.allTrialsSequence_RT_trial1InSeq)
        temp1=dataset.allTrialsSequence_RT_trial1InSeq{i};
        temp1=temp1(dataset.realrtpair_seq1{i}==1 & ismember(metadata.sessid(dataset.allTrialsSequence_isSeq{i}==1),whichsess)');
        temp2=dataset.allTrialsSequence_RT_trialiInSeq{i};
        temp2=temp2(dataset.realrtpair_seq1{i}==1 & ismember(metadata.sessid(dataset.allTrialsSequence_isSeq{i}==1),whichsess)');
        rtchanges=temp1-temp2;
        temp1_seq2=dataset.event_RT_trial1InSeq{i};
        temp1_seq2=temp1_seq2(dataset.realrtpair_seq2{i}==1 & ismember(metadata.sessid(dataset.event_isSeq{i}==1),whichsess)');
        temp2_seq2=dataset.event_RT_trialiInSeq{i};
        temp2_seq2=temp2_seq2(dataset.realrtpair_seq2{i}==1 & ismember(metadata.sessid(dataset.event_isSeq{i}==1),whichsess)');
        rtchanges_seq2=temp1_seq2-temp2_seq2;
        plotScatter(temp1,rtchanges,temp1_seq2,rtchanges_seq2,scatterJitter,scatterJitter,['Change in RT ' num2str(dataset.nInSequence(i)-1) ' trials later as a function of RT'],'RT trials 1','Change in RT',alpha);
    end
end

if plot_delta_rt_asFunc_rt_removeMeanRegression==true
    for i=1:1
        temp1=dataset.allTrialsSequence_RT_trial1InSeq{i};
        temp1=temp1(dataset.realrtpair_seq1{i}==1 & ismember(metadata.sessid(dataset.allTrialsSequence_isSeq{i}==1),whichsess)');
        temp2=dataset.allTrialsSequence_RT_trialiInSeq{i};
        temp2=temp2(dataset.realrtpair_seq1{i}==1 & ismember(metadata.sessid(dataset.allTrialsSequence_isSeq{i}==1),whichsess)');
        rtchanges=temp1-temp2;
        temp1_seq2=dataset.event_RT_trial1InSeq{i};
        temp1_seq2=temp1_seq2(dataset.realrtpair_seq2{i}==1 & ismember(metadata.sessid(dataset.event_isSeq{i}==1),whichsess)');
        temp2_seq2=dataset.event_RT_trialiInSeq{i};
        temp2_seq2=temp2_seq2(dataset.realrtpair_seq2{i}==1 & ismember(metadata.sessid(dataset.event_isSeq{i}==1),whichsess)');
        rtchanges_seq2=temp1_seq2-temp2_seq2;
        % do all trials first
        [diff2Dhist_alltrials,x,y]=removeMeanRegression(temp1,temp1,rtchanges,scatterJitter,alpha,'All trials',nBinsFor2Dhist);
        % do event trials second
        %removeMeanRegression(temp1_seq2,temp1_seq2,rtchanges_seq2,scatterJitter,alpha,'Event',nBinsFor2Dhist); % get bootstrapped distribution from event pdf
        [diff2Dhist_event,x,y]=removeMeanRegression(temp1,temp1_seq2,rtchanges_seq2,scatterJitter,alpha,'Event',{x; y}); % get bootstrapped distribution from all trials pdf
        figure();
        imagesc(x,y,diff2Dhist_event'-diff2Dhist_alltrials');
        set(gca,'YDir','normal');
        title(['Difference of all trials and event histograms']);
        xlabel('Reaction time trial 1');
        ylabel('Change in reaction times');
        if smoothSize>1
            K=ones(smoothSize);
            smoothMat=conv2(diff2Dhist_event'-diff2Dhist_alltrials',K,'same');
            figure();
            imagesc(x,y,smoothMat);
            set(gca,'YDir','normal');
            title(['Smoothed -- Difference of all trials and event histograms']);
            xlabel('Reaction time trial 1');
            ylabel('Change in reaction times');
        end
    end
end

if plot_dim2_delta_asFunc_rt==true
    for i=1:length(dataset.allTrialsSequence_RT_trial1InSeq)
        temp1=dataset.allTrialsSequence_RT_trial1InSeq{i};
        temp1=temp1(dataset.realrtpair_seq1{i}==1 & ismember(metadata.sessid(dataset.allTrialsSequence_isSeq{i}==1),whichsess)');
        temp2=dataset.dim2_rtchanges_allTrialsSequence{i};
        isThisSess=ismember(metadata.sessid(dataset.allTrialsSequence_isSeq{i}==1),whichsess)';
        temp2=temp2(isThisSess(dataset.realrtpair_seq1{i}==1));
        
        temp1_seq2=dataset.event_RT_trial1InSeq{i};
        temp1_seq2=temp1_seq2(dataset.realrtpair_seq2{i}==1 & ismember(metadata.sessid(dataset.event_isSeq{i}==1),whichsess)');
        temp2_seq2=dataset.dim2_rtchanges_event{i};
        isThisSess=ismember(metadata.sessid(dataset.event_isSeq{i}==1),whichsess)';
        temp2_seq2=temp2_seq2(isThisSess(dataset.realrtpair_seq2{i}==1));
        plotScatter(temp1,temp2,temp1_seq2,temp2_seq2,scatterJitter,scatterJitter,['Dim 2 of change in RT ' num2str(dataset.nInSequence(i)-1) ' trials later as a function of RT'],'RT trials 1','Change in RT',alpha);
    end 
end


end

function [diff2Dhist,x,y]=removeMeanRegression(rts,actual_first_rts,actual_rt_changes,scatterJitter,alpha,tit,nBinsFor2Dhist)

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
    curr_rts=rts(randsample(length(rts),2,true)); % with replacement
    delta_rts(i)=diff(curr_rts(end:-1:1));
    first_rt(i)=curr_rts(1);
    second_rt(i)=curr_rts(2);
end

% compare actual change in reaction times distribution to this bootstrapped
% distribution
% Scatter plot comparison
plotScatter(first_rt,delta_rts,actual_first_rts,actual_rt_changes,scatterJitter,scatterJitter,[tit ' comparing bootstrapped vs real rt pairs'],'RT trial 1','Change in RT',[alpha/10 alpha]);
% Heatmap comparison
[diff2Dhist,x,y]=compareWithHeatmaps(first_rt,delta_rts,actual_first_rts,actual_rt_changes,nBinsFor2Dhist,tit);

end

function [diff_n,x,y]=compareWithHeatmaps(x1,y1,x2,y2,nBinsPerDim,tit)

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
imagesc(bin_c{1},bin_c{2},n');
set(gca,'YDir','normal');
title([tit ' Bootstrap histogram']);
xlabel('Reaction time trial 1');
ylabel('Change in reaction times');

figure();
imagesc(bin_c{1},bin_c{2},n2');
set(gca,'YDir','normal');
title([tit ' Real events histogram']);
xlabel('Reaction time trial 1');
ylabel('Change in reaction times');

figure();
imagesc(bin_c{1},bin_c{2},diff_n');
set(gca,'YDir','normal');
title([tit ' Difference histogram']);
xlabel('Reaction time trial 1');
ylabel('Change in reaction times');

x=bin_c{1};
y=bin_c{2};


end

function plotScatter(RT_pairs1_x,RT_pairs1_y,RT_pairs2_x,RT_pairs2_y,jitter1,jitter2,tit,xlab,ylab,alpha)

if length(alpha)>1
    alpha2=alpha(2);
    alpha=alpha(1);
else
    alpha2=alpha;
end

figure();
if length(jitter1)==1
    useX=RT_pairs1_x+rand(size(RT_pairs1_x)).*jitter1;
    useY=RT_pairs1_y+rand(size(RT_pairs1_y)).*jitter1;
else
    useX=RT_pairs1_x+jitter1;
    useY=RT_pairs1_y+jitter1;
end
    
s=scatter(useX,useY,100,'k','filled');
set(gcf,'position',[10,10,1000,1000]);
pause;
m=get(s,'MarkerHandle');
%alpha=0.01;
m.FaceColorData=uint8(255*[0;0;0;alpha]);
xlabel(xlab);
ylabel(ylab);
title(tit);
% xlim([0 9.5]);
% ylim([0 9.5]);

hold on;
if length(jitter2)==1
    useX=RT_pairs2_x+rand(size(RT_pairs2_x)).*jitter2;
    useY=RT_pairs2_y+rand(size(RT_pairs2_y)).*jitter2;
else
    useX=RT_pairs2_x+jitter2;
    useY=RT_pairs2_y+jitter2;
end
s=scatter(useX,useY,100,'r','filled');
pause;
m=get(s,'MarkerHandle');
%alpha=0.01;
m.FaceColorData=uint8(255*[1;0;0;alpha2]);

end

function [x_backup]=plotHist(data1,data2,bins,tit,xlab)

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

function [x_backup]=plotCDF(data1,data2,bins,tit)

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