function plotBehaviorEventFx(dataset,alltbt)

plot_rawReaching=false;
plot_rt_pdf=true;
plot_rt_cdf=false;
plot_delta_rt_pdf=false;
plot_delta_rt_cdf=false;
plot_delta_rt_asFunc_rt=false;
plot_dim1_delta_asFunc_rt=false;
plot_dim2_delta_asFunc_rt=false;
plot_3D_dim1_dim2_asFunc_rt=false;
plot_delta_rt_asFunc_rt_removeMeanRegression=false;

histo_nbins=200; % number of bins for reaction time histogram
backup_histo_nbins=histo_nbins;
scatterJitter=0.5;
% alpha=0.01;
alpha=0.1;
%spotSize=1;
spotSize=10;
nBinsFor2Dhist=100;

% dataset contains ...
% Get raw reaching (all reaches)
% Get reaction times
% Get change in reaction times (non-corrected input distributions)
% Get change in reaction times (corrected input distributions)

% Plot raw reaching data
if plot_rawReaching==true
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
end

% Plot reaction times distribution
if plot_rt_pdf==true
    for i=1:length(dataset.allTrialsSequence_RT_trial1InSeq)
        histo_nbins=plotHist(dataset.allTrialsSequence_RT_trial1InSeq{i},dataset.allTrialsSequence_RT_trialiInSeq{i},histo_nbins,['Reaction Times all trials reference: trial 1 (black) vs ' num2str(dataset.nInSequence(i)-1) ' later (magenta)'],'RT (sec)');
    end
    temp=dataset.event_name;
    temp(regexp(temp,'_'))=' ';
    for i=1:length(dataset.event_RT_trial1InSeq)
        histo_nbins=plotHist(dataset.event_RT_trial1InSeq{i},dataset.event_RT_trialiInSeq{i},histo_nbins,['Reaction Times fx of ' temp ': trial 1 (black) vs ' num2str(dataset.nInSequence(i)-1) ' later (magenta)'],'RT (sec)');
    end
end

% Plot reaction times CDF
if plot_rt_cdf==true
    histo_nbins=backup_histo_nbins;
    for i=1:length(dataset.allTrialsSequence_RT_trial1InSeq)
        histo_nbins=plotCDF(dataset.allTrialsSequence_RT_trial1InSeq{i},dataset.allTrialsSequence_RT_trialiInSeq{i},histo_nbins,['CDF Reaction Times all trials reference: trial 1 (black) vs ' num2str(dataset.nInSequence(i)-1) ' later (magenta)']);
    end
    temp=dataset.event_name;
    temp(regexp(temp,'_'))=' ';
    for i=1:length(dataset.event_RT_trial1InSeq)
        histo_nbins=plotCDF(dataset.event_RT_trial1InSeq{i},dataset.event_RT_trialiInSeq{i},histo_nbins,['CDF Reaction Times fx of ' temp ': trial 1 (black) vs ' num2str(dataset.nInSequence(i)-1) ' later (magenta)']);
    end
end

% Plot change in reaction times
if plot_delta_rt_pdf==true
    histo_nbins=backup_histo_nbins;
    for i=1:length(dataset.allTrialsSequence_RT_trial1InSeq)
        histo_nbins=plotHist(dataset.dim1_rtchanges_allTrialsSequence{i},dataset.dim1_rtchanges_event{i},histo_nbins,['Dim 1 of change in RT ' num2str(dataset.nInSequence(i)-1) ' trials later, comparing reference vs ' temp],'Change in RT (sec)');
    end
    temp=dataset.event_name;
    temp(regexp(temp,'_'))=' ';
    for i=1:length(dataset.event_RT_trial1InSeq)
        histo_nbins=plotHist(dataset.dim2_rtchanges_allTrialsSequence{i},dataset.dim2_rtchanges_event{i},histo_nbins,['Dim 2 of change in RT ' num2str(dataset.nInSequence(i)-1) ' trials later, comparing reference vs ' temp],'Change in RT (sec)');
    end
end

% Plot CDF
if plot_delta_rt_cdf==true
    histo_nbins=backup_histo_nbins;
    for i=1:length(dataset.allTrialsSequence_RT_trial1InSeq)
        histo_nbins=plotCDF(dataset.dim1_rtchanges_allTrialsSequence{i},dataset.dim1_rtchanges_event{i},histo_nbins,['CDF Dim 1 of change in RT ' num2str(dataset.nInSequence(i)-1) ' trials later, comparing reference vs ' temp]);
    end
    temp=dataset.event_name;
    temp(regexp(temp,'_'))=' ';
    for i=1:length(dataset.event_RT_trial1InSeq)
        histo_nbins=plotCDF(dataset.dim2_rtchanges_allTrialsSequence{i},dataset.dim2_rtchanges_event{i},histo_nbins,['CDF Dim 2 of change in RT ' num2str(dataset.nInSequence(i)-1) ' trials later, comparing reference vs ' temp]);
    end
end

% Plot change in RT as a function of RT
if plot_delta_rt_asFunc_rt==true
    for i=1:length(dataset.allTrialsSequence_RT_trial1InSeq)
        temp1=dataset.allTrialsSequence_RT_trial1InSeq{i};
        temp2=dataset.allTrialsSequence_RT_trialiInSeq{i};
        rtchanges=temp1(dataset.realrtpair_seq1{i}==1)-temp2(dataset.realrtpair_seq1{i}==1);
        temp1_seq2=dataset.event_RT_trial1InSeq{i};
        temp2_seq2=dataset.event_RT_trialiInSeq{i};
        rtchanges_seq2=temp1_seq2(dataset.realrtpair_seq2{i}==1)-temp2_seq2(dataset.realrtpair_seq2{i}==1);
        plotScatter(temp1(dataset.realrtpair_seq1{i}==1),dataset.alldim_rtchanges_allTrialsSequence{i},temp1_seq2(dataset.realrtpair_seq2{i}==1),dataset.alldim_rtchanges_event{i},scatterJitter,scatterJitter,['Change in RT ' num2str(dataset.nInSequence(i)-1) ' trials later as a function of RT'],'RT trials 1','Change in RT',alpha);
        %plotScatter(temp1(dataset.realrtpair_seq1{i}==1),rtchanges,temp1_seq2(dataset.realrtpair_seq2{i}==1),rtchanges_seq2,scatterJitter,scatterJitter,['Change in RT ' num2str(dataset.nInSequence(i)-1) ' trials later as a function of RT'],'RT trials 1','Change in RT',alpha);
        %plotScatter(temp1(dataset.realrtpair_seq1{i}==1),dataset.dim1_rtchanges_allTrialsSequence{i}+dataset.dim2_rtchanges_allTrialsSequence{i},temp1_seq2(dataset.realrtpair_seq2{i}==1),dataset.dim1_rtchanges_event{i}+dataset.dim2_rtchanges_event{i},scatterJitter,scatterJitter,['Dim1+2 change in RT ' num2str(dataset.nInSequence(i)-1) ' trials later as a function of RT'],'RT trials 1','Change in RT',alpha);
    end
end
if plot_3D_dim1_dim2_asFunc_rt==true
    for i=1:length(dataset.allTrialsSequence_RT_trial1InSeq)
        temp1=dataset.allTrialsSequence_RT_trial1InSeq{i};
        temp2=dataset.allTrialsSequence_RT_trialiInSeq{i};
        rtchanges=temp1(dataset.realrtpair_seq1{i}==1)-temp2(dataset.realrtpair_seq1{i}==1);
        temp1_seq2=dataset.event_RT_trial1InSeq{i};
        temp2_seq2=dataset.event_RT_trialiInSeq{i};
        rtchanges_seq2=temp1_seq2(dataset.realrtpair_seq2{i}==1)-temp2_seq2(dataset.realrtpair_seq2{i}==1);
        plotScatter3D(temp1(dataset.realrtpair_seq1{i}==1),dataset.dim1_rtchanges_allTrialsSequence{i},dataset.dim2_rtchanges_allTrialsSequence{i},temp1_seq2(dataset.realrtpair_seq2{i}==1),dataset.dim1_rtchanges_event{i},dataset.dim2_rtchanges_event{i},scatterJitter,'3D','RT trial 1','Dim 1 delta RT','Dim 2 delta RT',spotSize)
    end
end
if plot_dim1_delta_asFunc_rt==true
    for i=1:length(dataset.allTrialsSequence_RT_trial1InSeq)
        temp1=dataset.allTrialsSequence_RT_trial1InSeq{i};
        temp1_seq2=dataset.event_RT_trial1InSeq{i};
        plotScatter(temp1(dataset.realrtpair_seq1{i}==1),dataset.dim1_rtchanges_allTrialsSequence{i},temp1_seq2(dataset.realrtpair_seq2{i}==1),dataset.dim1_rtchanges_event{i},scatterJitter,scatterJitter,['Dim 1 of change in RT ' num2str(dataset.nInSequence(i)-1) ' trials later as a function of RT'],'RT trials 1','Change in RT',alpha);
    end 
    rtBins=[0 2; 2 4; 4 8; 8 15];
    histo_nbins=backup_histo_nbins;
    for i=1:length(dataset.allTrialsSequence_RT_trial1InSeq)
        for j=1:size(rtBins,1)
            temp1=dataset.allTrialsSequence_RT_trial1InSeq{i};
            temp1_seq2=dataset.event_RT_trial1InSeq{i};
            dim1_all=dataset.dim1_rtchanges_allTrialsSequence{i};
            dim1_event=dataset.dim1_rtchanges_event{i};
            histo_nbins=plotCDF(dim1_all(temp1(dataset.realrtpair_seq1{i}==1)>=rtBins(j,1) & temp1(dataset.realrtpair_seq1{i}==1)<rtBins(j,2)),dim1_event(temp1_seq2(dataset.realrtpair_seq2{i}==1)>=rtBins(j,1) & temp1_seq2(dataset.realrtpair_seq2{i}==1)<rtBins(j,2)),histo_nbins,['CDF RT change for RTs ' num2str(dataset.nInSequence(i)-1) ' trials later, RTs less than ' num2str(rtBins(j,2))]);
        end
    end
end
if plot_dim2_delta_asFunc_rt==true
    for i=1:length(dataset.allTrialsSequence_RT_trial1InSeq)
        temp1=dataset.allTrialsSequence_RT_trial1InSeq{i};
        temp1_seq2=dataset.event_RT_trial1InSeq{i};
        plotScatter(temp1(dataset.realrtpair_seq1{i}==1),dataset.dim2_rtchanges_allTrialsSequence{i},temp1_seq2(dataset.realrtpair_seq2{i}==1),dataset.dim2_rtchanges_event{i},scatterJitter,scatterJitter,['Dim 2 of change in RT ' num2str(dataset.nInSequence(i)-1) ' trials later as a function of RT'],'RT trials 1','Change in RT',alpha);
    end 
    rtBins=[0 2; 2 4; 4 8; 8 15];
    histo_nbins=backup_histo_nbins;
    for i=1:length(dataset.allTrialsSequence_RT_trial1InSeq)
        for j=1:size(rtBins,1)
            temp1=dataset.allTrialsSequence_RT_trial1InSeq{i};
            temp1_seq2=dataset.event_RT_trial1InSeq{i};
            dim2_all=dataset.dim2_rtchanges_allTrialsSequence{i};
            dim2_event=dataset.dim2_rtchanges_event{i};
            histo_nbins=plotCDF(dim2_all(temp1(dataset.realrtpair_seq1{i}==1)>=rtBins(j,1) & temp1(dataset.realrtpair_seq1{i}==1)<rtBins(j,2)),dim2_event(temp1_seq2(dataset.realrtpair_seq2{i}==1)>=rtBins(j,1) & temp1_seq2(dataset.realrtpair_seq2{i}==1)<rtBins(j,2)),histo_nbins,['CDF RT change for RTs ' num2str(dataset.nInSequence(i)-1) ' trials later, RTs less than ' num2str(rtBins(j,2))]);
        end
    end
end

if plot_delta_rt_asFunc_rt_removeMeanRegression==true
    for i=1:1
        temp1=dataset.allTrialsSequence_RT_trial1InSeq{i};
        temp1=temp1(dataset.realrtpair_seq1{i}==1);
        temp2=dataset.allTrialsSequence_RT_trialiInSeq{i};
        temp2=temp2(dataset.realrtpair_seq1{i}==1);
        rtchanges=temp1-temp2;
        temp1_seq2=dataset.event_RT_trial1InSeq{i};
        temp1_seq2=temp1_seq2(dataset.realrtpair_seq2{i}==1);
        temp2_seq2=dataset.event_RT_trialiInSeq{i};
        temp2_seq2=temp2_seq2(dataset.realrtpair_seq2{i}==1);
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
plotScatter(first_rt,delta_rts,actual_first_rts,actual_rt_changes,scatterJitter,scatterJitter,[tit ' comparing bootstrapped vs real rt pairs'],'RT trial 1','Change in RT',alpha);
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

function plotScatter(RT_pairs1_x,RT_pairs1_y,RT_pairs2_x,RT_pairs2_y,jitter1,jitter2,tit,xlab,ylab,alpha)

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
m.FaceColorData=uint8(255*[1;0;0;alpha]);

end

function plotScatter3D(RT_pairs1_x,RT_pairs1_y,RT_pairs1_z,RT_pairs2_x,RT_pairs2_y,RT_pairs2_z,jitter1,tit,xlab,ylab,zlab,spotSize)

figure();
useX=RT_pairs1_x+rand(size(RT_pairs1_x)).*jitter1;
useY=RT_pairs1_y+rand(size(RT_pairs1_y)).*jitter1;
useZ=RT_pairs1_z+rand(size(RT_pairs1_z)).*jitter1;
    
s=scatter3(useX,useY,useZ,spotSize,'k','filled');
set(gcf,'position',[10,10,1000,1000]);
% pause;
% m=get(s,'MarkerHandle');
% alpha=0.1;
% m.FaceColorData=uint8(255*[0;0;0;alpha]);
xlabel(xlab);
ylabel(ylab);
zlabel(zlab);
title(tit);
% xlim([0 9.5]);
% ylim([0 9.5]);

hold on;
useX=RT_pairs2_x+rand(size(RT_pairs2_x)).*jitter1;
useY=RT_pairs2_y+rand(size(RT_pairs2_y)).*jitter1;
useZ=RT_pairs2_z+rand(size(RT_pairs2_z)).*jitter1;

s=scatter3(useX,useY,useZ,spotSize,'r','filled');
% pause;
% m=get(s,'MarkerHandle');
% alpha=0.1;
% m.FaceColorData=uint8(255*[1;0;0;alpha]);

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