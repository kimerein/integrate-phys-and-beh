function [returnThis,returnThisRef,returnThisEventTrial1]=plotBehaviorEventFx(varargin)

% Dim 2 conversion is 0.7885 sec in dim 2 is 1 sec real
% Dim 1 conversion is 0.611 sec in dim 1 is 1 sec real

if length(varargin)==3
    dataset=varargin{1};
    alltbt=varargin{2};
    ref=varargin{3};
    typeOfPlot=[];
    suppressPlotBehPlots=[];
elseif length(varargin)==4
    dataset=varargin{1};
    alltbt=varargin{2};
    ref=varargin{3};
    typeOfPlot=varargin{4}; % indicates the type of plot
    suppressPlotBehPlots=[];
elseif length(varargin)==5
    dataset=varargin{1};
    alltbt=varargin{2};
    ref=varargin{3};
    typeOfPlot=varargin{4}; % indicates the type of plot
    suppressPlotBehPlots=varargin{5};
end

returnThis=[];
returnThisRef=[];
returnThisEventTrial1=[];

if isempty(suppressPlotBehPlots)
    whetherToSuppressPlots(false);
else
    whetherToSuppressPlots(suppressPlotBehPlots);
end

issupp=whetherToSuppressPlots();
if issupp==true
    disp('Suppressing all plots in plotBehaviorEventFx.m.');
else
end

if isempty(typeOfPlot)
    % manually specify plot type here
    plot_rawReaching=false;
    plot_rawReaching_cdf=false;
    plot_rt_pdf=false; 
    plot_rt_cdf=true;
    reverse_rt_cdf=false;
    plot_delta_rt_pdf=false;
    plot_delta_rt_pdf_2D=false;
    plot_delta_rt_cdf_2D=false;
    plot_delta_rt_cdf=false; 
    plot_delta_rt_asFunc_rt=false;
    plot_dim1_delta_asFunc_rt=false;
    plot_dim2_delta_asFunc_rt=false;
    plot_3D_dim1_dim2_asFunc_rt=false;
    plot_delta_rt_asFunc_rt_removeMeanRegression=false;
    plot_earth_mover_wander=false;
else
    plot_rawReaching=false;
    plot_rawReaching_cdf=false;
    plot_rt_pdf=false; 
    plot_rt_cdf=false;
    reverse_rt_cdf=false;
    plot_delta_rt_pdf=false;
    plot_delta_rt_pdf_2D=false;
    plot_delta_rt_cdf_2D=false;
    plot_delta_rt_cdf=false;
    plot_delta_rt_asFunc_rt=false;
    plot_dim1_delta_asFunc_rt=false;
    plot_dim2_delta_asFunc_rt=false;
    plot_3D_dim1_dim2_asFunc_rt=false;
    plot_delta_rt_asFunc_rt_removeMeanRegression=false;
    plot_earth_mover_wander=false;
    
    switch typeOfPlot
        case 'plot_rawReaching'
            plot_rawReaching=true;
        case 'plot_rawReaching_cdf'
            plot_rawReaching_cdf=true;
        case 'plot_rt_pdf'
            plot_rt_pdf=true;
        case 'plot_rt_cdf'
            plot_rt_cdf=true;
        case 'reverse_rt_cdf'
            reverse_rt_cdf=true;
        case 'plot_delta_rt_pdf'
            plot_delta_rt_pdf=true;
        case 'plot_delta_rt_pdf_2D'
            plot_delta_rt_pdf_2D=true;
        case 'plot_delta_rt_cdf_2D'
            plot_delta_rt_cdf_2D=true;
        case 'plot_delta_rt_cdf'
            plot_delta_rt_cdf=true;
        case 'plot_delta_rt_asFunc_rt'
            plot_delta_rt_asFunc_rt=true;
        case 'plot_dim1_delta_asFunc_rt'
            plot_dim1_delta_asFunc_rt=true;
        case 'plot_dim2_delta_asFunc_rt'
            plot_dim2_delta_asFunc_rt=true;
        case 'plot_3D_dim1_dim2_asFunc_rt'
            plot_3D_dim1_dim2_asFunc_rt=true;
        case 'plot_delta_rt_asFunc_rt_removeMeanRegression'
            plot_delta_rt_asFunc_rt_removeMeanRegression=true;
        case 'plot_earth_mover_wander'
            plot_earth_mover_wander=true;
        otherwise 
            error('Do not recognize typeOfPlot parameter value');
    end
end

% histo_nbins=200; % number of bins for reaction time histogram
histo_nbins=[-4*12.4245:0.2510:4*12.4245];
% histo_nbins=[-4*12.4245-2*0.2510:4*0.2510:4*12.4245];
% histo_nbins=[-4*12.4245-0.75*0.2510:1.5*0.2510:4*12.4245];
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

% Plot raw reaching CDF
if plot_rawReaching_cdf==true
    timeStep=mode(diff(nanmean(alltbt.times,1)));
    timeBinsForReaching=0:timeStep:(size(dataset.rawReaching_allTrialsSequence_trial1InSeq{1},2)-1)*timeStep;
    % find cue ind
    [~,f]=nanmax(nanmean(alltbt.cueZone_onVoff,1));
    cueTime=timeBinsForReaching(f);
    for i=1:length(dataset.allTrialsSequence_RT_trial1InSeq)
        plotCDF_rawReaches(dataset.rawReaching_allTrialsSequence_trial1InSeq{i},dataset.rawReaching_allTrialsSequence_trialiInSeq{i},timeBinsForReaching,cueTime,['CDF Raw Reaches all trials reference: trial 1 (black) vs ' num2str(dataset.nInSequence(i)-1) ' later (red)']);
    end
    temp=dataset.event_name;
    temp(regexp(temp,'_'))=' ';
    for i=1:length(dataset.event_RT_trial1InSeq)
        plotCDF_rawReaches(dataset.rawReaching_event_trial1InSeq{i},dataset.rawReaching_event_trialiInSeq{i},timeBinsForReaching,cueTime,['CDF Raw Reaches fx of ' temp ': trial 1 (black) vs ' num2str(dataset.nInSequence(i)-1) ' later (red)']);
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
    if reverse_rt_cdf
        offset=14;
        multiple=-1;
    else
        offset=0;
        multiple=1;
    end
    histo_nbins=backup_histo_nbins;
    for i=1:length(dataset.allTrialsSequence_RT_trial1InSeq)
        [histo_nbins,returnThisRef]=plotCDF(offset+multiple*dataset.allTrialsSequence_RT_trial1InSeq{i},offset+multiple*dataset.allTrialsSequence_RT_trialiInSeq{i},histo_nbins,['CDF Reaction Times all trials reference: trial 1 (black) vs ' num2str(dataset.nInSequence(i)-1) ' later (magenta)']);
    end
    temp=dataset.event_name;
    temp(regexp(temp,'_'))=' ';
    for i=1:length(dataset.event_RT_trial1InSeq)
        [histo_nbins,returnThis,returnThisEventTrial1]=plotCDF(offset+multiple*dataset.event_RT_trial1InSeq{i},offset+multiple*dataset.event_RT_trialiInSeq{i},histo_nbins,['CDF Reaction Times fx of ' temp ': trial 1 (black) vs ' num2str(dataset.nInSequence(i)-1) ' later (magenta)']);
    end
end

% Plot change in reaction times
% previous minus current from line 80 of getPairedReactionTimes.m
if plot_delta_rt_pdf==true
    temp=dataset.event_name;
    temp(regexp(temp,'_'))=' ';
    histo_nbins=backup_histo_nbins;
    for i=1:length(dataset.allTrialsSequence_RT_trial1InSeq)
        histo_nbins=plotHist(dataset.alldim_rtchanges_allTrialsSequence{i},dataset.alldim_rtchanges_event{i},histo_nbins,['All dim of change in RT ' num2str(dataset.nInSequence(i)-1) ' trials later, comparing reference vs ' temp],'Change in RT (sec) previous minus current');
    end
    for i=1:length(dataset.allTrialsSequence_RT_trial1InSeq)
        histo_nbins=plotHist(dataset.dim1_rtchanges_allTrialsSequence{i},dataset.dim1_rtchanges_event{i},histo_nbins,['Dim 1 of change in RT ' num2str(dataset.nInSequence(i)-1) ' trials later, comparing reference vs ' temp],'Change in RT (sec) previous minus current');
    end
    for i=1:length(dataset.event_RT_trial1InSeq)
        [histo_nbins,returnThis_temp]=plotHist(dataset.dim2_rtchanges_allTrialsSequence{i},dataset.dim2_rtchanges_event{i},histo_nbins,['Dim 2 of change in RT ' num2str(dataset.nInSequence(i)-1) ' trials later, comparing reference vs ' temp],'Change in RT (sec) previous minus current');
    end
end

% Plot change in reaction times 2D smear
if plot_delta_rt_pdf_2D==true
    temp=dataset.event_name;
    temp(regexp(temp,'_'))=' ';
    histo_nbins=backup_histo_nbins;
    for i=1:length(dataset.allTrialsSequence_RT_trial1InSeq)
        plotHist_2Dsmear(dataset.dim1_rtchanges_event{i},dataset.dim2_rtchanges_event{i},histo_nbins,'2D smear','Change in RT (sec)');
    end
end

% Plot change in reaction times 2D smear CDF
if plot_delta_rt_cdf_2D==true
    temp=dataset.event_name;
    temp(regexp(temp,'_'))=' ';
    histo_nbins=backup_histo_nbins;
    for i=1:length(dataset.allTrialsSequence_RT_trial1InSeq)
        plotCDF_2Dsmear(dataset.dim1_rtchanges_allTrialsSequence{i},dataset.dim2_rtchanges_allTrialsSequence{i},dataset.dim1_rtchanges_event{i},dataset.dim2_rtchanges_event{i},histo_nbins,'2D smear','Change in RT (sec)');
    end
end

% Plot wander in dim 1 and 2
if plot_earth_mover_wander==true
    temp=dataset.event_name;
    temp(regexp(temp,'_'))=' ';
    histo_nbins=backup_histo_nbins;
    for i=1:length(dataset.allTrialsSequence_RT_trial1InSeq)
%     for i=1:1
        if isempty(ref)
            plot_earthmover_wander(dataset.dim1_rtchanges_allTrialsSequence{i},dataset.dim2_rtchanges_allTrialsSequence{i},dataset.dim1_rtchanges_event{i},dataset.dim2_rtchanges_event{i},histo_nbins);
        else
            plot_earthmover_wander(ref.dim1,ref.dim2,dataset.dim1_rtchanges_event{i},dataset.dim2_rtchanges_event{i},histo_nbins);
        end
    end
end

% Plot CDF
if plot_delta_rt_cdf==true
    histo_nbins=backup_histo_nbins;
    temp=dataset.event_name;
    temp(regexp(temp,'_'))=' ';
%     for i=1:length(dataset.allTrialsSequence_RT_trial1InSeq)
% %         histo_nbins=plotCDF(dataset.alldim_rtchanges_allTrialsSequence{i},dataset.alldim_rtchanges_event{i},histo_nbins,['CDF All Dim change in RT ' num2str(dataset.nInSequence(i)-1) ' trials later, comparing reference vs ' temp]);
%         if i==1
%             
%             rtBins=[6 9];
%             j=1;
%             temp1=dataset.allTrialsSequence_RT_trial1InSeq{i};
%             temp1_seq2=dataset.event_RT_trial1InSeq{i};
%             dimall_all=dataset.alldim_rtchanges_allTrialsSequence{i};
%             dimall_event=dataset.alldim_rtchanges_event{i};
%             temp_RTdiffs1=dimall_all(temp1(dataset.realrtpair_seq1{i}==1)>=rtBins(j,1) & temp1(dataset.realrtpair_seq1{i}==1)<rtBins(j,2));
%             temp_RTdiffsevent=dimall_event(temp1_seq2(dataset.realrtpair_seq2{i}==1)>=rtBins(j,1) & temp1_seq2(dataset.realrtpair_seq2{i}==1)<rtBins(j,2));
%             temp_RTdiffs1=temp_RTdiffs1(temp_RTdiffs1<=-5);
%             temp_RTdiffsevent=temp_RTdiffsevent(temp_RTdiffsevent<=-5);
%             histo_nbins=plotHist(temp_RTdiffs1,temp_RTdiffsevent,histo_nbins,['PDF RT change for RTs ' num2str(dataset.nInSequence(i)-1) ' trials later, RTs less than ' num2str(rtBins(j,2))],'Change in RT');
%             plotCDF(temp_RTdiffs1,temp_RTdiffsevent,histo_nbins,['CDF RT change for RTs ' num2str(dataset.nInSequence(i)-1) ' trials later, RTs less than ' num2str(rtBins(j,2))]);
%         end
%     end
    for i=1:length(dataset.allTrialsSequence_RT_trial1InSeq)
        [histo_nbins,returnThis_temp,returnThisRef_temp]=plotCDF(dataset.alldim_rtchanges_allTrialsSequence{i},dataset.alldim_rtchanges_event{i},histo_nbins,['CDF All Dim of change in RT ' num2str(dataset.nInSequence(i)-1) ' trials later, comparing reference vs ' temp]);
        if i==1
            returnThis=returnThis_temp;
            returnThisRef=returnThisRef_temp;
        end
    end
%     for i=1:length(dataset.allTrialsSequence_RT_trial1InSeq)
%         histo_nbins=plotCDF(dataset.dim1_rtchanges_allTrialsSequence{i},dataset.dim1_rtchanges_event{i},histo_nbins,['CDF Dim 1 of change in RT ' num2str(dataset.nInSequence(i)-1) ' trials later, comparing reference vs ' temp]);
%     end
%     for i=1:length(dataset.event_RT_trial1InSeq)
%         [histo_nbins]=plotCDF(dataset.dim2_rtchanges_allTrialsSequence{i},dataset.dim2_rtchanges_event{i},histo_nbins,['CDF Dim 2 of change in RT ' num2str(dataset.nInSequence(i)-1) ' trials later, comparing reference vs ' temp]);
%     end
end

% Plot change in RT as a function of RT
if plot_delta_rt_asFunc_rt==true
    for i=1:length(dataset.allTrialsSequence_RT_trial1InSeq)
%      for i=1:1
        temp1=dataset.allTrialsSequence_RT_trial1InSeq{i};
        temp2=dataset.allTrialsSequence_RT_trialiInSeq{i};
        rtchanges=temp1(dataset.realrtpair_seq1{i}==1)-temp2(dataset.realrtpair_seq1{i}==1);
        temp1_seq2=dataset.event_RT_trial1InSeq{i};
        temp2_seq2=dataset.event_RT_trialiInSeq{i};
        rtchanges_seq2=temp1_seq2(dataset.realrtpair_seq2{i}==1)-temp2_seq2(dataset.realrtpair_seq2{i}==1);
%         plotScatter(temp1(dataset.realrtpair_seq1{i}==1),dataset.alldim_rtchanges_allTrialsSequence{i},temp1_seq2(dataset.realrtpair_seq2{i}==1),dataset.alldim_rtchanges_event{i},scatterJitter,scatterJitter,['Change in RT ' num2str(dataset.nInSequence(i)-1) ' trials later as a function of RT'],'RT trials 1','Change in RT',alpha);
        %plotScatter(temp1(dataset.realrtpair_seq1{i}==1),rtchanges,temp1_seq2(dataset.realrtpair_seq2{i}==1),rtchanges_seq2,scatterJitter,scatterJitter,['Change in RT ' num2str(dataset.nInSequence(i)-1) ' trials later as a function of RT'],'RT trials 1','Change in RT',alpha);
        %plotScatter(temp1(dataset.realrtpair_seq1{i}==1),dataset.dim1_rtchanges_allTrialsSequence{i}+dataset.dim2_rtchanges_allTrialsSequence{i},temp1_seq2(dataset.realrtpair_seq2{i}==1),dataset.dim1_rtchanges_event{i}+dataset.dim2_rtchanges_event{i},scatterJitter,scatterJitter,['Dim1+2 change in RT ' num2str(dataset.nInSequence(i)-1) ' trials later as a function of RT'],'RT trials 1','Change in RT',alpha);
        if i==1
            figure();
            [diff_n,x,y,~,n]=compareWithHeatmaps(temp1(dataset.realrtpair_seq1{i}==1),dataset.alldim_rtchanges_allTrialsSequence{i},temp1_seq2(dataset.realrtpair_seq2{i}==1),dataset.alldim_rtchanges_event{i},{histo_nbins; histo_nbins},'Difference Ref v Event');
%             for j=1:size(diff_n,2)
%                 diff_n(:,j)=smooth(diff_n(:,j),5);
%             end
            returnThis=n;
            for j=1:size(n,1)
                n(j,:)=smooth(n(j,:),5);
            end
            temp=n;
            temp=temp(x>0 & x<2,y>-9 & y<4);
            temp=temp./nansum(temp(1:end));
            nSmooth=1;
            K=ones(nSmooth);
%             temp=conv2(log(diff_n-nanmin(diff_n(1:end))),K,'same');
            temp=conv2(log(n),K,'same');
            temp=temp(nSmooth:end-nSmooth,nSmooth:end-nSmooth);
            imagesc(x(nSmooth:end-nSmooth),y(nSmooth:end-nSmooth),temp');
            set(gca,'YDir','normal');
            title(['Heatmap event diff']);
            xlabel('Reaction time trial 1');
            ylabel('Change in reaction times');
            
            figure();
            imagesc(x,y,log(diff_n-nanmin(nanmin(diff_n)))');
            set(gca,'YDir','normal');

        end
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
%     rtBins=[0 2; 2 4; 4 8; 8 15];
    rtBins=[0 15];
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

function out=whetherToSuppressPlots(varargin)

persistent suppressPlots

if length(varargin)==1
    suppressPlots=varargin{1};
end

out=suppressPlots;

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

function [diff_n,x,y,n,n2]=compareWithHeatmaps(x1,y1,x2,y2,nBinsPerDim,tit)

suppPlots=whetherToSuppressPlots();

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
if suppPlots==false
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
end

x=bin_c{1};
y=bin_c{2};


end

function plotTimeseries(data1_mean,data1_se,color1,data2_mean,data2_se,color2,timeBins)

suppPlots=whetherToSuppressPlots();

plotAsCityscape=true;

data1_mean=nanmean(data1_mean,1);
data2_mean=nanmean(data2_mean,1);
data1_se=sqrt(nansum(data1_se.^2,1));
data2_se=sqrt(nansum(data2_se.^2,1));

if suppPlots==false
    figure();
end
%fill([timeBins fliplr(timeBins)],[data1_mean+data1_se fliplr(data1_mean-data1_se)],[0.5 0.5 0.5]);
%hold on;
if plotAsCityscape==true
    if suppPlots==false
        [n,x]=cityscape_hist(data1_mean,timeBins);
        plot(x,n,'Color',color1); hold on;
        [n,x]=cityscape_hist(data1_mean+data1_se,timeBins);
        plot(x,n,'Color',color1);
        [n,x]=cityscape_hist(data1_mean-data1_se,timeBins);
        plot(x,n,'Color',color1);
    end
else
    if suppPlots==false    
        plot(timeBins,data1_mean,'Color',color1); hold on;
        plot(timeBins,data1_mean+data1_se,'Color',color1);
        plot(timeBins,data1_mean-data1_se,'Color',color1);
    end
end

%fill([timeBins fliplr(timeBins)],[data2_mean+data2_se fliplr(data2_mean-data2_se)],[0.1 0.7 0.5]);
%hold on;
if plotAsCityscape==true
    if suppPlots==false   
        [n,x]=cityscape_hist(data2_mean,timeBins);
        plot(x,n,'Color',color2); hold on;
        [n,x]=cityscape_hist(data2_mean+data2_se,timeBins);
        plot(x,n,'Color',color2);
        [n,x]=cityscape_hist(data2_mean-data2_se,timeBins);
        plot(x,n,'Color',color2);
    end
else
    if suppPlots==false    
        plot(timeBins,data2_mean,'Color',color2); hold on;
        plot(timeBins,data2_mean+data2_se,'Color',color2);
        plot(timeBins,data2_mean-data2_se,'Color',color2);
    end
end

end

function plotScatter(RT_pairs1_x,RT_pairs1_y,RT_pairs2_x,RT_pairs2_y,jitter1,jitter2,tit,xlab,ylab,alpha)

suppPlots=whetherToSuppressPlots();

if suppPlots==false   
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

end

function plotScatter3D(RT_pairs1_x,RT_pairs1_y,RT_pairs1_z,RT_pairs2_x,RT_pairs2_y,RT_pairs2_z,jitter1,tit,xlab,ylab,zlab,spotSize)

suppPlots=whetherToSuppressPlots();

if suppPlots==false   
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
    disp('p-value of ranksum');
    disp(p);
end

end

function [x_backup,returnThis,returnThisRef]=plotCDF(data1,data2,bins,tit)

returnThis=[];

doKStest=true;

suppPlots=whetherToSuppressPlots();

[n,x]=histcounts(data1,bins);
x_backup=x;
x_mids=nanmean([x(1:end-1); x(2:end)],1);
cond1_cdf=accumulateDistribution(n);
if suppPlots==false
    figure();
    plot(x_mids,cond1_cdf./nanmax(cond1_cdf),'Color','k');
    xlabel('CDF');
    ylabel('Count');
    title(tit);
end
returnThisRef.x=x_mids;
returnThisRef.y=cond1_cdf./nanmax(cond1_cdf);

[n,x]=histcounts(data2,x_backup);
cond2_cdf=accumulateDistribution(n);
if suppPlots==false
    hold on;
    plot(x_mids,cond2_cdf./nanmax(cond2_cdf),'Color','r');
end
returnThis.x=x_mids;
returnThis.y=cond2_cdf./nanmax(cond2_cdf);

if suppPlots==false
    if doKStest==true
        [~,p]=kstest2(data1,data2);
        disp('kstest pval');
        disp(p);
    end
end

end

function plotCDF_rawReaches(data1,data2,timesteps,cueTime,tit)

% raw reaching data comes in as a timeseries (i.e., already a pdf)
% so simply accumulate distribution, selecting only time points after the
% cue

suppPlots=whetherToSuppressPlots();

if suppPlots==false
    cond1_cdf=accumulateDistribution(data1(timesteps>cueTime));
    figure();
    plot(timesteps(timesteps>cueTime),cond1_cdf./nanmax(cond1_cdf),'Color','k');
    xlabel('CDF');
    ylabel('Count');
    title(tit);
    
    hold on;
    cond2_cdf=accumulateDistribution(data2(timesteps>cueTime));
    plot(timesteps(timesteps>cueTime),cond2_cdf./nanmax(cond2_cdf),'Color','r');
end

end

function cdf=accumulateDistribution(data)

cdf=nan(size(data));
for i=1:length(data)
    cdf(i)=nansum(data(1:i));
end

end

function [x_backup,returnThis]=plotHist(data1,data2,bins,tit,xlab)

useLogForY=false;
suppPlots=whetherToSuppressPlots();

[n,x]=histcounts(data1,bins);
x_backup=x;
[n,x]=cityscape_hist(n,x);
if suppPlots==false
    figure();
    if useLogForY==true
        semilogy(x,n./nansum(n),'Color','k');
    else
        plot(x,n./nansum(n),'Color','k');
    end
    xlabel(xlab);
    ylabel('Count');
    title(tit);
end

[n,x]=histcounts(data2,x_backup);
[n,x]=cityscape_hist(n,x);
if suppPlots==false
    hold on;
    if useLogForY==true
        semilogy(x,n./nansum(n),'Color','r');
    else
        plot(x,n./nansum(n),'Color','r');
    end
    leg={'data1','data2'};
    legend(leg);
end

returnThis.x=x;
returnThis.y=n./nansum(n);

end

function [x_backup]=plotHist_2Dsmear(data1_x,data1_y,bins,tit,xlab)

useLogForY=false;
suppPlots=whetherToSuppressPlots();

[n,x]=histcounts(data1_x,bins);
x_mids=nanmean([x(1:end-1); x(2:end)],1);
x_backup=x;

[n2,x]=histcounts(data1_y,x_backup);
if suppPlots==false
    figure();
    if useLogForY==true
        imagesc(x_mids,x_mids,log(repmat(n,length(x)-1,1)+repmat(n2',1,length(x)-1)));
    else
        imagesc(x_mids,x_mids,repmat(n,length(x)-1,1)+repmat(n2',1,length(x)-1));
    end
    set(gca,'YDir','normal');
end

end

function [x_backup]=plotCDF_2Dsmear(data1_x,data1_y,data2_x,data2_y,bins,tit,xlab)

useLogForY=false;
suppPlots=whetherToSuppressPlots();

[n,x]=histcounts(data1_x,bins);
x_mids=nanmean([x(1:end-1); x(2:end)],1);
x_cdf=accumulateDistribution(n);
x_cdf=x_cdf./nanmax(x_cdf);
x_backup=x;
[n,x]=histcounts(data2_x,bins);
x_cdf_cond2=accumulateDistribution(n);
x_cdf_cond2=x_cdf_cond2./nanmax(x_cdf_cond2);

[n2,x]=histcounts(data1_y,x_backup);
y_cdf=accumulateDistribution(n2);
y_cdf=y_cdf./nanmax(y_cdf);
[n,x]=histcounts(data2_y,bins);
y_cdf_cond2=accumulateDistribution(n);
y_cdf_cond2=y_cdf_cond2./nanmax(y_cdf_cond2);

temp=repmat(x_cdf_cond2-x_cdf,length(x)-1,1)+repmat(y_cdf'-y_cdf_cond2',1,length(x)-1);
temp=temp(x_mids>0,x_mids<0);

if suppPlots==false
    figure();
    if useLogForY==true
        imagesc(x_mids,x_mids,log(temp));
    else
        imagesc(x_mids,x_mids,temp);
    end
    set(gca,'YDir','normal');
end

end

function plot_earthmover_wander(data1_x,data1_y,data2_x,data2_y,bins)

doBootstrap=true;
suppPlots=whetherToSuppressPlots();

if doBootstrap==true
    nReps=100;
    earth_mover_dim1=nan(1,nReps);
    earth_mover_dim2=nan(1,nReps);
    for i=1:nReps
        curr_data1_x=data1_x(randsample(1:length(data1_x),length(data1_x),true));
        curr_data1_y=data1_y(randsample(1:length(data1_y),length(data1_y),true));
        curr_data2_x=data2_x(randsample(1:length(data2_x),length(data2_x),true));
        curr_data2_y=data2_y(randsample(1:length(data2_y),length(data2_y),true));
        [earth_mover_dim1(i),earth_mover_dim2(i)]=sub_earthmover_wander(curr_data1_x,curr_data1_y,curr_data2_x,curr_data2_y,bins);
    end
    em_dim1=prctile(earth_mover_dim1,[5 50 95]);
    em_dim2=prctile(earth_mover_dim2,[5 50 95]);
    if suppPlots==false
        figure();
        scatter(em_dim1(2),em_dim2(2));
        hold on;
        line([em_dim1(1) em_dim1(3)],[em_dim2(2) em_dim2(2)]);
        line([em_dim1(2) em_dim1(2)],[em_dim2(1) em_dim2(3)]);
    end
else
    [earth_mover_dim1,earth_mover_dim2]=sub_earthmover_wander(data1_x,data1_y,data2_x,data2_y,bins);
    if suppPlots==false
        figure();
        scatter(earth_mover_dim1,earth_mover_dim2);
    end
end

end

function [earth_mover_dim1,earth_mover_dim2]=sub_earthmover_wander(data1_x,data1_y,data2_x,data2_y,bins)
        
% Dim 2 conversion is 0.7885 sec in dim 2 is 1 sec real
% Dim 1 conversion is 0.611 sec in dim 1 is 1 sec real

x_binsize=1/0.611;
y_binsize=1/0.7885;

[n,x]=histcounts(data1_x,bins);
x_mids=nanmean([x(1:end-1); x(2:end)],1);
x_cdf=accumulateDistribution(n);
x_cdf=x_cdf./nanmax(x_cdf);
x_backup=x;
[n,x]=histcounts(data2_x,bins);
x_cdf_cond2=accumulateDistribution(n);
x_cdf_cond2=x_cdf_cond2./nanmax(x_cdf_cond2);

[n2,x]=histcounts(data1_y,x_backup);
y_cdf=accumulateDistribution(n2);
y_cdf=y_cdf./nanmax(y_cdf);
[n,x]=histcounts(data2_y,bins);
y_cdf_cond2=accumulateDistribution(n);
y_cdf_cond2=y_cdf_cond2./nanmax(y_cdf_cond2);

earth_mover_dim1=nansum(abs(x_cdf-x_cdf_cond2))*x_binsize;
sign_of_dim1=sign(nanmean(x_cdf-x_cdf_cond2));
earth_mover_dim1=sign_of_dim1*earth_mover_dim1;

earth_mover_dim2=nansum(abs(y_cdf-y_cdf_cond2))*y_binsize;
sign_of_dim2=sign(nanmean(y_cdf-y_cdf_cond2));
earth_mover_dim2=sign_of_dim2*earth_mover_dim2;

end