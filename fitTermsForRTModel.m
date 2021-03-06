function fitTermsForRTModel(varargin)

% this is an incarnation of the formal model
% 
% in this model, I use Matlab's parameterization of gamma distribution, 
% i.e., shape parameter a and scale parameter b

if length(varargin)==6
    dataset=varargin{1};
    wrtDataset=varargin{2};
    alltbt=varargin{3};
    reachType=varargin{4};
    useCuedRateData=varargin{5};
    realDataset=varargin{6};
    
    % Set up binning for histograms
    bins=0:0.1:9; % am choosing 0.1 sec as resolution, because that is the resolution of reaches from low-speed video
    histo_nbins=[-9:0.1:9];
elseif length(varargin)==7
    dataset=varargin{1};
    wrtDataset=varargin{2};
    alltbt=varargin{3};
    reachType=varargin{4};
    useCuedRateData=varargin{5};
    realDataset=varargin{6};
    settings=varargin{7};
    
    % Set up binning for histograms
    bins=settings.bins;
    histo_nbins=settings.histo_nbins;
    n_trials_away=settings.n_trials_away;
    trialStartCountAt1=settings.trialStartCountAt1;
    useAllDatasetTrialsAsRef=settings.useAllDatasetTrialsAsRef;
    cueName=settings.cueName;
    baselineWindow=settings.baselineWindow;
    n_sems=settings.n_sems;
    trialDuration=settings.trialDuration;
    nTrialsAheadForRT1=settings.nTrialsAheadForRT1;
end

% Plot real data
i=n_trials_away;
temp1_seq2=dataset.event_RT_trial1InSeq{i};
x2=temp1_seq2(dataset.realrtpair_seq2{i}==1);
y2=dataset.alldim_rtchanges_event{i};
nBinsPerDim={bins; histo_nbins};
realData=plotRTbyRTchangePlot(x2,y2,nBinsPerDim,'Real data');

% Plot reference data
% Reference input
if useAllDatasetTrialsAsRef==false
    temp1_seq2=wrtDataset.event_RT_trial1InSeq{i};
    x2=temp1_seq2(wrtDataset.realrtpair_seq2{i}==1);
    y2=wrtDataset.alldim_rtchanges_event{i};
else
    % All trials as reference
    temp1_seq2=dataset.allTrialsSequence_RT_trial1InSeq{i};
    x2=temp1_seq2(dataset.realrtpair_seq1{i}==1);
    y2=dataset.alldim_rtchanges_allTrialsSequence{i};
end
nBinsPerDim={bins; histo_nbins};
referenceData=plotRTbyRTchangePlot(x2,y2,nBinsPerDim,'Reference data');

% Plot reaction times
plotRTbasics(dataset,n_trials_away,nBinsPerDim{1});

% To find change in rate with behavioral event, take reaches far away from
% the cue, so that only contribution is change in Poisson random variable 
% (reaches) or exponential distribution change (RTs)
% Check for when this data set's reach PDF returns to baseline after cue
% Within n_sems SEMs of mean
timeInd=getReturnToBaseline(dataset,alltbt,cueName,n_sems,baselineWindow,n_trials_away,false);
withinRange=[timeInd histo_nbins(end)-3]; % in sec

% Real data rate term fit
rateFits=getRateFits(alltbt,dataset,reachType,useCuedRateData,withinRange,cueName,trialDuration,trialStartCountAt1);
checkRates(rateFits,dataset,nBinsPerDim,n_trials_away);
temp=makeRateTerm(nBinsPerDim,rateFits,dataset,n_trials_away);
prediction.rate_term=temp;

% Reference data rate term fit
rateFits=getRateFits(alltbt,wrtDataset,'miss',useCuedRateData,withinRange,cueName,trialDuration,trialStartCountAt1);
temp=makeRateTerm(nBinsPerDim,rateFits,wrtDataset,n_trials_away);
prediction.reference_rate_term=temp;

% RTM term
% temp=fitRTMtermToOtherReaches_local(dataset,nBinsPerDim,bins,1);
% prediction.rtm_term=temp;

% After removing component of pdf dependent on reaching rate, are left with
% component of pdf dependent on 
% Knowledge of relationship between cue and pellet

% Solve for gamma_cued
temp=solveForCuedReaching(dataset,realDataset,prediction.rate_term,nBinsPerDim,n_trials_away,alltbt,cueName);
prediction.cued_term=temp;

% Check ability of rate and cued terms combined together to predict real data
checkTermsTogether(dataset,prediction.rate_term,prediction.cued_term,nBinsPerDim,n_trials_away);

% Get cued reaching rate 
plotOnlyRTbins=[];
% plotOnlyRTbins=[0.075 4]; % in seconds; if plotting, include only the cued reaching for first trial RTs in this range; will plot all bins if is empty
prediction.cued_term=getRawCuedReachingFromRTpdf(nBinsPerDim,prediction.cued_term,dataset,n_trials_away,false,plotOnlyRTbins);
title('Cued reaching second trial in pair');

% Get cued reaching rate for first trial in pair
getRTsFromTbt=true;
firstTrial=firstTrial_solveForCuedReaching(dataset,prediction.rate_term,nBinsPerDim,n_trials_away,alltbt,cueName,baselineWindow,nTrialsAheadForRT1,realDataset,getRTsFromTbt);
firstTrial=getRawCuedReachingFromRTpdf(nBinsPerDim,firstTrial,dataset,n_trials_away,false,plotOnlyRTbins);
title('Cued reaching first trial in pair');

% Compare cued RTs across trial pair
compareCuedRTs(firstTrial,prediction.cued_term,dataset,n_trials_away,nBinsPerDim);

% Fit RPE
n_update_steps=trialStartCountAt1-1;
sameRTforEachTrial=false;
rpe_term=fitRPEtermFromData(dataset,nBinsPerDim,firstTrial,prediction.cued_term,n_trials_away,n_update_steps,sameRTforEachTrial);

return

% RPE term from theory
n_update_steps=1;
alpha=0.1;
temp=fitRPEtermFromTheory_local(dataset,nBinsPerDim,false,n_update_steps,alpha,true);
prediction.rpe_term=temp;

% RPE term by subtraction of wrtDataset
temp=fitRPEtermToEarlyReaches_local(dataset,wrtDataset,nBinsPerDim,1);
prediction.reference_rpe_term=temp;

plotRTbyRTchange(realData-prediction.rate_term.rt_change_pdfs,nBinsPerDim,'Real Data - Rate');
plotRTbyRTchange(referenceData-prediction.reference_rate_term.rt_change_pdfs,nBinsPerDim,'Reference Data - Rate');
plotRTbyRTchange((realData-prediction.rate_term.rt_change_pdfs)-(referenceData-prediction.reference_rate_term.rt_change_pdfs),nBinsPerDim,'(Real-Rate)-(Reference-Rate)');

realData_minusRate=realData-prediction.rate_term.rt_change_pdfs;
referenceData_minusRate=referenceData-prediction.reference_rate_term.rt_change_pdfs;

plotRTbyRTchange(realData-referenceData,nBinsPerDim,'Real Data Minus reference data');

temp=(realData-prediction.rate_term.rt_change_pdfs)-(referenceData-prediction.reference_rate_term.rt_change_pdfs);
temp=make2Dpdf(temp);
smoothSize=50;
K=ones(smoothSize);
smoothMat=conv2(temp,K,'same');
plotRTbyRTchange(temp,nBinsPerDim,'RPE term smoothed');

% Get raw data shift
% takeRTchanges=[-1 9];
% takeFirstRTs=[0 3];
% takeRTchanges=[-0.15 9]; % Good for comparison with reach but no touch
% takeFirstRTs=[0 1.5]; % Good for comparison with reach but no touch
% takeRTchanges=[-0.15 6]; 
% takeFirstRTs=[0 7.5]; 
takeRTchanges=[-0.1 9]; 
takeFirstRTs=[0 3]; 
% takeRTchanges=[0 0.5]; 
% takeFirstRTs=[0 1]; 
zonemask=getZoneMask(takeFirstRTs,takeRTchanges,bins,histo_nbins,'square'); 
% zonemask=getZoneMask(takeFirstRTs,takeRTchanges,bins,histo_nbins,'angled'); 
zonemask(zonemask==0)=nan; 
figure(); 
imagesc(bins,histo_nbins,zonemask'); 
set(gca,'YDir','normal');
tempRef=referenceData.*zonemask;
tempRef=nanmean(tempRef,1); % collapse over first RTs to get RT changes, weighted by pdf
tempData=realData.*zonemask;
tempData=nanmean(tempData,1); % collapse over first RTs to get RT changes, weighted by pdf
em=get_earthmover_1D(histo_nbins,tempRef,nansum(wrtDataset.realrtpair_seq2{1}==1),histo_nbins,tempData,nansum(dataset.realrtpair_seq2{i}==1),histo_nbins,false);

% Get rate-subtracted data shifts
takeFirstRTs=[0 9];
takeRTchanges=[-9 9];
realData_minusRate=make2Dpdf(realData_minusRate);
referenceData_minusRate=make2Dpdf(referenceData_minusRate);
tempRef=referenceData_minusRate(bins>=takeFirstRTs(1) & bins<=takeFirstRTs(2),histo_nbins>=takeRTchanges(1) & histo_nbins<=takeRTchanges(2));
tempRef=nanmean(tempRef,1); % collapse over first RTs to get RT changes, weighted by pdf
rtChangesToUse=histo_nbins(histo_nbins>=takeRTchanges(1) & histo_nbins<=takeRTchanges(2));
tempData=realData_minusRate(bins>=takeFirstRTs(1) & bins<=takeFirstRTs(2),histo_nbins>=takeRTchanges(1) & histo_nbins<=takeRTchanges(2));
tempData=nanmean(tempData,1); % collapse over first RTs to get RT changes, weighted by pdf
i=1;
em=get_earthmover_1D(rtChangesToUse,tempRef,nansum(wrtDataset.realrtpair_seq2{i}==1),rtChangesToUse,tempData,nansum(dataset.realrtpair_seq2{i}==1),histo_nbins,false);

% Get dim1, dim2 approach

end

function compareCuedRTs(firstTrial,cued_term,dataset,n_trials_away,nBinsPerDim)

rt_pdf = @(rts,bins) histcounts(rts,bins)./nansum(histcounts(rts,bins));
bin_centers = @(bins) nanmean([bins(1:end-1); bins(2:end)],1);

bins=nBinsPerDim{1};
bins=bins-(bins(2)-bins(1))/2;
bins=[bins bins(end)+(bins(2)-bins(1))/2];

rts1=dataset.event_RT_trial1InSeq{n_trials_away};
rts=dataset.event_RT_trialiInSeq{n_trials_away};

trials1_cued_rt_pdf_outs=nanmean(repmat(rt_pdf(rts1,bins)',1,size(firstTrial.rt_pdf_outs,2)).*firstTrial.rt_pdf_outs,1);
trials1_cued_rt_pdf_outs=trials1_cued_rt_pdf_outs./nansum(trials1_cued_rt_pdf_outs);
trials2_cued_rt_pdf_outs=nanmean(repmat(rt_pdf(rts,bins)',1,size(cued_term.rt_pdf_outs,2)).*cued_term.rt_pdf_outs,1);
trials2_cued_rt_pdf_outs=trials2_cued_rt_pdf_outs./nansum(trials2_cued_rt_pdf_outs);

[n,x]=cityscape_hist(trials1_cued_rt_pdf_outs,bin_centers(bins));
figure();
plot(x,n,'Color','k');
hold on;
xlabel('Time (seconds)');
ylabel('Count');
title('Comparing cued reaction times for first vs second trial in pair');
[n,x]=cityscape_hist(trials2_cued_rt_pdf_outs,bin_centers(bins));
plot(x,n,'Color','r');
leg={'trial1','trial2'};
legend(leg);

cond1_cdf=accumulateDistribution(trials1_cued_rt_pdf_outs(bin_centers(bins)<3));
bincents=bin_centers(bins);
figure();
plot(bincents(bin_centers(bins)<3),cond1_cdf./nanmax(cond1_cdf),'Color','k');
hold on;
xlabel('CDF');
ylabel('Count');

cond2_cdf=accumulateDistribution(trials2_cued_rt_pdf_outs(bin_centers(bins)<3));
plot(bincents(bin_centers(bins)<3),cond2_cdf./nanmax(cond2_cdf),'Color','r');

end

function plotRTbasics(dataset,n_trials_away,histo_nbins)

plot_rt_pdf=true;
plot_rt_cdf=true;

% Plot reaction times distribution
if plot_rt_pdf==true
    for i=n_trials_away
        histo_nbins=plotHist(dataset.allTrialsSequence_RT_trial1InSeq{i},dataset.allTrialsSequence_RT_trialiInSeq{i},histo_nbins,['Reaction Times all trials reference: trial 1 (black) vs ' num2str(dataset.nInSequence(i)-1) ' later (magenta)'],'RT (sec)');
    end
    temp=dataset.event_name;
    temp(regexp(temp,'_'))=' ';
    for i=n_trials_away
        histo_nbins=plotHist(dataset.event_RT_trial1InSeq{i},dataset.event_RT_trialiInSeq{i},histo_nbins,['Reaction Times fx of ' temp ': trial 1 (black) vs ' num2str(dataset.nInSequence(i)-1) ' later (magenta)'],'RT (sec)');
    end
end

% Plot reaction times CDF
if plot_rt_cdf==true
%     histo_nbins=backup_histo_nbins;
    for i=n_trials_away
        histo_nbins=plotCDF(dataset.allTrialsSequence_RT_trial1InSeq{i},dataset.allTrialsSequence_RT_trialiInSeq{i},histo_nbins,['CDF Reaction Times all trials reference: trial 1 (black) vs ' num2str(dataset.nInSequence(i)-1) ' later (magenta)']);
    end
    temp=dataset.event_name;
    temp(regexp(temp,'_'))=' ';
    for i=n_trials_away
        histo_nbins=plotCDF(dataset.event_RT_trial1InSeq{i},dataset.event_RT_trialiInSeq{i},histo_nbins,['CDF Reaction Times fx of ' temp ': trial 1 (black) vs ' num2str(dataset.nInSequence(i)-1) ' later (magenta)']);
    end
end

end

function timeInd=getReturnToBaseline(dataset,alltbt,cueName,n_sems,baselineWindow,n_trials_away,suppressOutput)

timeStep=mode(diff(nanmean(alltbt.times,1)));
timeBinsForReaching=0:timeStep:(size(dataset.rawReaching_event_trial1InSeq{n_trials_away},2)-1)*timeStep;

% find cue
cueInd=find(nanmean(alltbt.(cueName),1)>0,1,'first');
cueInd=cueInd+floor(0.25/timeStep); % wait at least 250 ms after cue
baselineInds=timeBinsForReaching>baselineWindow(1) & timeBinsForReaching<baselineWindow(2);

% find when reach distribution returns to baseline after cue
% within n_sems of baseline
reachMeans=dataset.rawReaching_event_trial1InSeq{n_trials_away};
reachSEMs=dataset.se_rawReaching_event_trial1InSeq{n_trials_away};
baseline=nanmean(nanmean(reachMeans(:,baselineInds),1),2);
semAtBaseline=nanmean(sqrt(nansum(reachSEMs(:,baselineInds).^2,1)));
reachesAfterCue=nanmean(reachMeans(:,cueInd:end),1);
firstAtBaseAfterCue=find(reachesAfterCue<=baseline+n_sems*semAtBaseline,1,'first');
timeInd=timeBinsForReaching(cueInd+firstAtBaseAfterCue-1);

if suppressOutput==false
    i=n_trials_away;
    plotTimeseries(dataset.rawReaching_event_trial1InSeq{i},dataset.se_rawReaching_event_trial1InSeq{i},'k',dataset.rawReaching_event_trialiInSeq{i},dataset.se_rawReaching_event_trialiInSeq{i},'m',timeBinsForReaching,baseline+n_sems*semAtBaseline,[timeInd baseline+n_sems*semAtBaseline]);
    temp=dataset.event_name;
    temp(regexp(temp,'_'))=' ';
    title(['Fx of ' temp ' first trial (black) vs trial ' num2str(dataset.nInSequence(i)-1) ' later (magenta)']);
    %legend({'me+-se','first trial','','','me+-se',['trial ' num2str(dataset.nInSequence(i)-1) ' later'],'',''});
    legend({'first trial','','',['trial ' num2str(dataset.nInSequence(i)-1) ' later'],'','','baseline','returnToBase'});
end
      
end

function plotTimeseries(data1_mean,data1_se,color1,data2_mean,data2_se,color2,timeBins,lineVal,scatterVal)

plotAsCityscape=true;

data1_mean=nanmean(data1_mean,1);
data2_mean=nanmean(data2_mean,1);
data1_se=sqrt(nansum(data1_se.^2,1));
data2_se=sqrt(nansum(data2_se.^2,1));

figure();
%fill([timeBins fliplr(timeBins)],[data1_mean+data1_se fliplr(data1_mean-data1_se)],[0.5 0.5 0.5]);
%hold on;
if plotAsCityscape==true
    [n,x]=cityscape_hist(data1_mean,timeBins);
    plot(x,n,'Color',color1); hold on;
    [n,x]=cityscape_hist(data1_mean+data1_se,timeBins);
    plot(x,n,'Color',color1);
    [n,x]=cityscape_hist(data1_mean-data1_se,timeBins);
    plot(x,n,'Color',color1);
else
    plot(timeBins,data1_mean,'Color',color1); hold on;
    plot(timeBins,data1_mean+data1_se,'Color',color1);
    plot(timeBins,data1_mean-data1_se,'Color',color1);
end

%fill([timeBins fliplr(timeBins)],[data2_mean+data2_se fliplr(data2_mean-data2_se)],[0.1 0.7 0.5]);
%hold on;
if plotAsCityscape==true
    [n,x]=cityscape_hist(data2_mean,timeBins);
    plot(x,n,'Color',color2); hold on;
    [n,x]=cityscape_hist(data2_mean+data2_se,timeBins);
    plot(x,n,'Color',color2);
    [n,x]=cityscape_hist(data2_mean-data2_se,timeBins);
    plot(x,n,'Color',color2);
else
    plot(timeBins,data2_mean,'Color',color2); hold on;
    plot(timeBins,data2_mean+data2_se,'Color',color2);
    plot(timeBins,data2_mean-data2_se,'Color',color2);
end

line([timeBins(1) timeBins(end)],[lineVal lineVal],'Color','c');
scatter(scatterVal(1),scatterVal(2));

end

function [zonemask,y_bincents]=getZoneMask(rt_range,rt_change_range,bincents,y_bincents,maskType)

if strcmp(maskType,'square')
    zonemask=zeros(length(bincents),length(y_bincents));
    zonemask(bincents>=rt_range(1) & bincents<=rt_range(2),y_bincents>=rt_change_range(1) & y_bincents<=rt_change_range(2))=1;
else
    zonemask=nan(length(bincents),length(y_bincents));
    x_on=bincents>=rt_range(1) & bincents<=rt_range(2);
    y_on=y_bincents>=rt_change_range(1) & y_bincents<=rt_change_range(2);
    y_on=y_on(y_bincents<=0);
    zonemask=nan(length(bincents),length(y_bincents));
    for i=1:length(bincents)
        startat=bincents(i);
        endInY=find(y_bincents>=startat,1,'first');
        zonemask(find(bincents>=startat,1,'first'),endInY:-1:endInY-length(y_on)+1)=fliplr(y_on);
    end
    zonemask(bincents>=0,:)=zonemask(bincents>=0,:).*repmat(x_on,size(zonemask,2),1)';
    
    y_bincents=y_bincents(y_bincents>=rt_change_range(1) & y_bincents<=rt_change_range(2));
    y_bincents=y_bincents(y_bincents<=0);
end

end

function data=make2Dpdf(data)

data=data-nanmin(data(1:end)); % non-negative
data=data./nansum(nansum(data)); % integrates to 1

end

function differ=subtractPDFfromPDF(pdf1,pdf2)

% Scale pdf2 such that pdf2 explains maximum amount of pdf1, but pdf1-pdf2
% is all non-negative
% These are assumed to be 2D pdfs

% First, ensure that nothing below zero
pdf1=pdf1-nanmin(nanmin(pdf1));
pdf2=pdf2-nanmin(nanmin(pdf2));

% Find max pdf2
[ma,maind]=nanmax(pdf2(1:end));
% Find max pdf1
[ma1,maind1]=nanmax(pdf1(1:end));
% Max of both is 1
pdf2=pdf2./ma;
pdf1=pdf1./ma1;
% Scale down pdf2 such that max of pdf2 is the value of pdf1 at this
% location
[ma,maind]=nanmax(pdf2(1:end));
% Get pdf1 value at this location
temp=pdf1(1:end);
pdf1At=temp(maind);
% Scale pdf2 such that pdf2 max is pdf1At
pdf2=(pdf2./ma)*pdf1At;
% Subtract
differ=pdf1-pdf2;
% Renormalize difference
differ=differ-nanmin(nanmin(differ));
differ=differ./nansum(nansum(differ));

end

function earth_mover_dim1=get_earthmover_1D(dataRef_x,dataRef_weights,nRefSamps,data2_x,data2_weights,nDataSamps,bins,suppressOutput) 

doBootstrap=true;

% Set all weight nans to zeros
dataRef_weights(isnan(dataRef_weights))=0;
data2_weights(isnan(data2_weights))=0;

if doBootstrap==true
    nReps=100;
    earth_mover_dim1=nan(1,nReps);
    for i=1:nReps
        curr_data1_x=randsample(dataRef_x,nRefSamps,true,dataRef_weights);
        curr_data2_x=randsample(data2_x,nDataSamps,true,data2_weights);
        [earth_mover_dim1(i)]=sub_earthmover_1D(curr_data1_x,curr_data2_x,bins);
    end
    if suppressOutput==0
        plotCDF(curr_data1_x,curr_data2_x,bins,'ref black vs real red CDF');
        em_dim1=prctile(earth_mover_dim1,[5 50 95]);
        figure();
        scatter(1,em_dim1(2));
        hold on;
        line([1 1],[em_dim1(1) em_dim1(3)]);
    end
else
    [earth_mover_dim1]=sub_earthmover_1D(data1_x,data2_x,bins);
    if suppressOutput==0
        figure();
        scatter(1,earth_mover_dim1);
    end
end

end

function [x_backup,returnThis]=plotCDF(data1,data2,bins,tit)

returnThis=[];

doKStest=false;

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
returnThis.x=x_mids;
returnThis.y=cond2_cdf./nanmax(cond2_cdf);

if doKStest==true
    [~,p]=kstest2(data1,data2);
    disp('kstest pval');
    disp(p);    
end

end

function plot_earthmover_wander(data1_x,data1_y,data2_x,data2_y,bins)

doBootstrap=true;

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
    figure();
    scatter(em_dim1(2),em_dim2(2));
    hold on;
    line([em_dim1(1) em_dim1(3)],[em_dim2(2) em_dim2(2)]);
    line([em_dim1(2) em_dim1(2)],[em_dim2(1) em_dim2(3)]);
else
    [earth_mover_dim1,earth_mover_dim2]=sub_earthmover_wander(data1_x,data1_y,data2_x,data2_y,bins);
    figure();
    scatter(earth_mover_dim1,earth_mover_dim2);
end

end

function [earth_mover_dim1]=sub_earthmover_1D(data1_x,data2_x,bins)
        
x_binsize=1;

[n,x]=histcounts(data1_x,bins);
x_mids=nanmean([x(1:end-1); x(2:end)],1);
x_cdf=accumulateDistribution(n);
x_cdf=x_cdf./nanmax(x_cdf);
x_backup=x;
[n,x]=histcounts(data2_x,bins);
x_cdf_cond2=accumulateDistribution(n);
x_cdf_cond2=x_cdf_cond2./nanmax(x_cdf_cond2);

differ=x_cdf-x_cdf_cond2; % where ref is greater than data, this is shift toward improvement, so positive
earth_mover_dim1=nansum(differ)*x_binsize;

% earth_mover_dim1=nansum(abs(x_cdf-x_cdf_cond2))*x_binsize;
% sign_of_dim1=sign(nanmean(x_cdf-x_cdf_cond2));
% earth_mover_dim1=sign_of_dim1*earth_mover_dim1;

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

function isLong=getFractionLongerTrials(alltbt,trialDuration)

disp('calculating fraction longer trials');

trialTimeStarts=nan(1,size(alltbt.timesfromarduino,1));
trialTimeEnds=nan(1,size(alltbt.timesfromarduino,1));
for i=1:size(alltbt.timesfromarduino,1)
    temp=alltbt.timesfromarduino(i,:);
    f=find(~isnan(temp),1,'first');
    if ~isempty(f)
        trialTimeStarts(i)=temp(f);
    else
        trialTimeStarts(i)=nan;
    end
    tempflip=fliplr(temp);
    f=find(~isnan(tempflip),1,'first');
    if ~isempty(f)
        trialTimeEnds(i)=tempflip(f);
    else
        trialTimeEnds(i)=nan;
    end
end
trialLengths=trialTimeEnds-trialTimeStarts;

[n,xout]=hist(trialLengths,100);
figure();
plot(xout,n);
title('Histogram of trial lengths');
xlabel('Seconds');

breakLength=19;
disp(['Using ' num2str(breakLength) ' sec as separation of short and long trials']);

isLong=trialLengths>breakLength;

end

function rateFits=getRateFits(alltbt,dataset,reachType,useCuedRateData,withinRange,cueName,trialDuration,trialStartCountAt1)

% Non-cued approach
% withinRange=[2 6];  % take reaches in this range, time is in seconds from onset of cue
% choose a start time for withinRange well after cue
nIndsToTake=200;
timeStep=mode(diff(nanmean(alltbt.times,1)));
cueInd=find(nanmean(alltbt.(cueName),1)>0,1,'first');
withinRange_inds=[cueInd+ceil(withinRange(1)/timeStep) cueInd+ceil(withinRange(2)/timeStep)];
temp=alltbt.all_reachBatch;
anyReachAligned=alignToEvent(temp,withinRange_inds,alltbt.all_reachBatch,nIndsToTake);
temp=alltbt.reachBatch_drop_reachStarts;
touchAligned=alignToEvent(temp,withinRange_inds,alltbt.all_reachBatch,nIndsToTake);
temp=alltbt.reachBatch_miss_reachStarts;
reachNoTouchAligned=alignToEvent(temp,withinRange_inds,alltbt.all_reachBatch,nIndsToTake);
temp=alltbt.reachBatch_success_reachStarts;
successTouchAligned=alignToEvent(temp,withinRange_inds,alltbt.all_reachBatch,nIndsToTake);
temp=alltbt.reachBatch_drop_reachStarts+alltbt.reachBatch_success_reachStarts;
anyTouchAligned=alignToEvent(temp,withinRange_inds,alltbt.all_reachBatch,nIndsToTake);
% Zero out first reach inds
anyReachAligned(:,1)=0;
touchAligned(:,1)=0;
reachNoTouchAligned(:,1)=0;
successTouchAligned(:,1)=0;
anyTouchAligned(:,1)=0;
figure();
plot(0:timeStep:timeStep*(length(nanmean(touchAligned,1))-1),nanmean(touchAligned,1)./timeStep,'Color','k');
hold on;
plot(0:timeStep:timeStep*(length(nanmean(reachNoTouchAligned,1))-1),nanmean(reachNoTouchAligned,1)./timeStep,'Color','r');
plot(0:timeStep:timeStep*(length(nanmean(successTouchAligned,1))-1),nanmean(successTouchAligned,1)./timeStep,'Color','g');
plot(0:timeStep:timeStep*(length(nanmean(anyTouchAligned,1))-1),nanmean(anyTouchAligned,1)./timeStep,'Color','b');
plot(0:timeStep:timeStep*(length(nanmean(anyReachAligned,1))-1),nanmean(anyReachAligned,1)./timeStep,'Color','c');
legend({'drop','miss','success','any touch','any reach'});

% there are 2 kinds of breaks between trials
% 1. 0 sec break and 2. a wheel turn without a cue, in which case there are
% 9.5 sec (or trial duration) between the end of one trial and the
% beginning of the next
% get the proportion of trials that are type 1 vs. type 2
isLong=getFractionLongerTrials(alltbt,trialDuration);
fractionIsLong=nansum(isLong)/length(isLong);
% for each rate as a function of time, take weighted average given
% proportion trial types 1 vs. 2
afterALongTimeLambda=nanmean(nanmean(anyReachAligned(:,end),1)./timeStep); % rate of Poisson (uncued) reaching after a long time (at end of trial)

rateFromNoncuedData.x=0:timeStep:timeStep*(length(nanmean(touchAligned,1))-1);
if strcmp(reachType,'drop')
    rateFromNoncuedData.y=(1-fractionIsLong)*nanmean(touchAligned,1)./timeStep + fractionIsLong.*afterALongTimeLambda.*ones(1,length(rateFromNoncuedData.x)); % reach rate in seconds
elseif strcmp(reachType,'miss')
    rateFromNoncuedData.y=(1-fractionIsLong)*nanmean(reachNoTouchAligned,1)./timeStep + fractionIsLong.*afterALongTimeLambda.*ones(1,length(rateFromNoncuedData.x));
elseif strcmp(reachType,'success')
    rateFromNoncuedData.y=(1-fractionIsLong)*nanmean(successTouchAligned,1)./timeStep + fractionIsLong.*afterALongTimeLambda.*ones(1,length(rateFromNoncuedData.x));
elseif strcmp(reachType,'anyTouch')
    rateFromNoncuedData.y=(1-fractionIsLong)*nanmean(anyTouchAligned,1)./timeStep + fractionIsLong.*afterALongTimeLambda.*ones(1,length(rateFromNoncuedData.x));
elseif strcmp(reachType,'anyReach')
    rateFromNoncuedData.y=(1-fractionIsLong)*nanmean(anyReachAligned,1)./timeStep + fractionIsLong.*afterALongTimeLambda.*ones(1,length(rateFromNoncuedData.x));
end
rateFromNoncuedData.y(3:-1:1)=rateFromNoncuedData.y(3);

figure();
title('After adjusting for fraction of long trials');
plot(rateFromNoncuedData.x,rateFromNoncuedData.y,'Color','k');
hold on;
line([0 timeStep*(length(nanmean(anyTouchAligned,1))-1)],[afterALongTimeLambda afterALongTimeLambda]);

% Fit rates using cued approach
if useCuedRateData==true
    [temp1,temp2]=fitRateTermToLateReaches_local(dataset,1); % last argument is suppressOutput
    rateFromCuedData.x=nan(1,length(temp2));
    for i=1:length(temp2)
        rateFromCuedData.x(i)=nanmean(temp2{i}); % x in terms of previous trial RT
    end
    rateFromCuedData.y=temp1; % y is rate as a function of previous trial RT
    rateFits=rateFromCuedData;
    rateFits.afterALongTimeLambda=afterALongTimeLambda;
else
    % Now adjust rateFits such that x is in terms of previous trial RT
    % and y is a rate as a function of previous trial RT
    rateFromNoncuedData.x=trialDuration-rateFromNoncuedData.x;
    rateFits=rateFromNoncuedData;
    rateFits.afterALongTimeLambda=afterALongTimeLambda;
end

if ~isempty(trialStartCountAt1)
    if trialStartCountAt1>2.5
        % means we are comparing trials separated by more than one full
        % trial in all cases
        % and I don't know when the RT was on the previous trial, so just
        % average
        rateFits.y=nanmean(rateFits.y).*ones(size(rateFits.y));
    end
end

end

function rpe_term=fitRPEtermFromData(dataset,nBinsPerDim,firstTrial,cued_term,n_trials_away,n_update_steps,sameRTforEachTrial)

% RPE term
% How the cued reaching rate changes as a result of behavior event

% Reach rate pdf as proxy for value prediction
% Value function is CDF
% Idea is that mouse reaches when mouse thinks there's a pellet
pdf = @(reas,bins) histcounts(reas,bins)./nansum(histcounts(reas,bins));
bin_centers = @(bins) nanmean([bins(1:end-1); bins(2:end)],1);

bins=nBinsPerDim{1};
bins=bins-(bins(2)-bins(1))/2;
bins=[bins bins(end)+(bins(2)-bins(1))/2];

rts=dataset.event_RT_trial1InSeq{n_trials_away};
alpha_unit=1;
% As a function of different previous reaction times
try_curr_rts=bin_centers(bins);
rt_change_pdfs=nan(length(try_curr_rts),length(bin_centers(unique([-fliplr(bins) bins]))));
rt_pdf_outs=nan(length(try_curr_rts),length(bin_centers(bins)));
for i=1:length(try_curr_rts)
    raw_pdf_out=nanmean(repmat(pdf(rts,bins)',1,size(firstTrial.raw_reaches,2)).*firstTrial.raw_reaches,1);
    rt_pdf_out=pdf(rts,bins);
    % do updates on the value function from raw reaching
    if sameRTforEachTrial==true
        for j=1:n_update_steps % if RT is same for each trial
            % put in alpha equals 1, then get real shift, compare to fit
            % alpha
            [rt_pdf_out,raw_pdf_out]=update_RPE_term_withRawReaching(rt_pdf_out,bins,try_curr_rts(i),alpha_unit,true,raw_pdf_out); 
        end
    else
        for j=1:n_update_steps % if RT is randomly sampled from RT distribution for each trial
            if j==1 
                [rt_pdf_out,raw_pdf_out]=update_RPE_term_withRawReaching(rt_pdf_out,bins,try_curr_rts(i),alpha_unit,true,raw_pdf_out);
            else
                [rt_pdf_out,raw_pdf_out]=update_RPE_term_withRawReaching(rt_pdf_out,bins,randsample(try_curr_rts,1,true,rt_pdf_out),alpha_unit,true,raw_pdf_out);
            end
        end
    end
    % assume that RPE update causes an instantaneous increase in rate of reaching,
    % followed immediately by a decrease
    % then update to RTs is same as update to the raw rate of reaching
    rt_pdf_outs(i,:)=rt_pdf_out;
    [rt_change_pdfs(i,:),rt_change_bins]=getRTchange_fromCurrAndUpdatedRTpdfs(pdf(rts,bins),rt_pdf_out,bins,true,[bins(i) bins(i+1)]);
end
rt_change_pdfs=rt_change_pdfs./nansum(nansum(rt_change_pdfs));
rpe_term.theory.rt_change_pdfs=rt_change_pdfs;
rpe_term.theory.raw_pdf_out=raw_pdf_out;
rpe_term.theory.rt_pdf_outs=rt_pdf_outs;
plotRTbyRTchange(rpe_term.theory.rt_change_pdfs,nBinsPerDim,'RPE term from theory');

% Get real change in reach rate
cuedRateUpdate=cued_term.raw_reaches-firstTrial.raw_reaches;

% Get change in RT from data
for i=1:length(try_curr_rts)
    [rt_change_pdfs(i,:),rt_change_bins]=getRTchange_fromCurrAndUpdatedRTpdfs(firstTrial.rt_pdf_outs(i,:),cued_term.rt_pdf_outs(i,:),bins,true,[bins(i) bins(i+1)]);
end
rt_change_pdfs=rt_change_pdfs./nansum(nansum(rt_change_pdfs));
rpe_term.real.rt_change_pdfs=rt_change_pdfs;
rpe_term.real.rt_change_bins=rt_change_bins;
plotRTbyRTchange(rpe_term.real.rt_change_pdfs,nBinsPerDim,'RPE term from data');

% Compare
doComparison=1;
if doComparison==1
    % no matter how many update steps, most recent update is always at bincents==0
    % rpe_term.theory.rt_change_pdfs(3,91)
    rtChangeCents=bin_centers(rt_change_bins);
    [~,mi]=nanmin(abs(rtChangeCents));
    % integrate over all previous RTs up to maxPrevRT
    maxPrevRT=3; % in seconds 
    maxInd=ceil(maxPrevRT/(bins(2)-bins(1)));
    alpha_integral_theory=nansum(rpe_term.theory.rt_change_pdfs(2:maxInd,mi),1);
    % however, there may be some slop in real reaction times
    allowedSlop=1; % in indices
    disp(['Allowing ' num2str((bins(2)-bins(1))*allowedSlop) ' seconds of RPE slop +/- for real data']);
    alpha_integral_data=nanmean(nansum(rpe_term.real.rt_change_pdfs(2:maxInd,mi-allowedSlop:mi+allowedSlop),1),2);
    fitAlpha=alpha_integral_data/alpha_integral_theory;
    disp(['Real learning rate, alpha, is ' num2str(fitAlpha)]);
    
    % Simulate
    for i=1:length(try_curr_rts)
        raw_pdf_out=nanmean(repmat(pdf(rts,bins)',1,size(firstTrial.raw_reaches,2)).*firstTrial.raw_reaches,1);
        rt_pdf_out=pdf(rts,bins);
        % do updates on the value function from raw reaching
        if sameRTforEachTrial==true
            for j=1:n_update_steps % if RT is same for each trial
                % put in alpha equals 1, then get real shift, compare to fit
                % alpha
                [rt_pdf_out,raw_pdf_out]=update_RPE_term_withRawReaching(rt_pdf_out,bins,try_curr_rts(i),fitAlpha,true,raw_pdf_out);
            end
        else
            for j=1:n_update_steps % if RT is randomly sampled from RT distribution for each trial
                if j==1
                    [rt_pdf_out,raw_pdf_out]=update_RPE_term_withRawReaching(rt_pdf_out,bins,try_curr_rts(i),fitAlpha,true,raw_pdf_out);
                else
                    [rt_pdf_out,raw_pdf_out]=update_RPE_term_withRawReaching(rt_pdf_out,bins,randsample(try_curr_rts,1,true,rt_pdf_out),fitAlpha,true,raw_pdf_out);
                end
            end
        end
        % assume that RPE update causes an instantaneous increase in rate of reaching,
        % followed immediately by a decrease
        % then update to RTs is same as update to the raw rate of reaching
        rt_pdf_outs(i,:)=rt_pdf_out;
        [rt_change_pdfs(i,:),rt_change_bins]=getRTchange_fromCurrAndUpdatedRTpdfs(pdf(rts,bins),rt_pdf_out,bins,true,[bins(i) bins(i+1)]);
    end
    rt_change_pdfs=rt_change_pdfs./nansum(nansum(rt_change_pdfs));
    rpe_term.theory_wFitAlpha.rt_change_pdfs=rt_change_pdfs;
    rpe_term.theory_wFitAlpha.raw_pdf_out=raw_pdf_out;
    rpe_term.theory_wFitAlpha.rt_pdf_outs=rt_pdf_outs;
    rpe_term.fitAlpha=fitAlpha;
    plotRTbyRTchange(rpe_term.theory_wFitAlpha.rt_change_pdfs,nBinsPerDim,'RPE term from theory with fit alpha');
end

end

function rpe_term=fitRPEtermFromTheory_local(dataset,nBinsPerDim,sameRTforEachTrial,n_update_steps,alpha,suppressOutput)

% RT pdf as proxy for value prediction
rt_pdf = @(rts,bins) histcounts(rts,bins)./nansum(histcounts(rts,bins));
bin_centers = @(bins) nanmean([bins(1:end-1); bins(2:end)],1);

i=1;
rts=dataset.event_RT_trial1InSeq{i};

% As a function of different current reaction times
bins=nBinsPerDim{1};
bins=bins-(bins(2)-bins(1))/2;
bins=[bins bins(end)+(bins(2)-bins(1))/2];
try_curr_rts=bin_centers(bins);
rt_change_pdfs=nan(length(try_curr_rts),length(bin_centers(unique([-fliplr(bins) bins]))));
rt_pdf_outs=nan(length(try_curr_rts),length(bin_centers(bins)));
for i=1:length(try_curr_rts)
    rt_pdf_out=rt_pdf(rts,bins);
    if sameRTforEachTrial==true
        for j=1:n_update_steps % if RT is same for each trial
            rt_pdf_out=update_RPE_term(rt_pdf_out,bins,try_curr_rts(i),alpha,suppressOutput);
        end
    else
        for j=1:n_update_steps % if RT is randomly sampled from RT distribution for each trial
            if j==1 
                rt_pdf_out=update_RPE_term(rt_pdf_out,bins,try_curr_rts(i),alpha,suppressOutput);
            else
                rt_pdf_out=update_RPE_term(rt_pdf_out,bins,randsample(try_curr_rts,1,true,rt_pdf_out),alpha,suppressOutput);
            end
        end
    end
    rt_pdf_outs(i,:)=rt_pdf_out;
    [rt_change_pdfs(i,:),rt_change_bins]=getRTchange_fromCurrAndUpdatedRTpdfs(rt_pdf(rts,bins),rt_pdf_out,bins,true,[bins(i) bins(i+1)]);
end
rt_change_pdfs=rt_change_pdfs./nansum(nansum(rt_change_pdfs));

rpe_term.rt_change_pdfs=rt_change_pdfs;
rpe_term.rt_change_bins=rt_change_bins;
rpe_term.rt_pdf_outs=rt_pdf_outs;
plotRTbyRTchange(rpe_term.rt_change_pdfs,nBinsPerDim,'RPE term from theory');

end

function [rt_pdf_out,raw_pdf_out]=update_RPE_term_withRawReaching(rt_pdf,bins,curr_rt,alpha,suppressOutput,curr_rawReach)

% pdf gives cdf
cdf = @(rt_pdf) accumulateDistribution(rt_pdf)./nanmax(accumulateDistribution(rt_pdf));

% cdf gives 1-cdf
rpe = @(rt_pdf) 1-cdf(rt_pdf);

% Plot
bin_centers = @(bins) nanmean([bins(1:end-1); bins(2:end)],1);
if suppressOutput==false
    pause;
    figure(); plot(bin_centers(bins),curr_rawReach,'Color','k'); title('PDF');
    figure(); plot(bin_centers(bins),cdf(curr_rawReach),'Color','k'); title('CDF');
    figure(); plot(bin_centers(bins),rpe(curr_rawReach),'Color','k'); title('RPE');
end
    
% define current RT
% curr_rt=1; % in seconds

% define learning rate
% alpha=0.03;

% 1-cdf gives update to value prediction
curr_touch = @(bins,curr_rt) find(curr_rt>=bins(1:end-1) & curr_rt<bins(2:end));
rpe_per_bin=rpe(curr_rawReach);
rpe_at_touch=zeros(size(rpe_per_bin));
rpe_at_touch(1:length(rpe_at_touch)>=curr_touch(bins,curr_rt))=rpe_per_bin(curr_touch(bins,curr_rt));
value_update = @(rt_pdf,curr_touch,rpe_at_touch,alpha) cdf(rt_pdf)+alpha.*rpe_at_touch;

% Update to value prediction gives update to pdf (both raw reaching and RT
% pdf, same update size)
rt_pdf_update = @(value_update) [0 diff(value_update)];
rt_pdf_out=rt_pdf_update(value_update(rt_pdf,curr_touch,rpe_at_touch,alpha));
raw_pdf_out=rt_pdf_update(value_update(curr_rawReach,curr_touch,rpe_at_touch,alpha));

% Plot
if suppressOutput==false
    pause;
    figure(); plot(bin_centers(bins),rpe_at_touch,'Color','k'); title('RPE at touch');
    figure(); plot(bin_centers(bins),value_update(rt_pdf,curr_touch,rpe_at_touch,alpha),'Color','k'); title('Updated value estimate');
    figure(); plot(bin_centers(bins),rt_pdf_out,'Color','k'); title('Update RT PDF');
end

end

function [rt_pdf_out,bins,whereDidUpdate]=update_RPE_term(rt_pdf,bins,curr_rt,alpha,suppressOutput)

% RT pdf gives RT cdf
rt_cdf = @(rt_pdf) accumulateDistribution(rt_pdf)./nanmax(accumulateDistribution(rt_pdf));

% RT cdf gives 1-cdf
rt_rpe = @(rt_pdf) 1-rt_cdf(rt_pdf);

% Plot
bin_centers = @(bins) nanmean([bins(1:end-1); bins(2:end)],1);
if suppressOutput==false
    pause;
    figure(); plot(bin_centers(bins),rt_pdf,'Color','k'); title('RT PDF');
    figure(); plot(bin_centers(bins),rt_cdf(rt_pdf),'Color','k'); title('RT CDF');
    figure(); plot(bin_centers(bins),rt_rpe(rt_pdf),'Color','k'); title('RT RPE');
end
    
% define current RT
% curr_rt=1; % in seconds

% define learning rate
% alpha=0.03;

% 1-cdf gives update to value prediction
curr_touch = @(bins,curr_rt) find(curr_rt>=bins(1:end-1) & curr_rt<bins(2:end));
rpe_per_bin=rt_rpe(rt_pdf);
rpe_at_touch=zeros(size(rpe_per_bin));
rpe_at_touch(1:length(rpe_at_touch)>=curr_touch(bins,curr_rt))=rpe_per_bin(curr_touch(bins,curr_rt));
value_update = @(rt_pdf,curr_touch,rpe_at_touch,alpha) rt_cdf(rt_pdf)+alpha.*rpe_at_touch;

% Update to value prediction gives update to RT pdf
rt_pdf_update = @(value_update) [0 diff(value_update)];
rt_pdf_out=rt_pdf_update(value_update(rt_pdf,curr_touch,rpe_at_touch,alpha));

% Plot
if suppressOutput==false
    pause;
    figure(); plot(bin_centers(bins),rpe_at_touch,'Color','k'); title('RPE at touch');
    figure(); plot(bin_centers(bins),value_update(rt_pdf,curr_touch,rpe_at_touch,alpha),'Color','k'); title('Updated value estimate');
    figure(); plot(bin_centers(bins),rt_pdf_out,'Color','k'); title('Update RT PDF');
end

whereDidUpdate=rpe_at_touch;

end

function cdf=accumulateDistribution(data)

cdf=nan(size(data));
for i=1:length(data)
    cdf(i)=nansum(data(1:i));
end

end

function rpe_term=fitRPEtermToEarlyReaches_local(dataset,wrtDataset,nBinsPerDim,suppressOutput)

% RPErange_prevRT=[0 2];
% RPErange_deltaRT=[-2 9];

i=1;
temp1_seq2=wrtDataset.event_RT_trial1InSeq{i};
% temp2_seq2=wrtDataset.event_RT_trialiInSeq{i};
% rtchanges_seq2=temp1_seq2(wrtDataset.realrtpair_seq2{i}==1)-temp2_seq2(wrtDataset.realrtpair_seq2{i}==1);
x2=temp1_seq2(wrtDataset.realrtpair_seq2{i}==1);
y2=wrtDataset.alldim_rtchanges_event{i};
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
if suppressOutput==0
    figure();
    imagesc(bin_c{1},bin_c{2},log(n)');
    set(gca,'YDir','normal');
end

i=1;
temp1_seq2=dataset.event_RT_trial1InSeq{i};
% temp2_seq2=dataset.event_RT_trialiInSeq{i};
% rtchanges_seq2=temp1_seq2(dataset.realrtpair_seq2{i}==1)-temp2_seq2(dataset.realrtpair_seq2{i}==1);
x2=temp1_seq2(dataset.realrtpair_seq2{i}==1);
y2=dataset.alldim_rtchanges_event{i};
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
if suppressOutput==0
    figure();
    imagesc(bin_c{1},bin_c{2},log(n_data)');
    set(gca,'YDir','normal');
end

rt_change_pdfs=n_data-n;
rt_change_pdfs=rt_change_pdfs./nansum(nansum(rt_change_pdfs));
rpe_term.rt_change_pdfs=rt_change_pdfs;
smoothSize=1;
K=ones(smoothSize);
smoothMat=conv2(rt_change_pdfs,K,'same');
plotRTbyRTchange(smoothMat,nBinsPerDim,'RPE term');

% refShiftRPE=nanmean(n(bin_c{1}>=RPErange_prevRT(1) & bin_c{1}<RPErange_prevRT(2),bin_c{2}>=RPErange_deltaRT(1) & bin_c{2}<RPErange_deltaRT(2)),1);
% disp(nanmean(refShiftRPE));
% dataShiftRPE=nanmean(n_data(bin_c{1}>=RPErange_prevRT(1) & bin_c{1}<RPErange_prevRT(2),bin_c{2}>=RPErange_deltaRT(1) & bin_c{2}<RPErange_deltaRT(2)),1);
% disp(nanmean(dataShiftRPE));

% refShiftRPE=nanmean(n(:,bin_c{2}>RPErange(1) & bin_c{2}<=RPErange(2)),2);
% dataShiftRPE=nanmean(n_data(:,bin_c{2}>RPErange(1) & bin_c{2}<=RPErange(2)),2);

% figure();
% temp=bin_c{1};
% temp=temp(bin_c{2}>=RPErange_deltaRT(1) & bin_c{2}<RPErange_deltaRT(2));
% plot(temp,refShiftRPE./nansum(refShiftRPE),'Color','k');
% hold on;
% plot(temp,dataShiftRPE./nansum(refShiftRPE),'Color','r');

end

function rtm_term=fitLocalCorrelationStructure_local(dataset,nBinsPerDim,bins,suppressOutput)

% RTM = regression to the mean

% checkRange_prevRT=[6 8];
% checkRange_deltaRT=[-20 9];

i=1;
temp1_seq2=dataset.event_RT_trial1InSeq{i};
% temp2_seq2=dataset.event_RT_trialiInSeq{i};
% rtchanges_seq2=temp1_seq2(dataset.realrtpair_seq2{i}==1)-temp2_seq2(dataset.realrtpair_seq2{i}==1);
x2=temp1_seq2(dataset.realrtpair_seq2{i}==1);
y2=dataset.alldim_rtchanges_event{i};
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
if suppressOutput==0
    figure();
    imagesc(bin_c{1},bin_c{2},log(n)');
    set(gca,'YDir','normal');
end

% RTM -- shuffled
i=1;
temp1_seq2=dataset.event_RT_trial1InSeq{i};
temp2_seq2=dataset.event_RT_trialiInSeq{i};
% shuffle the second reaction time
temp2_seq2=temp2_seq2(randperm(length(temp2_seq2)));
x2=temp1_seq2;
y2=temp1_seq2-temp2_seq2;
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
if suppressOutput==0
    figure();
    imagesc(bin_c{1},bin_c{2},log(n_data)');
    set(gca,'YDir','normal');
end

rt_pdf_outs=nan(length(bins)-1,length(bins)-1);
for i=1:length(bins)-1
    n=histcounts(y2(x2>=bins(i) & x2<=bins(i+1)),bins);
    n=n./nansum(n);
    rt_pdf_outs(i,:)=n;
end

rtm_term.rt_change_pdfs=n_data;
rtm_term.rt_change_bins=bin_c{2};
rtm_term.rt_pdf_outs=rt_pdf_outs;
plotRTbyRTchange(rtm_term.rt_change_pdfs,nBinsPerDim,'RTM term');

% realShiftRPE=nanmean(n(bin_c{1}>=checkRange_prevRT(1) & bin_c{1}<checkRange_prevRT(2),bin_c{2}>=checkRange_deltaRT(1) & bin_c{2}<checkRange_deltaRT(2)),1);
% disp(nanmean(realShiftRPE));
% shuffleShiftRPE=nanmean(n_data(bin_c{1}>=checkRange_prevRT(1) & bin_c{1}<checkRange_prevRT(2),bin_c{2}>=checkRange_deltaRT(1) & bin_c{2}<checkRange_deltaRT(2)),1);
% disp(nanmean(shuffleShiftRPE));
% 
% figure();
% temp=bin_c{1};
% temp=temp(bin_c{2}>=checkRange_deltaRT(1) & bin_c{2}<checkRange_deltaRT(2));
% plot(temp,realShiftRPE./nansum(realShiftRPE),'Color','k');
% hold on;
% plot(temp,shuffleShiftRPE./nansum(shuffleShiftRPE),'Color','r');
% legend({'real', 'shuffled'});

end

function n_data=plotRTbyRTchangePlot(x2,y2,nBinsPerDim,tit)

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
xlabel('Reaction time prev trial');
ylabel('Change in reaction times');
title(tit);

end

function plotRTbyRTchange(n_data,bin_c,tit)

smoothSize=1;
K=ones(smoothSize);
smoothMat=conv2(n_data,K,'same');
figure();
if any(smoothMat(1:end)<0)
    imagesc(bin_c{1},bin_c{2},smoothMat');
%     smoothMat=smoothMat-nanmin(smoothMat(1:end));
else
    imagesc(bin_c{1},bin_c{2},log(smoothMat)');
end
set(gca,'YDir','normal');
xlabel('Reaction time prev trial');
ylabel('Change in reaction times');
title(tit);

end

function cued_term=getRawCuedReachingFromRTpdf(nBinsPerDim,cued_term,dataset,n_trials_away,suppressOutput,plotRTrange)

pdf = @(rts,bins) histcounts(rts,bins)./nansum(histcounts(rts,bins));
bin_centers = @(bins) nanmean([bins(1:end-1); bins(2:end)],1);

bins=nBinsPerDim{1};
bins=bins-(bins(2)-bins(1))/2;
bins=[bins bins(end)+(bins(2)-bins(1))/2];

% note that cued reach rate at time t_i is cuedLambdas at time t_i x
% cued_shape_params at time t_i
cued_term.raw_reaches=cued_term.cuedLambdas.*cued_term.cued_shape_params;

% plot 
if suppressOutput==false
    bincents=bin_centers(bins);
    figure();
    imagesc(bincents,bincents,cued_term.raw_reaches);
    title('Cued raw reaching as a function of prev trial RT');
    xlabel('Time (sec)');
    ylabel('Prev trial RT (sec)');
    
    % Plot expected raw cued reaching rate across all previous trial RTs
    % by taking average weighted by frequency of each trial type
    rts=dataset.event_RT_trial1InSeq{n_trials_away};
    rt_pdf=pdf(rts,bins);
    figure();
    if ~isempty(plotRTrange)
        % take only prev rt bins in this range
        temp=repmat(rt_pdf',1,size(cued_term.raw_reaches,2)).*cued_term.raw_reaches;
        takeBins=bincents>=plotRTrange(1) & bincents<=plotRTrange(2);
        plot(bincents,nanmean(temp(takeBins,:),1));
    else    
        plot(bincents,nanmean(repmat(rt_pdf',1,size(cued_term.raw_reaches,2)).*cued_term.raw_reaches,1));
    end
    xlabel('Time (sec)');
    ylabel('Cued reach rate');
end

end

function [rt_pdf_out,best_a,rt_pdfs_combo_out,cuedLambdas]=solve_for_cued_term(bins,rate_term,empirical_rts2,empirical_rts1,dataset,realDataset,alltbt,n_trials_away,cueName)

rt_pdf = @(rts,bins) histcounts(rts,bins)./nansum(histcounts(rts,bins));
bin_centers = @(bins) nanmean([bins(1:end-1); bins(2:end)],1);

offsettime=0.0075; % this should be a very small number
bins=bins+offsettime; % avoid infinities at gampdf when t=0
bincents=bin_centers(bins);
rt_pdf_out=nan(length(bincents),length(bincents));
rt_pdfs_combo_out=nan(length(bincents),length(bincents));
best_a=nan(length(bincents),length(bincents)); 
cuedLambdas=nan(length(bincents),length(bincents)); 
% estimate the proportion of cued to non-cued reaching using raw reaches
timeStep=mode(diff(nanmean(alltbt.times,1)));
cueInd=find(nanmean(alltbt.(cueName),1)>0,1,'first');
rawReaches=realDataset.rawReaching_event_trialiInSeq{n_trials_away};
% then get rts from trial 1 from real dataset
realData_rts1=realDataset.event_RT_trial1InSeq{n_trials_away};
if length(realData_rts1)~=size(rawReaches,1)
    if size(rawReaches,1)>1
        rawReaches=nanmean(rawReaches,1);
    end
    warning('problem estimating cued vs uncued rates');
end
totalRate.x=0+offsettime:timeStep:timeStep*(length(cueInd:size(rawReaches,2))-1)+offsettime;
for i=1:length(bincents)
    % get current rate term RT pdf
    ratePDF=rate_term.rt_pdf_outs(i,:);
    % get lambdas from rate term
    lambdas=rate_term.lambdas(i,:); % note that lambda is 1/b when switch to Matlab's gamma pdf parameterization w b
    % get empirical RT pdf
    empiricalPDF=rt_pdf(empirical_rts2(empirical_rts1>=bins(i) & empirical_rts1<=bins(i+1)),bins);
    % get estimate of cued and uncued rates
    % be sure that have enough trials with data to get a good estimate of
    % the total reach rate 
    if nansum(realData_rts1>=bins(i) & realData_rts1<=bins(i+1))<200
        % expand bin size
        bigger_i=i-3:i+3;
        bigger_i=bigger_i(bigger_i>0 & bigger_i<=length(bins)-1);
        curr_rawReaches=nanmean(rawReaches(realData_rts1>=bins(bigger_i(1)) & realData_rts1<=bins(bigger_i(end)+1),:),1);
    else
        curr_rawReaches=nanmean(rawReaches(realData_rts1>=bins(i) & realData_rts1<=bins(i+1),:),1);
    end
    curr_rawReaches=curr_rawReaches(cueInd:end);
    totalRate.y=curr_rawReaches./timeStep;
    % consider moment-by-moment
    % then Gamma_empirical(1+a,1/(lambda+beta)) % know lambda and can estimate beta, need to fit a
    % thus, for each time point, solve for parameter a such that
    % Gamma_empirical(1+a,1/(lambda+beta)) = empiricalPDF(i,time point)
    try_a=0.01:0.02:10; % a is in units of lambda
    disp('Fitting cued component of reaching');
    for j=1:length(bincents)
        % each time point
        currtime=bincents(j);
        % lambda at this time delay is lambdas(j)
        % get total reach rate at this time delay
        [~,mi]=nanmin(abs(totalRate.x-currtime));
        currTotRate=totalRate.y(mi);
        estimateCuedRate=currTotRate-lambdas(j); % call this beta, the estimated cued rate
        % note that estimateCuedRate doesn't make sense < 0
        % moreover, as estimateCuedRate gets very small, scale parameter of
        % gamma gets very big, which produces a flat curve near zero
        if estimateCuedRate<0
            estimateCuedRate=0;
        end
        cuedLambdas(i,j)=estimateCuedRate; % this is an estimate of the rate for the cued process
        y=empiricalPDF(j); % Gamma_empirical at this time delay
        errs=nan(1,length(try_a));
        for k=1:length(try_a)
            curr_a=try_a(k);
            b=1/(lambdas(j)+estimateCuedRate); % note that this is just 1/currTotRate
            curr_y_guess=(1/((b.^(1+curr_a))*gamma(1+curr_a)))*(currtime^((1+curr_a)-1))*(exp(-currtime/b)); % gamma pdf
            errs(k)=(y-curr_y_guess).^2;
        end
        [~,mi]=nanmin(abs(errs));
        best_a(i,j)=try_a(mi);
        % therefore gampdf for combined cued and uncued is
        temp=gampdf(currtime,1+best_a(i,j),b);
        if isnan(y)
            rt_pdfs_combo_out(i,j)=nan;
        else
            rt_pdfs_combo_out(i,j)=temp;
        end
        % therefore gampdf for cued component specifically is 
        if estimateCuedRate==0
            temp=0;
        else
            temp=gampdf(currtime,best_a(i,j),1/estimateCuedRate);
        end
        if isnan(y)
            rt_pdf_out(i,j)=nan; % will fill in later with nearest neighbor data
        else
            rt_pdf_out(i,j)=temp;
        end
    end
end

% fill in missing data in rt_pdf_out
rt_pdf_out=fillIn_wNearestRows(rt_pdf_out);
% for i=1:10
%     rt_pdf_out=fillIn2Dnans_wNearestNeighbor(rt_pdf_out);
%     if all(~isnan(rt_pdf_out(1:end)))
%         break
%     end
% end

for i=1:size(rt_pdf_out,1)
    if all(rt_pdf_out(i,:)==0)
        rt_pdf_out(i,1)=0.000001;
    end
end

end

function filledin_data=fillIn_wNearestRows(data)

data(isnan(data))=0;

filledin_data=data;
for i=1:size(data,1)
    if all(data(i,:)<0.01)
        % fill in with nearest populated row
        % direction 1
        dir1_val=nan;
        dir1_distance=nan;
        for j=i+1:size(data,1)
            if ~all(data(j,:)<0.01)
                dir1_val=data(j,:);
                dir1_distance=j-i;
                break
            end
        end
        % direction 2
        dir2_val=nan;
        dir2_distance=nan;
        for j=i-1:-1:1
            if ~all(data(j,:)<0.01)
                dir2_val=data(j,:);
                dir2_distance=i-j;
                break
            end
        end
        % take value from min distance
        [themin,mi]=nanmin([dir1_distance dir2_distance]);
        if isnan(themin)
            % couldn't find neighboring row with data
            continue
        else
            if mi==1
                filledin_data(i,:)=dir1_val;
            elseif mi==2
                filledin_data(i,:)=dir2_val;
            end
        end
    end     
end

end

function filledin_data=fillIn2Dnans_wNearestNeighbor(data)

% assume a continuous 2D function, where some data is missing, indicated by
% presence of nans

% fill in data at nans with the value of the nearest data-populated
% neighbor, in either dimension

filledin_data=data;
for i=1:size(data,1)
    for j=1:size(data,2)
        if isnan(data(i,j))
            % look in both dimensions for nearest filled in data
            % choose data that is closest 
            % Direction 1
            dir1_val=nan;
            dir1_distance=nan;
            for k=j+1:size(data,2)
                if ~isnan(data(i,k))
                    dir1_val=data(i,k);
                    dir1_distance=k-j;
                    break
                end
            end
            % Direction 2
            dir2_val=nan;
            dir2_distance=nan;
            for k=j-1:-1:1
                if ~isnan(data(i,k))
                    dir2_val=data(i,k);
                    dir2_distance=j-k;
                    break
                end
            end
            % Direction 3
            dir3_val=nan;
            dir3_distance=nan;
            for k=i+1:size(data,1)
                if ~isnan(data(k,j))
                    dir3_val=data(k,j);
                    dir3_distance=k-i;
                    break
                end
            end
            % Direction 4
            dir4_val=nan;
            dir4_distance=nan;
            for k=i-1:-1:1
                if ~isnan(data(k,j))
                    dir4_val=data(k,j);
                    dir4_distance=i-k;
                    break
                end
            end
            % take value from min distance
            [themin,mi]=nanmin([dir1_distance dir2_distance dir3_distance dir4_distance]);
            if isnan(themin)
                % couldn't find neighbor with data
                continue
            else
                temp=[dir1_val dir2_val dir3_val dir4_val];
                filledin_data(i,j)=temp(mi);
            end
        end
    end
end

end

function checkTermsTogether(dataset,rate_term,cued_term,nBinsPerDim,n_trials_away)

bin_centers = @(bins) nanmean([bins(1:end-1); bins(2:end)],1);

bins=nBinsPerDim{1};
bins=bins-(bins(2)-bins(1))/2;
bins=[bins bins(end)+(bins(2)-bins(1))/2];
try_curr_rts=bin_centers(bins);

% Compare predictions to real data
j=0:0.5:8.5;
for i=1:length(j)
   if j(i)<3
       rtBins{i}=[j(i) j(i)+1]; 
   elseif j(i)<6
       rtBins{i}=[j(i) j(i)+1.5]; 
   else
       rtBins{i}=[j(i) j(i)+3]; 
   end
   nowRTbins=rtBins{i};
   temp1_seq2=dataset.event_RT_trial1InSeq{n_trials_away};
   b=temp1_seq2(dataset.realrtpair_seq2{n_trials_away}==1);
   dimall_event=dataset.alldim_rtchanges_event{n_trials_away};
   temp_RTdiffsevent=dimall_event(b>=nowRTbins(1,1) & b<nowRTbins(1,2));
   dataToFit=temp_RTdiffsevent;
   firstRTs=b(b>=nowRTbins(1,1) & b<nowRTbins(1,2));
   histo_nbins=[-4*12.4245:0.2510:4*12.4245];
   [~,mif]=nanmin(abs(try_curr_rts-nanmean(firstRTs)));
   % Sample from uncued process, let these RTs constitute a fraction of
   % total RTs
   % Sample from cued process, let these RTs constitute a fraction of total
   % RTs 
   % Get the probability of RT coming from each process (cued or uncued)
   % after integrating the probabilities as a function of time over all times (for
   % this prev RT)
   % This only works with rt_pdf_outs if have maintained scale of
   % rt_pdf_outs in terms of raw probabilities
   uncued_integral=nansum(rate_term.rt_pdf_outs(i,:));
   cued_integral=nansum(cued_term.rt_pdf_outs(i,:));
   fractionUncued=uncued_integral/(uncued_integral+cued_integral);
   pred_uncued=randsample(bin_centers(bins)+(bins(2)-bins(1)),floor(length(firstRTs)*fractionUncued),true,rate_term.rt_pdf_outs(mif,:));
   pred_cued=randsample(bin_centers(bins)+(bins(2)-bins(1)),ceil(length(firstRTs)*(1-fractionUncued)),true,cued_term.rt_pdf_outs(mif,:));
   pred=[pred_uncued pred_cued]; 
   pred=pred(1:length(firstRTs));
   prediction=firstRTs-pred; 
   [~,~,fHist]=plotHist(dataToFit,prediction,histo_nbins,['Checking fits for rate and cued term together for firstRTs ' num2str(nanmean(firstRTs))],'Change in RT (sec)');
   pause;
   close(fHist);
end

end

function rts=getRTsFromNeighboringTrials(alltbt,cueName,nTrialsAheadForRT1,isSequence)

if size(alltbt.(cueName),1)~=length(isSequence)
    error('alltbt size does not match sequence length');
end
indsOfCurrRTs=find(isSequence==1);
indsOfCurrRTs=indsOfCurrRTs+nTrialsAheadForRT1; % if -1, will take previous trial before sequence starts

cueInd=find(nanmean(alltbt.(cueName),1)>0,1,'first');
timeStep=mode(diff(nanmean(alltbt.times,1)));

reaches=alltbt.all_reachBatch;
rts=nan(1,length(indsOfCurrRTs));
for i=1:length(indsOfCurrRTs)
    currInd=indsOfCurrRTs(i);
    if currInd<1 || currInd>size(reaches,1)
        continue
    end
    % get RT
    temp=reaches(currInd,:);
    aftercue_temp=temp(cueInd:end);
    f=find(aftercue_temp>0.05,1,'first');
    if ~isempty(f)
        rts(i)=(f-1)*timeStep;
    end
end

end

function cued_term=firstTrial_solveForCuedReaching(dataset,rateFits,nBinsPerDim,n_trials_away,alltbt,cueName,baselineWindow,nTrialsAheadForRT1,realDataset,getRTsFromTbt)

% uncued reach rate in this case is just the rate from raw reaches from
% before the cue
timeStep=mode(diff(nanmean(alltbt.times,1)));
timeBinsForReaching=0:timeStep:(size(dataset.rawReaching_event_trial1InSeq{n_trials_away},2)-1)*timeStep;
baselineInds=timeBinsForReaching>baselineWindow(1) & timeBinsForReaching<baselineWindow(2);
rawReaches=nanmean(dataset.rawReaching_event_trial1InSeq{n_trials_away},1)./timeStep; % in terms of rates

% lambdas=nanmean(rawReaches(baselineInds)); % uncued reach rate
% use uncued rate from prior calculations
lambdas=rateFits.afterALongTimeLambda;

rt_pdf = @(rts,bins) histcounts(rts,bins)./nansum(histcounts(rts,bins));
bin_centers = @(bins) nanmean([bins(1:end-1); bins(2:end)],1);

rts=dataset.event_RT_trial1InSeq{n_trials_away};

bins=nBinsPerDim{1};
bins=bins-(bins(2)-bins(1))/2;
bins=[bins bins(end)+(bins(2)-bins(1))/2];

offsettime=0.0075; % this should be a very small number
bins=bins+offsettime; % avoid infinities at gampdf when t=0
bincents=bin_centers(bins);
rt_pdf_out=nan(length(bincents),length(bincents));
rt_pdfs_combo_out=nan(length(bincents),length(bincents));
best_a=nan(length(bincents),length(bincents)); 
cuedLambdas=nan(length(bincents),length(bincents)); 
% estimate the proportion of cued to non-cued reaching using raw reaches
cueInd=find(nanmean(alltbt.(cueName),1)>0,1,'first');
totalRate.x=0+offsettime:timeStep:timeStep*(length(cueInd:size(rawReaches,2))-1)+offsettime;
totalRate.y=rawReaches(cueInd:end);
if isempty(nTrialsAheadForRT1)
    for i=1:1 % note that there is no conditionality to the RT pdf in this case, thus all rows are the same
        % get empirical RT pdf
        empiricalPDF=rt_pdf(rts,bins);
        % consider moment-by-moment
        % then Gamma_empirical(1+a,1/(lambda+beta)) % know lambda and can estimate beta, need to fit a
        % thus, for each time point, solve for parameter a such that
        % Gamma_empirical(1+a,1/(lambda+beta)) = empiricalPDF(i,time point)
        try_a=0.01:0.02:10; % a is in units of lambda
        disp('Fitting cued component of reaching');
        for j=1:length(bincents)
            % each time point
            currtime=bincents(j);
            % lambda at this time delay is lambdas
            % get total reach rate at this time delay
            [~,mi]=nanmin(abs(totalRate.x-currtime));
            currTotRate=totalRate.y(mi);
            estimateCuedRate=currTotRate-lambdas; % call this beta, the estimated cued rate
            % note that estimateCuedRate doesn't make sense < 0
            % moreover, as estimateCuedRate gets very small, scale parameter of
            % gamma gets very big, which produces a flat curve near zero
            if estimateCuedRate<0
                estimateCuedRate=0;
            end
            cuedLambdas(i,j)=estimateCuedRate; % this is an estimate of the rate for the cued process
            y=empiricalPDF(j); % Gamma_empirical at this time delay
            errs=nan(1,length(try_a));
            for k=1:length(try_a)
                curr_a=try_a(k);
                b=1/(lambdas+estimateCuedRate); % note that this is just 1/currTotRate
                curr_y_guess=(1/((b.^(1+curr_a))*gamma(1+curr_a)))*(currtime^((1+curr_a)-1))*(exp(-currtime/b)); % gamma pdf
                errs(k)=(y-curr_y_guess).^2;
            end
            [~,mi]=nanmin(abs(errs));
            best_a(i,j)=try_a(mi);
            % therefore gampdf for combined cued and uncued is
            temp=gampdf(currtime,1+best_a(i,j),b);
            if isnan(y)
                rt_pdfs_combo_out(i,j)=nan;
            else
                rt_pdfs_combo_out(i,j)=temp;
            end
            % therefore gampdf for cued component specifically is
            if estimateCuedRate==0
                temp=0;
            else
                temp=gampdf(currtime,best_a(i,j),1/estimateCuedRate);
            end
            if isnan(y)
                rt_pdf_out(i,j)=0;
            else
                rt_pdf_out(i,j)=temp;
            end
        end
    end
    % copy first row to all other rows
    rt_pdf_out(2:end,:)=repmat(rt_pdf_out(1,:),size(rt_pdf_out,1)-1,1);
    rt_pdfs_combo_out(2:end,:)=repmat(rt_pdfs_combo_out(1,:),size(rt_pdf_out,1)-1,1);
    best_a(2:end,:)=repmat(best_a(1,:),size(rt_pdf_out,1)-1,1);
    cuedLambdas(2:end,:)=repmat(cuedLambdas(1,:),size(rt_pdf_out,1)-1,1);
else
    if getRTsFromTbt==true
        rts1=getRTsFromNeighboringTrials(alltbt,cueName,nTrialsAheadForRT1,realDataset.event_isSeq{n_trials_away});
        rawReaches=realDataset.rawReaching_event_trial1InSeq{n_trials_away};
        realData_rts1=rts1;
    else
        rts1=dataset.event_RT_trialiInSeq{n_trials_away};
        rawReaches=realDataset.rawReaching_event_trial1InSeq{n_trials_away};
        % shuffle locally to remove close neighbors' correlation structure while
        % preserving correlation structures over longer times
        % shuffle rts1
        inds=shuffleLocally(1:length(rts1),20);
        rts1=rts1(inds);
        % shuffle realData_rts1
        realData_rts1=realDataset.event_RT_trialiInSeq{n_trials_away};
        inds=shuffleLocally(1:length(realData_rts1),20);
        realData_rts1=realData_rts1(inds);
    end
    for i=1:length(bincents)
        % get empirical RT pdf
        empiricalPDF=rt_pdf(rts(rts1>=bins(i) & rts1<=bins(i+1)),bins);
        % get estimate of cued and uncued rates
        % be sure that have enough trials with data to get a good estimate of
        % the total reach rate
        if nansum(realData_rts1>=bins(i) & realData_rts1<=bins(i+1))<200
            % expand bin size
            bigger_i=i-3:i+3;
            bigger_i=bigger_i(bigger_i>0 & bigger_i<=length(bins)-1);
            curr_rawReaches=nanmean(rawReaches(realData_rts1>=bins(bigger_i(1)) & realData_rts1<=bins(bigger_i(end)+1),:),1);
        else
            curr_rawReaches=nanmean(rawReaches(realData_rts1>=bins(i) & realData_rts1<=bins(i+1),:),1);
        end
        curr_rawReaches=curr_rawReaches(cueInd:end);
        totalRate.y=curr_rawReaches./timeStep;
        % consider moment-by-moment
        % then Gamma_empirical(1+a,1/(lambda+beta)) % know lambda and can estimate beta, need to fit a
        % thus, for each time point, solve for parameter a such that
        % Gamma_empirical(1+a,1/(lambda+beta)) = empiricalPDF(i,time point)
        try_a=0.01:0.02:10; % a is in units of lambda
        disp('Fitting cued component of reaching');
        for j=1:length(bincents)
            % each time point
            currtime=bincents(j);
            % lambda at this time delay is lambdas(j)
            % get total reach rate at this time delay
            [~,mi]=nanmin(abs(totalRate.x-currtime));
            currTotRate=totalRate.y(mi);
            estimateCuedRate=currTotRate-lambdas; % call this beta, the estimated cued rate
            % note that estimateCuedRate doesn't make sense < 0
            % moreover, as estimateCuedRate gets very small, scale parameter of
            % gamma gets very big, which produces a flat curve near zero
            if estimateCuedRate<0
                estimateCuedRate=0;
            end
            cuedLambdas(i,j)=estimateCuedRate; % this is an estimate of the rate for the cued process
            y=empiricalPDF(j); % Gamma_empirical at this time delay
            errs=nan(1,length(try_a));
            for k=1:length(try_a)
                curr_a=try_a(k);
                b=1/(lambdas+estimateCuedRate); % note that this is just 1/currTotRate
                curr_y_guess=(1/((b.^(1+curr_a))*gamma(1+curr_a)))*(currtime^((1+curr_a)-1))*(exp(-currtime/b)); % gamma pdf
                errs(k)=(y-curr_y_guess).^2;
            end
            [~,mi]=nanmin(abs(errs));
            best_a(i,j)=try_a(mi);
            % therefore gampdf for combined cued and uncued is
            temp=gampdf(currtime,1+best_a(i,j),b);
            if isnan(y)
                rt_pdfs_combo_out(i,j)=nan;
            else
                rt_pdfs_combo_out(i,j)=temp;
            end
            % therefore gampdf for cued component specifically is
            if estimateCuedRate==0
                temp=0;
            else
                temp=gampdf(currtime,best_a(i,j),1/estimateCuedRate);
            end
            if isnan(y)
                rt_pdf_out(i,j)=0;
            else
                rt_pdf_out(i,j)=temp;
            end
        end
    end
end

for i=1:size(rt_pdf_out,1)
    if all(rt_pdf_out(i,:)==0)
        rt_pdf_out(i,1)=0.000001;
    end
end

cued_term.rt_pdf_outs=rt_pdf_out;
cued_term.cued_shape_params=best_a;
cued_term.cuedLambdas=cuedLambdas;

end

function inds=shuffleLocally(inds,maxDistanceAllowed)

for i=1:maxDistanceAllowed:length(inds)
    if i+maxDistanceAllowed-1>length(inds)
        currinds=i:length(inds);
        f=randperm(length(currinds));
        inds(i:length(inds))=currinds(f);
    else
        currinds=i:i+maxDistanceAllowed-1;
        f=randperm(length(currinds));
        inds(i:i+maxDistanceAllowed-1)=currinds(f);
    end
end

end

function cued_term=solveForCuedReaching(dataset,realDataset,rate_term,nBinsPerDim,n_trials_away,alltbt,cueName)

% for each value of the previous RT (i.e., RT on trial 1), get the RTs in 
% real data set that are not explained by uncued rate

rt_pdf = @(rts,bins) histcounts(rts,bins)./nansum(histcounts(rts,bins));
bin_centers = @(bins) nanmean([bins(1:end-1); bins(2:end)],1);

rts=dataset.event_RT_trial1InSeq{n_trials_away};
rts2=dataset.event_RT_trialiInSeq{n_trials_away};

bins=nBinsPerDim{1};
bins=bins-(bins(2)-bins(1))/2;
bins=[bins bins(end)+(bins(2)-bins(1))/2];

% As a function of different current reaction times
try_curr_rts=bin_centers(bins);
rt_change_pdfs=nan(length(try_curr_rts),length(bin_centers(unique([-fliplr(bins) bins]))));
[rt_pdf_outs,cued_shape_params,~,cuedLambdas]=solve_for_cued_term(bins,rate_term,rts2,rts,dataset,realDataset,alltbt,n_trials_away,cueName);
for i=1:length(try_curr_rts)
    rt_pdf_out=rt_pdf_outs(i,:);
    [rt_change_pdfs(i,:),rt_change_bins]=getRTchange_fromCurrAndUpdatedRTpdfs(rt_pdf(rts,bins),rt_pdf_out,bins,true,[bins(i) bins(i+1)]);
end
rt_change_pdfs=rt_change_pdfs./nansum(nansum(rt_change_pdfs));

cued_term.rt_change_pdfs=rt_change_pdfs;
cued_term.rt_change_bins=rt_change_bins;
cued_term.rt_pdf_outs=rt_pdf_outs;
cued_term.cued_shape_params=cued_shape_params;
cued_term.cuedLambdas=cuedLambdas; % note that in Matlab's parameterization of gamma distribution, second parameter b is 1/lambda

plotRTbyRTchange(rt_change_pdfs,nBinsPerDim,'Cued term');

end

function rate_term=makeRateTerm(nBinsPerDim,rateFits,dataset,n_trials_away)

rt_pdf = @(rts,bins) histcounts(rts,bins)./nansum(histcounts(rts,bins));
bin_centers = @(bins) nanmean([bins(1:end-1); bins(2:end)],1);

i=n_trials_away;
rts=dataset.event_RT_trial1InSeq{i};

bins=nBinsPerDim{1};
bins=bins-(bins(2)-bins(1))/2;
bins=[bins bins(end)+(bins(2)-bins(1))/2];

% As a function of different current reaction times
try_curr_rts=bin_centers(bins);
rt_change_pdfs=nan(length(try_curr_rts),length(bin_centers(unique([-fliplr(bins) bins]))));
[rt_pdf_outs,lambdas]=update_rate_term(bins,rateFits,dataset);
for i=1:length(try_curr_rts)
    rt_pdf_out=rt_pdf_outs(i,:);
    [rt_change_pdfs(i,:),rt_change_bins]=getRTchange_fromCurrAndUpdatedRTpdfs(rt_pdf(rts,bins),rt_pdf_out,bins,true,[bins(i) bins(i+1)]);
end
rt_change_pdfs=rt_change_pdfs./nansum(nansum(rt_change_pdfs)); % renormalize here

rate_term.rt_change_pdfs=rt_change_pdfs;
rate_term.rt_change_bins=rt_change_bins;
rate_term.rt_pdf_outs=rt_pdf_outs;
rate_term.lambdas=lambdas;
rate_term.afterALongTimeLambda=rateFits.afterALongTimeLambda;

nBinsPerDim={try_curr_rts; bin_centers(rt_change_bins)};
plotRTbyRTchange(rt_change_pdfs,nBinsPerDim,'Rate term');

end

function [rt_change_pdf,rt_change_bins]=getRTchange_fromCurrAndUpdatedRTpdfs(curr_rt_pdf,updated_rt_pdf,bins,suppressOutput,curr_rt)

% Get the distribution of reaction time changes expected given the input
% reaction time distribution (curr_rt_pdf) and the reaction time
% distribution on the next trial (updated_rt_pdf)

% This is the difference of these two random variables
sample_n=30000;
bin_centers = @(bins) nanmean([bins(1:end-1); bins(2:end)],1);
rts1=randsample(bin_centers(bins),sample_n,true,curr_rt_pdf);
rts2=randsample(bin_centers(bins),sample_n,true,updated_rt_pdf);
rt_diff=rts1-rts2;
[n,x]=histcounts(rts1,bins);
if suppressOutput==false
    [n,x]=cityscape_hist(n,x);
    figure();
    plot(x,n./nansum(n),'Color','b');
    hold on;
    [n,x]=cityscape_hist(curr_rt_pdf,bins);
    plot(x,n./nansum(n),'Color','r');
    legend({'Resample RT PDF','RT PDF'});
    xlabel('Reaction times (sec)');
    ylabel('Count');
    title('Reaction time PDF');
end
if ~isempty(curr_rt)
    % will take only rt pairs for which first rt is curr_rt
    rt_diff=rt_diff(rts1>curr_rt(1) & rts1<=curr_rt(2));
else 
    % will take all rt pairs for differences
end
[n,x]=histcounts(rt_diff,unique([-fliplr(bins) bins]));
rt_change_pdf=n;
rt_change_bins=x;

if suppressOutput==false
    [n,x]=cityscape_hist(n,x);
    figure();
    plot(x,n./nansum(n),'Color','k');
    xlabel('Reaction times (sec)');
    ylabel('Count');
    title('Change in reaction time');
end

end

function [rt_pdf_out,lambdas]=update_rate_term(bins,fits,dataset)

bin_centers = @(bins) nanmean([bins(1:end-1); bins(2:end)],1);

bincents=bin_centers(bins);
rt_pdf_out=nan(length(bincents),length(bincents));
lambdas=nan(length(bincents),length(bincents));
for i=1:length(bincents)
    k=1;
    n=nan(1,length(bincents));
    l=nan(1,length(bincents));
    % first read off starting lambda (reach rate), given that reaction time
    % of first trial was bincents(i)
    [~,mi]=min(abs(fits.x-bincents(i)));
    currLambda=fits.y(mi);
    currtime=0;
    n(k)=currLambda*exp(-currLambda*currtime);
    l(k)=currLambda;
    k=k+1;
    % then model rest of curve as lambda (reach rate) returns to baseline
    currtime=currtime+(bins(2)-bins(1));
    for k=k:length(n)
        mi=mi+1;
        if mi>length(fits.y)
            currLambda=fits.y(end);
        else
            currLambda=fits.y(mi);
        end
        currtime=currtime+(bins(2)-bins(1));
        n(k)=currLambda*exp(-currLambda*currtime);
        l(k)=currLambda;
    end
    rt_pdf_out(i,:)=n;
    lambdas(i,:)=l;
end

end

function checkRates(rateFits,dataset,nBinsPerDim,n_trials_away)

rt_pdf = @(rts,bins) histcounts(rts,bins)./nansum(histcounts(rts,bins));
bin_centers = @(bins) nanmean([bins(1:end-1); bins(2:end)],1);

i=n_trials_away;
rts=dataset.event_RT_trial1InSeq{i};

bins=nBinsPerDim{1};
bins=bins-(bins(2)-bins(1))/2;
bins=[bins bins(end)+(bins(2)-bins(1))/2];

% Get RT pdfs on next trial as a function of different current reaction times
% Predicted from rateFits
try_curr_rts=bin_centers(bins);
rt_change_pdfs=nan(length(try_curr_rts),length(bin_centers(unique([-fliplr(bins) bins]))));
[rt_pdf_outs,lambdas]=update_rate_term(bins,rateFits,dataset);
for i=1:length(try_curr_rts)
    rt_pdf_out=rt_pdf_outs(i,:);
    [rt_change_pdfs(i,:),rt_change_bins]=getRTchange_fromCurrAndUpdatedRTpdfs(rt_pdf(rts,bins),rt_pdf_out,bins,true,[bins(i) bins(i+1)]);
end
rt_change_pdfs=rt_change_pdfs./nansum(nansum(rt_change_pdfs)); 

rate_term.rt_change_pdfs=rt_change_pdfs;
rate_term.rt_change_bins=rt_change_bins;
rate_term.rt_pdf_outs=rt_pdf_outs;
rate_term.lambdas=lambdas;

% Compare predictions to real data
j=0:0.5:8.5;
for i=1:length(j)
   if j(i)<3
       rtBins{i}=[j(i) j(i)+1]; 
   elseif j(i)<6
       rtBins{i}=[j(i) j(i)+1.5]; 
   else
       rtBins{i}=[j(i) j(i)+3]; 
   end
   nowRTbins=rtBins{i};
   temp1_seq2=dataset.event_RT_trial1InSeq{n_trials_away};
   b=temp1_seq2(dataset.realrtpair_seq2{n_trials_away}==1);
   dimall_event=dataset.alldim_rtchanges_event{n_trials_away};
   temp_RTdiffsevent=dimall_event(b>=nowRTbins(1,1) & b<nowRTbins(1,2));
   dataToFit=temp_RTdiffsevent;
   firstRTs=b(b>=nowRTbins(1,1) & b<nowRTbins(1,2));
   histo_nbins=[-4*12.4245:0.2510:4*12.4245];
   [~,mif]=nanmin(abs(try_curr_rts-nanmean(firstRTs)));
   pred=randsample(bin_centers(bins),length(firstRTs),true,rt_pdf_outs(mif,:));
   prediction=firstRTs-pred;
   [~,~,fHist]=plotHist(dataToFit,prediction,histo_nbins,['Checking rate fits for first RT ' num2str(nanmean(firstRTs))],'Change in RT (sec)');
   pause;
   close(fHist);
end

end

function [lambdas,rtBins]=fitRateTermToLateReaches_local(dataset,suppressOutput)

j=0:0.5:8.5;
for i=1:length(j)
   if j(i)<3
       rtBins{i}=[j(i) j(i)+1]; 
   elseif j(i)<6
       rtBins{i}=[j(i) j(i)+1.5]; 
   else
       rtBins{i}=[j(i) j(i)+3]; 
   end
   lambdas(i)=getTau(dataset,rtBins{i},suppressOutput);
   disp('Fitting rate');
end

end

function bestTau=getTau(dataset,rtBins,suppressOutput)

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
% [histo_nbins,RTchangeHist]=plotHist(temp_RTdiffs1,temp_RTdiffsevent,histo_nbins,['CDF RT change for RTs ' num2str(dataset.nInSequence(i)-1) ' trials later, RTs less than ' num2str(rtBins(j,2))],'Change in RT');

temp=dataset.event_name;
temp(regexp(temp,'_'))=' ';
a=temp1(dataset.realrtpair_seq1{i}==1);
b=temp1_seq2(dataset.realrtpair_seq2{i}==1);
if suppressOutput==0
    [histo_nbins,RTpdf]=plotHist(a(a>=rtBins(j,1) & a<=rtBins(j,2)),b(b>=rtBins(j,1) & b<=rtBins(j,2)),histo_nbins,['Reaction Times fx of ' temp ': trial 1 (black) vs ' num2str(dataset.nInSequence(i)-1) ' later (magenta)'],'RT (sec)');
end

bestTau=fitSingleTau(temp_RTdiffsevent,b(b>=rtBins(j,1) & b<=rtBins(j,2)),[rtBins(2)-maxTrialLength rtBins(2)],histo_nbins,suppressOutput);
if suppressOutput==0
    disp(bestTau);
    pause;
    close all;
end

% [slope_early,slope_late]=fitTwoTaus(dimall_all(a>=rtBins(j,1) & a<rtBins(j,2)),histo_nbins,[]);
% disp(slope_early);
% disp(slope_late);
% x=0:0.001:20;
% y=randsample(x,length(temp_RTdiffs1),true,-(exp(-abs(slope_late).*x)-exp(-abs(slope_early).*x)));
% plotHist(a(a>=rtBins(j,1) & a<=rtBins(j,2))-dimall_all(a>=rtBins(j,1) & a<rtBins(j,2)),y,histo_nbins,'Fitting tau','RT (sec)');
% 
% [slope_early,slope_late]=fitTwoTaus(dimall_event(b>=rtBins(j,1) & b<rtBins(j,2)),histo_nbins,[]);
% disp(slope_early);
% disp(slope_late);
% x=0:0.001:20;
% y=randsample(x,length(temp_RTdiffs1),true,-(exp(-abs(slope_late).*x)-exp(-abs(slope_early).*x)));
% plotHist(b(b>=rtBins(j,1) & b<=rtBins(j,2))-dimall_event(b>=rtBins(j,1) & b<rtBins(j,2)),y,histo_nbins,'Fitting tau','RT (sec)');

end

function [slope_early,slope_late]=fitTwoTaus(dataToFit,histo_nbins,forceEarlySlope)

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
if ~isempty(forceEarlySlope)
    slope_early=forceEarlySlope;
end
slope_late=nanmean(diff(log_n_data(lateIndBegin:end))./diff(x_log(lateIndBegin:end)));

figure();
plot(x_log,log_n_data,'Color','k');
hold on;
plot(x_log(~isnan(log_n_data)),slope_early.*x_log(~isnan(log_n_data))+nanmean(log_n_data(1:earlyIndEnd)),'Color','r');
plot(x_log(~isnan(log_n_data)),slope_late.*x_log(~isnan(log_n_data))+log_n_data(end)-slope_late.*x_log(end),'Color','g');

end

function bestTau=fitSingleTau(dataToFit,firstRTs,fitRange,histo_nbins,suppressOutput)

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
%     differ=nansum(abs(n_data(xhist>=fitRange(1) & xhist<=fitRange(2))-n(xhist>=fitRange(1) & xhist<=fitRange(2))));
    differ=nansum(sqrt((n_data(xhist>=fitRange(1) & xhist<=fitRange(2))-n(xhist>=fitRange(1) & xhist<=fitRange(2))).^2));
    differs(i)=differ;
end
[~,mi]=nanmin(differs);
bestTau=tauRange(mi);
y=randsample(x,length(firstRTs),true,exp(-bestTau.*x));
if suppressOutput==0
    plotHist(dataToFit,firstRTs-y,histo_nbins,'Fitting tau','Change in RT (sec)');
end

end
    
function [x_backup,returnThis,f]=plotHist(data1,data2,bins,tit,xlab)

useLogForY=false;

[n,x]=histcounts(data1,bins);
x_backup=x;
[n,x]=cityscape_hist(n,x);
f=figure();
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

function alignedOut=alignToEvent(event,withinRange_inds,whatToAlign,howMuchToTake)

event=event(:,withinRange_inds(1):withinRange_inds(2));

eventInds=nan(1,size(event,1));
for i=1:size(event,1)
    if any(event(i,:)>0.5,2)
        f=find(event(i,:)>0.5,1,'first');
        eventInds(i)=f;
    end
end
eventInds=eventInds+withinRange_inds(1)-1;

j=1;
alignedOut=nan(sum(~isnan(eventInds)),howMuchToTake);
for i=1:size(event,1)
    if ~isnan(eventInds(i))
        alignedOut(j,:)=whatToAlign(i,eventInds(i):eventInds(i)+howMuchToTake-1);
        j=j+1;
    end
end

end