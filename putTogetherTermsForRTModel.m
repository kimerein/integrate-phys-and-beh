function em=putTogetherTermsForRTModel(dataset,wrtDataset,alltbt,reachType,useCuedRateData)

% bins=0:0.035:9;
% % histo_nbins=[-1*12.4245:0.2510:1*12.4245];
% histo_nbins=[-9-0.035/2:0.035:9];

bins=0:0.1:9;
% histo_nbins=[-1*12.4245:0.2510:1*12.4245];
% histo_nbins=[-9-0.1/2:0.1:9];
histo_nbins=[-9:0.1:9];

% Plot real data
i=1;
temp1_seq2=dataset.event_RT_trial1InSeq{i};
x2=temp1_seq2(dataset.realrtpair_seq2{i}==1);
y2=dataset.alldim_rtchanges_event{i};
nBinsPerDim={bins; histo_nbins};
realData=plotRTbyRTchangePlot(x2,y2,nBinsPerDim,'Real data');

% Plot reference data
% Reference input
temp1_seq2=wrtDataset.event_RT_trial1InSeq{1};
x2=temp1_seq2(wrtDataset.realrtpair_seq2{1}==1);
y2=wrtDataset.alldim_rtchanges_event{1};
% All trials as reference
% temp1_seq2=dataset.allTrialsSequence_RT_trial1InSeq{1};
% x2=temp1_seq2(dataset.realrtpair_seq1{1}==1);
% y2=dataset.alldim_rtchanges_allTrialsSequence{1};
nBinsPerDim={bins; histo_nbins};
referenceData=plotRTbyRTchangePlot(x2,y2,nBinsPerDim,'Reference data');

% Real data rate term fit
rateFits=getRateFits(alltbt,dataset,reachType,useCuedRateData);
temp=makeRateTerm(nBinsPerDim,rateFits,dataset);
prediction.rate_term=temp;

% Reference data rate term fit
rateFits=getRateFits(alltbt,wrtDataset,'miss',useCuedRateData);
temp=makeRateTerm(nBinsPerDim,rateFits,wrtDataset);
prediction.reference_rate_term=temp;

% RTM term
% temp=fitRTMtermToOtherReaches_local(dataset,nBinsPerDim,bins,1);
% prediction.rtm_term=temp;

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

function rateFits=getRateFits(alltbt,dataset,reachType,useCuedRateData)

% Non-cued approach
withinRange=[2 6];  % take reaches in this range, time is in seconds from onset of cue
% choose a start time for withinRange well after cue
nIndsToTake=200;
timeStep=mode(diff(nanmean(alltbt.times,1)));
cueInd=find(nanmean(alltbt.cueZone_onVoff,1)>0,1,'first');
withinRange_inds=[cueInd+ceil(withinRange(1)/timeStep) cueInd+ceil(withinRange(2)/timeStep)];
temp=alltbt.reachBatch_drop_reachStarts;
touchAligned=alignToEvent(temp,withinRange_inds,alltbt.all_reachBatch,nIndsToTake);
temp=alltbt.reachBatch_miss_reachStarts;
reachNoTouchAligned=alignToEvent(temp,withinRange_inds,alltbt.all_reachBatch,nIndsToTake);
temp=alltbt.reachBatch_success_reachStarts;
successTouchAligned=alignToEvent(temp,withinRange_inds,alltbt.all_reachBatch,nIndsToTake);
temp=alltbt.reachBatch_drop_reachStarts+alltbt.reachBatch_success_reachStarts;
anyTouchAligned=alignToEvent(temp,withinRange_inds,alltbt.all_reachBatch,nIndsToTake);
% Zero out first reach inds
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
legend({'drop','miss','success','any touch'});
rateFromNoncuedData.x=0:timeStep:timeStep*(length(nanmean(touchAligned,1))-1);
if strcmp(reachType,'drop')
    rateFromNoncuedData.y=nanmean(touchAligned,1)./timeStep;
elseif strcmp(reachType,'miss')
    rateFromNoncuedData.y=nanmean(reachNoTouchAligned,1)./timeStep;
elseif strcmp(reachType,'success')
    rateFromNoncuedData.y=nanmean(successTouchAligned,1)./timeStep;
elseif strcmp(reachType,'anyTouch')
    rateFromNoncuedData.y=nanmean(anyTouchAligned,1)./timeStep;
end
rateFromNoncuedData.y(3:-1:1)=rateFromNoncuedData.y(3);
% Fit rates using cued approach
if useCuedRateData==true
    [temp1,temp2]=fitRateTermToLateReaches_local(dataset,1); % last argument is suppressOutput
    rateFromCuedData.x=nan(1,length(temp2));
    for i=1:length(temp2)
        rateFromCuedData.x(i)=nanmean(temp2{i});
    end
    rateFromCuedData.y=temp1;
    rateFits=rateFromCuedData;
end
% Generate rate term
rateFits=rateFromNoncuedData;

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

function [rt_pdf_out,bins]=update_RPE_term(rt_pdf,bins,curr_rt,alpha,suppressOutput)

% RT pdf gives RT cdf
rt_cdf = @(rt_pdf) accumulateDistribution(rt_pdf)./nanmax(accumulateDistribution(rt_pdf));

% RT cdf gives 1-cdf
rt_rpe = @(rt_pdf) 1-rt_cdf(rt_pdf);

% Plot
bin_centers = @(bins) nanmean([bins(1:end-1); bins(2:end)],1);
if suppressOutput==false
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
    figure(); plot(bin_centers(bins),rpe_at_touch,'Color','k'); title('RPE at touch');
    figure(); plot(bin_centers(bins),value_update(rt_pdf,curr_touch,rpe_at_touch,alpha),'Color','k'); title('Updated value estimate');
    figure(); plot(bin_centers(bins),rt_pdf_out,'Color','k'); title('Update RT PDF');
end

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

function rate_term=makeRateTerm(nBinsPerDim,rateFits,dataset)

rt_pdf = @(rts,bins) histcounts(rts,bins)./nansum(histcounts(rts,bins));
bin_centers = @(bins) nanmean([bins(1:end-1); bins(2:end)],1);

i=1;
rts=dataset.event_RT_trial1InSeq{i};

bins=nBinsPerDim{1};
bins=bins-(bins(2)-bins(1))/2;
bins=[bins bins(end)+(bins(2)-bins(1))/2];

% As a function of different current reaction times
try_curr_rts=bin_centers(bins);
rt_change_pdfs=nan(length(try_curr_rts),length(bin_centers(unique([-fliplr(bins) bins]))));
rt_pdf_outs=update_rate_term(bins,rateFits,dataset);
for i=1:length(try_curr_rts)
    rt_pdf_out=rt_pdf_outs(i,:);
    [rt_change_pdfs(i,:),rt_change_bins]=getRTchange_fromCurrAndUpdatedRTpdfs(rt_pdf(rts,bins),rt_pdf_out,bins,true,[bins(i) bins(i+1)]);
end
rt_change_pdfs=rt_change_pdfs./nansum(nansum(rt_change_pdfs));

rate_term.rt_change_pdfs=rt_change_pdfs;
rate_term.rt_change_bins=rt_change_bins;
rate_term.rt_pdf_outs=rt_pdf_outs;

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

function rt_pdf_out=update_rate_term(bins,fits,dataset)

bin_centers = @(bins) nanmean([bins(1:end-1); bins(2:end)],1);
histo_nbins=[-4*12.4245:0.2510:4*12.4245];

% i=1;
% temp1_seq2=dataset.event_RT_trial1InSeq{i};
% b=temp1_seq2(dataset.realrtpair_seq2{i}==1);

bincents=bin_centers(bins);
x=0:0.001:20;
rt_pdf_out=nan(length(bincents),length(bincents));
for i=1:length(bincents)
    [~,mi]=min(abs(fits.x-bincents(i)));
    currLambda=fits.y(mi);
%     firstRTs=b(b>=bins(i) & b<=bins(i+1));
    y=randsample(x,10000,true,exp(-currLambda.*x));
%     [n,x]=histcounts(firstRTs-y,histo_nbins); % for RT changes
    [n,x]=histcounts(y,bins);
    n=n./nansum(n);
    rt_pdf_out(i,:)=n;
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
    differ=nansum(abs(n_data(xhist>=fitRange(1) & xhist<=fitRange(2))-n(xhist>=fitRange(1) & xhist<=fitRange(2))));
    differs(i)=differ;
end
[~,mi]=nanmin(differs);
bestTau=tauRange(mi);
y=randsample(x,length(firstRTs),true,exp(-bestTau.*x));
if suppressOutput==0
    plotHist(dataToFit,firstRTs-y,histo_nbins,'Fitting tau','RT (sec)');
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