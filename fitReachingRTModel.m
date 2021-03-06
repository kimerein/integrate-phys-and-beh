function prediction=fitReachingRTModel(dataset,alltbt,metadata,trialTypes,rateTermParams,behaviorEvent,ITI)

bins=0:0.035:9;
rts=dataset.allTrialsSequence_RT_trial1InSeq{1};
temp=dataset.event_RT_trial1InSeq{1};
temp=temp(dataset.realrtpair_seq2{1}==1);
temp2_seq2=dataset.event_RT_trialiInSeq{1};
temp2_seq2=temp2_seq2(dataset.realrtpair_seq2{1}==1);
rtchanges_seq2=temp-temp2_seq2;
% rts=temp;
smoothSize=1;
sameRTforEachTrial=false; % if true, will assume same reaction time for all subsequent trials, otherwise will sample randomly from current rt pdf
n_update_steps=1; % How many update steps (e.g., trials) 

% Rate-based term and RPE-based term
% Each of these makes a prediction about how the reaction time distribution
% will update, given a particular input distribution and behavioral event

% Get prediction for change in the reaction time distribution, as a
% function of the input reaction time distribution and current reaction time



% Rate-based term

% Get effects of behavioral event on rate of reaching
fits=fitGammaAandBeta(alltbt,rateTermParams,true);
if isempty(rateTermParams)
    save('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\mouse summary data 20190726\gammafits.mat','fits');
end

rt_pdf = @(rts,bins) histcounts(rts,bins)./nansum(histcounts(rts,bins));
bin_centers = @(bins) nanmean([bins(1:end-1); bins(2:end)],1);

% As a function of different current reaction times
try_curr_rts=bin_centers(bins);
rt_change_pdfs=nan(length(try_curr_rts),length(bin_centers(unique([-fliplr(bins) bins]))));
rt_pdf_outs=nan(length(try_curr_rts),length(bin_centers(bins)));
for i=1:length(try_curr_rts)
    rt_pdf_out=rt_pdf(rts,bins);
    if sameRTforEachTrial==true
        for j=1:n_update_steps % if RT is same for each trial
            rt_pdf_out=update_rate_term(rt_pdf_out,bins,try_curr_rts(i),ITI,behaviorEvent,true,fits);
        end
    else
        for j=1:n_update_steps % if RT is randomly sampled from RT distribution for each trial
            if j==1
                rt_pdf_out=update_rate_term(rt_pdf_out,bins,try_curr_rts(i),ITI,behaviorEvent,true,fits); 
            else
                rt_pdf_out=update_rate_term(rt_pdf_out,bins,randsample(try_curr_rts,1,true,rt_pdf_out),ITI,behaviorEvent,true,fits);
            end
        end
    end
    rt_pdf_outs(i,:)=rt_pdf_out;
    [rt_change_pdfs(i,:),rt_change_bins]=getRTchange_fromCurrAndUpdatedRTpdfs(rt_pdf(rts,bins),rt_pdf_out,bins,true,[bins(i) bins(i+1)]); % if last argument is empty, will plot change distribution for all current reaction times
end

prediction.rate_term.rt_change_pdfs=rt_change_pdfs;
prediction.rate_term.rt_change_bins=rt_change_bins;
prediction.rate_term.rt_pdf_outs=rt_pdf_outs;

% Plot result
K=ones(smoothSize);
smoothMat=conv2(rt_change_pdfs',K,'same');
figure();
imagesc(try_curr_rts,bin_centers(rt_change_bins),log(smoothMat));
set(gca,'YDir','normal');
title(['Rate term -- Change in RT as a function of RT (smoothed)']);
xlabel('Reaction time trial 1');
ylabel('Change in reaction times');



% RPE-based term

% Learning rate, alpha
% alpha=0.03;
% alpha=0.055;
alpha=0.5;
% alpha=0.034;
% alpha=0;

temp1=dataset.allTrialsSequence_RT_trial1InSeq{1};
temp1=temp1(dataset.realrtpair_seq1{1}==1);
temp1_seq2=dataset.event_RT_trial1InSeq{1};
temp1_seq2=temp1_seq2(dataset.realrtpair_seq2{1}==1);
temp2_seq2=dataset.event_RT_trialiInSeq{1};
temp2_seq2=temp2_seq2(dataset.realrtpair_seq2{1}==1);
rtchanges_seq2=temp1_seq2-temp2_seq2;
fitAlpha(temp1_seq2,temp1,temp1_seq2,rtchanges_seq2,{bin_centers(bins), bin_centers(rt_change_bins)});
% fitAlphaAfterRateSubtract(temp1_seq2,temp1,temp1_seq2,rtchanges_seq2,{bin_centers(bins), bin_centers(rt_change_bins)},prediction);

% RT pdf as proxy for value prediction
rt_pdf = @(rts,bins) histcounts(rts,bins)./nansum(histcounts(rts,bins));
bin_centers = @(bins) nanmean([bins(1:end-1); bins(2:end)],1);

% As a function of different current reaction times
try_curr_rts=bin_centers(bins);
rt_change_pdfs=nan(length(try_curr_rts),length(bin_centers(unique([-fliplr(bins) bins]))));
rt_pdf_outs=nan(length(try_curr_rts),length(bin_centers(bins)));
for i=1:length(try_curr_rts)
    rt_pdf_out=rt_pdf(rts,bins);
    if sameRTforEachTrial==true
        for j=1:n_update_steps % if RT is same for each trial
            rt_pdf_out=update_RPE_term(rt_pdf_out,bins,try_curr_rts(i),alpha,true);
        end
    else
        for j=1:n_update_steps % if RT is randomly sampled from RT distribution for each trial
            if j==1 
                rt_pdf_out=update_RPE_term(rt_pdf_out,bins,try_curr_rts(i),alpha,true);
            else
                rt_pdf_out=update_RPE_term(rt_pdf_out,bins,randsample(try_curr_rts,1,true,rt_pdf_out),alpha,true);
            end
        end
    end
    rt_pdf_outs(i,:)=rt_pdf_out;
    [rt_change_pdfs(i,:),rt_change_bins]=getRTchange_fromCurrAndUpdatedRTpdfs(rt_pdf(rts,bins),rt_pdf_out,bins,true,[bins(i) bins(i+1)]);
end

prediction.rpe_term.rt_change_pdfs=rt_change_pdfs;
prediction.rpe_term.rt_change_bins=rt_change_bins;
prediction.rpe_term.rt_pdf_outs=rt_pdf_outs;

% Plot result
K=ones(smoothSize);
smoothMat=conv2(rt_change_pdfs',K,'same');
figure();
imagesc(try_curr_rts,bin_centers(rt_change_bins),log(smoothMat));
set(gca,'YDir','normal');
title(['RPE term -- Change in RT as a function of RT (smoothed)']);
xlabel('Reaction time trial 1');
ylabel('Change in reaction times');

% Plot how RPE update accumulates over multiple trials, given different RTs
% on different trials
% If assume non-changing underlying RT pdf 
% Note that this is the same as increasing alpha, the learning rate
% update_matrix=rt_pdf_outs-repmat(rt_pdf(rts,bins),size(rt_pdf_outs,1),1);
% n_update_steps=50;
% rt_change_pdfs=nan(length(try_curr_rts),length(bin_centers(unique([-fliplr(bins) bins]))));
% for i=1:length(try_curr_rts)
%     [rt_change_pdfs(i,:),rt_change_bins]=getRTchange_fromCurrAndUpdatedRTpdfs(rt_pdf(rts,bins),rt_pdf(rts,bins)+n_update_steps.*update_matrix(i,:),bins,true,[bins(i) bins(i+1)]);
% end
% if smoothSize>1
%     K=ones(smoothSize);
%     smoothMat=conv2(rt_change_pdfs',K,'same');
%     figure();
%     imagesc(try_curr_rts,bin_centers(rt_change_bins),smoothMat);
%     set(gca,'YDir','normal');
%     title(['RPE term -- multiple update']);
%     xlabel('Reaction time trial 1');
%     ylabel('Change in reaction times');
% end
return

% Plot sum
A=1;
B=1;
K=ones(smoothSize);
smoothMat=conv2(A.*prediction.rpe_term.rt_change_pdfs'+B.*prediction.rate_term.rt_change_pdfs',K,'same');
figure();
imagesc(try_curr_rts,bin_centers(rt_change_bins),smoothMat);
set(gca,'YDir','normal');
title(['Sum of RPE and rate terms -- Change in RT as a function of RT (smoothed)']);
xlabel('Reaction time trial 1');
ylabel('Change in reaction times');



% Feed RPE-updated RT pdf into rate adjustment
try_curr_rts=bin_centers(bins);
rt_change_pdfs=nan(length(try_curr_rts),length(bin_centers(unique([-fliplr(bins) bins]))));
for i=1:length(try_curr_rts)
    rt_pdf_out=update_RPE_term(rt_pdf(rts,bins),bins,try_curr_rts(i),alpha,true);
    rt_pdf_out=update_rate_term(rt_pdf_out,bins,try_curr_rts(i),ITI,behaviorEvent,true,fits);
    [rt_change_pdfs(i,:),rt_change_bins]=getRTchange_fromCurrAndUpdatedRTpdfs(rt_pdf(rts,bins),rt_pdf_out,bins,true,[bins(i) bins(i+1)]);
end

prediction.rpe_then_rate_term.rt_change_pdfs=rt_change_pdfs;
prediction.rpe_then_rate_term.rt_change_bins=rt_change_bins;

% Plot RPE adjustment then rate adjustment, combined
K=ones(smoothSize);
smoothMat=conv2(prediction.rpe_then_rate_term.rt_change_pdfs',K,'same');
figure();
imagesc(try_curr_rts,bin_centers(rt_change_bins),smoothMat);
set(gca,'YDir','normal');
title(['RPE adjustment before rate adjustment -- Change in RT as a function of RT (smoothed)']);
xlabel('Reaction time trial 1');
ylabel('Change in reaction times');



% Plot after removing regression to the mean
rt_change_pdfs=prediction.rpe_term.rt_change_pdfs;
rt_change_bins=prediction.rpe_term.rt_change_bins;
% [diff2Dhist,x,y,n,n2]=removeMeanRegression(rts,rts,randsample(bin_centers(rt_change_bins),length(rts),true,rt_pdf(rts,bins)*rt_change_pdfs),{bin_centers(bins), bin_centers(rt_change_bins)});
[diff2Dhist,x,y,n,n2]=removeMeanRegression(rts,temp1_seq2,rtchanges_seq2,{bin_centers(bins), bin_centers(rt_change_bins)});
prediction.regress_to_mean.rt_change_pdfs=n;
prediction.regress_to_mean.rt_change_bins=rt_change_bins;
prediction.from_theory_rpe.rt_change_pdfs=n2;
prediction.from_theory_rpe.rt_change_bins=rt_change_bins;
prediction.rt_bins=x;
prediction.rt_change_bins=y;
K=ones(smoothSize);
smoothMat=conv2(diff2Dhist'-nanmin(nanmin(diff2Dhist')),K,'same');
figure();
imagesc(x,y,log(smoothMat(1:end-1,1:end-1)));
set(gca,'YDir','normal');
title(['Rate/RPE term -- regression to mean removed']);
xlabel('Reaction time trial 1');
ylabel('Change in reaction times');

end

function fitAlphaAfterRateSubtract(rts1,rts2,actual_first_rts,actual_rt_changes,binsFor2Dhist,ratePred)
% pass in a dataset for quantifying RPE learning rate, alpha
% alpha is the learning rate per trial
% may vary as a function of learning stage, etc.
% but here just approximate a single learning rate 

RPE_slice_at=[0 0.25]; % bin of delta_rts where expect to see RPE
nonRPE_slice_at=[-14 -0.5]; % bin of delta_rts where do not expect to see RPE

% fit regression to the mean to the slowed reaction time pairs
% fit putative RPE component to the faster reaction time pairs

% get what we expect from regression to the mean
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

% make heatmap of regression to the mean
% make heatmap of real, presumably RPE-containing data
[diff_n,x,y,regToMean,realData]=compareWithHeatmaps(first_rt,delta_rts,actual_first_rts,actual_rt_changes,binsFor2Dhist);

regToMean=ratePred.rate_term.rt_change_pdfs;
regToMean=regToMean./nansum(nansum(regToMean));

figure();
temp2=realData-regToMean;
temp2=temp2-nanmin(nanmin(temp2));
imagesc(x,y,log(temp2'));
set(gca,'YDir','normal');

nonRPE_regToMean=getPDFat(x,y,regToMean,nonRPE_slice_at);
nonRPE_realData=getPDFat(x,y,realData,nonRPE_slice_at);
scaleFac=fitNonRPEtoData(nonRPE_regToMean,nonRPE_realData);

figure();
plot(x,scaleFac.*nonRPE_regToMean,'Color','k');
hold on;
plot(x,nonRPE_realData,'Color','r');
legend({'Regress to mean','Real data'});
title('Scale non-RPE components');

RPE_regToMean=getPDFat(x,y,regToMean,RPE_slice_at);
RPE_realData=getPDFat(x,y,realData,RPE_slice_at);

figure();
plot(x,scaleFac.*RPE_regToMean,'Color','k');
hold on;
plot(x,RPE_realData,'Color','r');
legend({'Regress to mean','Real data'});
title('RPE components, scaled');

pause;

% Fit RT pdf to difference, in order to get alpha
% Be sure that RT pdf sums to one
putRPE=RPE_realData-scaleFac.*RPE_regToMean; 
[n,x_centers]=histcounts(actual_first_rts,0:x(2)-x(1):9);
n=n./nansum(n);
scaleFacForAlpha=fitNonRPEtoData(n,putRPE');
figure();
plot(x,scaleFacForAlpha.*n,'Color','k');
hold on;
plot(x,putRPE,'Color','r');
legend({'Fit RT pdf','Putative RPE'});
title('Guess alpha');
disp('alpha');
disp(scaleFacForAlpha);

end

function fitAlpha(rts1,rts2,actual_first_rts,actual_rt_changes,binsFor2Dhist)
% pass in a dataset for quantifying RPE learning rate, alpha
% alpha is the learning rate per trial
% may vary as a function of learning stage, etc.
% but here just approximate a single learning rate 

RPE_slice_at=[0 0.25]; % bin of delta_rts where expect to see RPE
nonRPE_slice_at=[-14 -0.5]; % bin of delta_rts where do not expect to see RPE

% fit regression to the mean to the slowed reaction time pairs
% fit putative RPE component to the faster reaction time pairs

% get what we expect from regression to the mean
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

% make heatmap of regression to the mean
% make heatmap of real, presumably RPE-containing data
[diff_n,x,y,regToMean,realData]=compareWithHeatmaps(first_rt,delta_rts,actual_first_rts,actual_rt_changes,binsFor2Dhist);

nonRPE_regToMean=getPDFat(x,y,regToMean,nonRPE_slice_at);
nonRPE_realData=getPDFat(x,y,realData,nonRPE_slice_at);
scaleFac=fitNonRPEtoData(nonRPE_regToMean,nonRPE_realData);

figure();
plot(x,scaleFac.*nonRPE_regToMean,'Color','k');
hold on;
plot(x,nonRPE_realData,'Color','r');
legend({'Regress to mean','Real data'});
title('Scale non-RPE components');

RPE_regToMean=getPDFat(x,y,regToMean,RPE_slice_at);
RPE_realData=getPDFat(x,y,realData,RPE_slice_at);

figure();
plot(x,scaleFac.*RPE_regToMean,'Color','k');
hold on;
plot(x,RPE_realData,'Color','r');
legend({'Regress to mean','Real data'});
title('RPE components, scaled');

% Fit RT pdf to difference, in order to get alpha
% Be sure that RT pdf sums to one
putRPE=RPE_realData-scaleFac.*RPE_regToMean;
[n,x_centers]=histcounts(actual_first_rts,0:x(2)-x(1):9);
n=n./nansum(n);
scaleFacForAlpha=fitNonRPEtoData(n,putRPE');
figure();
plot(x,scaleFacForAlpha.*n,'Color','k');
hold on;
plot(x,putRPE,'Color','r');
legend({'Fit RT pdf','Putative RPE'});
title('Guess alpha');
disp('alpha');
disp(scaleFacForAlpha);

end

function scaleFac=fitNonRPEtoData(regToMean,realData)

tryRange=0.001:0.001:50;
% tryRange=0.00001:0.00001:1;
cost=nan(1,length(tryRange));
for i=1:length(tryRange)
    currRange=tryRange(i);
    cost(i)=nansum((currRange.*regToMean-realData).^2); % least sum of squares
%     cost(i)=nansum(abs(currRange.*regToMean-realData)); 
end
[~,mi]=nanmin(cost);
scaleFac=tryRange(mi);
disp('scale fac mins and max');
disp(num2str([min(tryRange) max(tryRange)]));
disp('scale fac');
disp(scaleFac);

end

function slice=getPDFat(x,y,pdf2D,y_bin)

slice=nanmean(pdf2D(:,y>=y_bin(1) & y<=y_bin(2)),2);

end

function [diff2Dhist,x,y,n,n2]=removeMeanRegression(rts,actual_first_rts,actual_rt_changes,binsFor2Dhist)

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
% Heatmap comparison
[diff2Dhist,x,y,n,n2]=compareWithHeatmaps(first_rt,delta_rts,actual_first_rts,actual_rt_changes,binsFor2Dhist);

end

function [diff_n,x,y,n,n2]=compareWithHeatmaps(x1,y1,x2,y2,nBinsPerDim)

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
title(['Regression to the mean only']);
xlabel('Reaction time trial 1');
ylabel('Change in reaction times');

figure();
imagesc(bin_c{1},bin_c{2},log(n2'));
set(gca,'YDir','normal');
title(['From RT pairs with Rate/RPE']);
xlabel('Reaction time trial 1');
ylabel('Change in reaction times');

figure();
imagesc(bin_c{1},bin_c{2},log(diff_n'-nanmin(nanmin(diff_n'))));
set(gca,'YDir','normal');
title(['RT pairs with Rate/RPE minus regression to the mean']);
xlabel('Reaction time trial 1');
ylabel('Change in reaction times');

x=bin_c{1};
y=bin_c{2};

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

function update_rt_pdf=update_rate_term(rt_pdf,bins,curr_rt,curr_ITI,curr_event,suppressOutput,fits)

% After a long time (ITI), effect of touching pellet on previous trial is
% primarily a change in rate of reaching
% Model this as a multiplication of pdf(reaching) by pdf(wait time)

% Mapping from curr_rt to rate parameter of pdf(wait time):
% suppression of reaching wears off after the mouse touches pellet, 
% so b, the inverse rate parameter, will decrease as time delay from touch
% increases, i.e., shorter curr_rt means higher b

bin_centers = @(bins) nanmean([bins(1:end-1); bins(2:end)],1);
pdf_wait = @(bins,a,b) gampdf(bin_centers(bins),a,b);

% get a,b for pdf(wait time) from curr_rt
[a,b]=getGammaAandBeta(curr_rt,curr_ITI,curr_event,fits);

normFac=calibrateRateUpdate(fits,pdf_wait);

update_rt_pdf = normFac*rt_pdf.*pdf_wait(bins,a,b);

if suppressOutput==false
    figure(); 
    plot(bin_centers(bins),rt_pdf,'Color','k'); 
    hold on;
    plot(bin_centers(bins),update_rt_pdf,'Color','r');
    legend({'before rate update','after rate update'});
    title('Update rate based on pellet touch');     
end

end

function normFac=calibrateRateUpdate(fits,pdf_wait)

% after a very long time, should be no effect of previous reach on the
% current rate of reaching 
% (note RPE is considered elsewhere)
% set normalization factors such that this is true

bins=0:0.01:1000;
curr_ITI=nanmax(fits.nSecondsAfterTouch);
[a,b]=getGammaAandBeta(0,curr_ITI,'drop',fits);
rt_pdf=exp(-bins(2:end));
update_rt_pdf = rt_pdf.*pdf_wait(bins,a,b);
% update_rt_pdf should have integral that matches integral of rt_pdf
% because we have waited a very long time
% need to multiply integral of update_rt_pdf by normFac to make this true
rt_pdf_integral=nansum(rt_pdf);
update_rt_pdf_integral=nansum(update_rt_pdf);
normFac=rt_pdf_integral/update_rt_pdf_integral;

% figure();
% plot(rt_pdf,'Color','k');
% hold on;
% plot(update_rt_pdf,'Color','r');
% plot(normFac*update_rt_pdf,'Color','b');

end

function [a,b]=getGammaAandBeta(curr_rt,curr_ITI,curr_event,fits)

% curr_rt is the reaction time on the current trial
% I am trying to predict the reaction time on the next trial
% this will depend on the ITI between the current and next trials
% that inter-trial interval (ITI) is given in curr_ITI
% the time until the next cue is the ITI minus the curr_rt

timeTilCue=curr_ITI-curr_rt;
if timeTilCue<=fits.nSecondsAfterTouch(1)
    f=1;
else
    f=find(fits.nSecondsAfterTouch(1:end-1)<timeTilCue & fits.nSecondsAfterTouch(2:end)>=timeTilCue,1,'first');
end
if isempty(f)
    error('Rate effects not defined for this reaction time / ITI in getGammaAandBeta');
end

switch curr_event
    case 'drop'
        % shape, a
        a=fits.shapeParam_dropPellet(f);
        % inverse rate, b
        b=fits.rateParam_dropPellet(f);
    case 'miss'
        % shape, a
        a=fits.shapeParam_reachButNoTouch(f);
        % inverse rate, b
        b=fits.rateParam_reachButNoTouch(f);
    case 'success'
        % shape, a
        a=fits.shapeParam_successTouch(f);
        % inverse rate, b
        b=fits.rateParam_successTouch(f);
    otherwise
        error('Do not recognize curr_event in getGammaAandBeta');
end

end

function fits=fitGammaAandBeta(alltbt,rateTermParams,suppressOutput)

withinRange=[2 6];  % take reaches in this range, time is in seconds from onset of cue
% choose a start time for withinRange well after cue
nIndsToTake=200;
smoothFitParams=false; % will smooth fitted params over time
useFitsToFits=true; % will use fits for how params change over time

if isempty(rateTermParams)
    timeStep=mode(diff(nanmean(alltbt.times,1)));
    cueInd=find(nanmean(alltbt.cueZone_onVoff,1)>0,1,'first');
    
    withinRange_inds=[cueInd+ceil(withinRange(1)/timeStep) cueInd+ceil(withinRange(2)/timeStep)];

    temp=alltbt.reachBatch_drop_reachStarts;
    touchAligned=alignToEvent(temp,withinRange_inds,alltbt.all_reachBatch,nIndsToTake);
    
    temp=alltbt.reachBatch_miss_reachStarts;
    reachNoTouchAligned=alignToEvent(temp,withinRange_inds,alltbt.all_reachBatch,nIndsToTake);
    
    temp=alltbt.reachBatch_success_reachStarts;
    successTouchAligned=alignToEvent(temp,withinRange_inds,alltbt.all_reachBatch,nIndsToTake);
    
    % Zero out first reach inds
    normForTouch=size(touchAligned,1);
    touchAligned(:,1)=0;
    normForNoTouch=size(reachNoTouchAligned,1);
    reachNoTouchAligned(:,1)=0;
    normForSuccessTouch=size(successTouchAligned,1);
    successTouchAligned(:,1)=0;
    
    % Have already grouped reaches into batches, so each batch includes
    % multiple "paw still out" attempts
    % Animal's max rate of executing reach batches is about 2 Hz
    % This makes sense, because animal take about 150 ms to get paw from the
    % perch to the pellet zone, and each reach batch begins when paw returns to
    % perch
    % 0.5/0.035=14.2857 for a single bin mode
    scaleUp=14.2857;
    
    % Get distributions of waits
    indsAfter=1:floor(0.1/timeStep):nIndsToTake;
    waitTimes=nan(length(indsAfter),size(touchAligned,1));
    waitTimes_noTouch=nan(length(indsAfter),size(reachNoTouchAligned,1));
    waitTimes_successTouch=nan(length(indsAfter),size(successTouchAligned,1));
    count_waits=nan(length(indsAfter),length(0:timeStep:size(touchAligned,2)*timeStep)-1);
    count_waits_noTouch=nan(length(indsAfter),length(0:timeStep:size(reachNoTouchAligned,2)*timeStep)-1);
    count_waits_successTouch=nan(length(indsAfter),length(0:timeStep:size(successTouchAligned,2)*timeStep)-1);
    nSecondsAfterTouch=indsAfter.*timeStep;
    for i=1:length(indsAfter)
        waitTimes(i,:)=getWaitTimes(touchAligned(:,indsAfter(i):end),timeStep);
        n=histcounts(waitTimes(i,:),0:timeStep:size(touchAligned,2)*timeStep);
        count_waits(i,:)=(n./normForTouch)*scaleUp;
        waitTimes_noTouch(i,:)=getWaitTimes(reachNoTouchAligned(:,indsAfter(i):end),timeStep);
        n=histcounts(waitTimes_noTouch(i,:),0:timeStep:size(reachNoTouchAligned,2)*timeStep);
        count_waits_noTouch(i,:)=(n./normForNoTouch)*scaleUp;
        waitTimes_successTouch(i,:)=getWaitTimes(successTouchAligned(:,indsAfter(i):end),timeStep);
        n=histcounts(waitTimes_successTouch(i,:),0:timeStep:size(successTouchAligned,2)*timeStep);
        count_waits_successTouch(i,:)=(n./normForSuccessTouch)*scaleUp;
    end
    
    % count_waits is the expected distribution of wait times before reach, given
    % no cue, aligned to pellet touch
    % count_waits_noTouch is the expected distribution of wait times before
    % reach, given no cue, aligned to reach without pellet touch
    % assume that, after reach but no pellet touch on previous trial, the
    % distribution of RTs on next trial is already multiplied by
    % count_waits_noTouch
    % so, to account for rate change associated with pellet touch, need to
    % divide out count_waits_noTouch and multiply by count_waits, or,
    % equivalently, multiply by count_waits/count_waits_noTouch
    
    % fit gamma distribution to each count_waits and count_waits_noTouch
    % plot change in gamma rate parameter as the time delay from touch increases
    % then can use this relationship to predict gamma rate parameter at some
    % arbitrary time delay
    
    [nSecondsAfterTouch,rateParam,shapeParam]=getRateAndShapeForWaiting(nSecondsAfterTouch,count_waits,timeStep);
    [nSecondsAfterTouch,rateParam_noTouch,shapeParam_noTouch]=getRateAndShapeForWaiting(nSecondsAfterTouch,count_waits_noTouch,timeStep);
    [nSecondsAfterTouch,rateParam_successTouch,shapeParam_successTouch]=getRateAndShapeForWaiting(nSecondsAfterTouch,count_waits_successTouch,timeStep);
    
    % note that later time points are harder to fit, if curve is flattening
    % out (also losing points as take longer nSecondsAfterTouch);
    % therefore, take trend from fits to lower nSecondsAfterTouch
    % this makes an assumption that there is a consistent progression of rate
    % change as effects of reach wear off, but this makes sense
    % assume that after touch_no_pellet and touch_with_pellet, rate param
    % converges to the same baseline value
    convergeToRate=45.25; % from previous grid search
    convergeToShape=0.755; % from previous grid search
    % max trial length is 9 seconds, last possible reach in this batch is at 6
    % s, so after about 3 s, quality of curves degrades
    convergeAfterThisTime=3.885; % in seconds
    shapeParam(nSecondsAfterTouch>convergeAfterThisTime)=convergeToShape;
    shapeParam_noTouch(nSecondsAfterTouch>convergeAfterThisTime)=convergeToShape;
    shapeParam_successTouch(nSecondsAfterTouch>convergeAfterThisTime)=convergeToShape;
    rateParam(nSecondsAfterTouch>convergeAfterThisTime)=convergeToRate;
    rateParam_noTouch(nSecondsAfterTouch>convergeAfterThisTime)=convergeToRate;
    rateParam_successTouch(nSecondsAfterTouch>convergeAfterThisTime)=convergeToRate;
else
    nSecondsAfterTouch=rateTermParams.nSecondsAfterTouch;
    shapeParam=rateTermParams.shapeParam_dropPellet;
    shapeParam_noTouch=rateTermParams.shapeParam_reachButNoTouch;
    shapeParam_successTouch=rateTermParams.shapeParam_successTouch;
    rateParam=rateTermParams.rateParam_dropPellet;
    rateParam_noTouch=rateTermParams.rateParam_reachButNoTouch;
    rateParam_successTouch=rateTermParams.rateParam_successTouch;
end

if useFitsToFits==true
    x=nSecondsAfterTouch(1):0.01:30; 
    y=lognpdf(x,0.6,1);
    fitted=0.755-0.45*y;
    if suppressOutput==false
        figure();
        plot(nSecondsAfterTouch,shapeParam_noTouch,'Color','k');
        hold on;
        plot(x,fitted,'Color','r');
        legend({'raw data','fit'});
        title('Shape parameter after reach but no touch');
    end
    shapeParam_noTouch=fitted;
    
    y=lognpdf(x,0,1);
    fitted=0.755+0.45*y;
    if suppressOutput==false
        figure();
        plot(nSecondsAfterTouch,shapeParam,'Color','k');
        hold on;
        plot(x,fitted,'Color','r');
        legend({'raw data','fit'});
        title('Shape parameter after touch and drop');
    end
    shapeParam=fitted;
    
    y1=lognpdf(x,0,1); y2=lognpdf(x,0.1,0.25); y3=lognpdf(x,1.1,0.25);
    fitted=0.755+0.2*y1-0.15*y2+0.8*y3;
    if suppressOutput==false
        figure();
        plot(nSecondsAfterTouch,shapeParam_successTouch,'Color','k');
        hold on;
        plot(x,fitted,'Color','r');
        legend({'raw data','fit'});
        title('Shape parameter after touch with eat');
    end
    shapeParam_successTouch=fitted;
    
    y=lognpdf(x,0,15);
    fitted=45.25-90*y;
    if suppressOutput==false
        figure();
        plot(nSecondsAfterTouch,rateParam_noTouch,'Color','k');
        hold on;
        plot(x,fitted,'Color','r');
        legend({'raw data','fit'});
        title('Inverse rate parameter after reach but no touch');
    end
    rateParam_noTouch=fitted;
    
    y=lognpdf(x,0.4,1.2);
    fitted=45.25-45*y;
    if suppressOutput==false
        figure();
        plot(nSecondsAfterTouch,rateParam,'Color','k');
        hold on;
        plot(x,fitted,'Color','r');
        legend({'raw data','fit'});
        title('Inverse rate parameter after touch and drop');
    end
    rateParam=fitted;
    
    y2=lognpdf(x,0.4,1.2); y1=lognpdf(x,0.25,0.35); y3=lognpdf(x,1.5,1); y4=lognpdf(x,0,20);
    fitted=45.25-10*y2+85*y1-150*y3+150*y4;
    if suppressOutput==false
        figure();
        plot(nSecondsAfterTouch,rateParam_successTouch,'Color','k');
        hold on;
        plot(x,fitted,'Color','r');
        legend({'raw data','fit'});
        title('Inverse rate parameter after touch with eat');
    end
    rateParam_successTouch=fitted;
    nSecondsAfterTouch=x;
end
    
% reach but no pellet touch
% shape parameter is something like
% x=0.1:0.1:7; y=lognpdf(x,0.6,1);
% plot(x,0.755-0.45*y,'Color','g');

% touch pellet, drop pellet
% shape parameter is something like
% x=0.1:0.1:7; y=lognpdf(x,0,1);
% plot(x,0.755+0.45*y,'Color','b');

% grabs and eats pellet
% shape parameter is something like
% x=0.1:0.1:7; y1=lognpdf(x,0,1); y2=lognpdf(x,0.1,0.25); y3=lognpdf(x,1.1,0.25);
% plot(x,0.755+0.2*y1-0.15*y2+0.8*y3,'Color','r');

% reach but no pellet touch
% inverse rate parameter is something like 
% x=0.1:0.1:7; y=lognpdf(x,0,15);
% plot(x,45.25-90*y);

% touch pellet, drop pellet
% inverse rate parameter is something like
% x=0.1:0.1:7; y=lognpdf(x,0.4,1.2);
% plot(x,45.25-45*y);

% grabs and eats pellet
% inverse rate parameter is something like
% x=0.1:0.1:7; y2=lognpdf(x,0.4,1.2); y1=lognpdf(x,0.25,0.35); y3=lognpdf(x,1.5,1); y4=lognpdf(x,0,20);
% plot(x,45.25-10*y2+85*y1-150*y3+150*y4);

% smooth fits? note that this will cause slight time shift
if smoothFitParams==true
    smoothBin=5;
    shapeParam=smooth(shapeParam,smoothBin);
    shapeParam_noTouch=smooth(shapeParam_noTouch,smoothBin);
    shapeParam_successTouch=smooth(shapeParam_successTouch,smoothBin);
    rateParam=smooth(rateParam,smoothBin);
    rateParam_noTouch=smooth(rateParam_noTouch,smoothBin);
    rateParam_successTouch=smooth(rateParam_successTouch,smoothBin);
end

if suppressOutput==false
    figure();
    plot(nSecondsAfterTouch,shapeParam,'Color','k');
    hold on;
    plot(nSecondsAfterTouch,shapeParam_noTouch,'Color','r');
    plot(nSecondsAfterTouch,shapeParam_successTouch,'Color','g');
    legend({'touches pellet','reach but no pellet touch','grabs and eats pellet'});
    title('Shape parameter of gamma fit to waits before next reach');
    
    figure();
    plot(nSecondsAfterTouch,rateParam,'Color','k');
    hold on;
    plot(nSecondsAfterTouch,rateParam_noTouch,'Color','r');
    plot(nSecondsAfterTouch,rateParam_successTouch,'Color','g');
    legend({'touches pellet','reach but no pellet touch','grabs and eats pellet'});
    title('Inverse rate parameter of gamma fit to waits before next reach');
end

fits.nSecondsAfterTouch=nSecondsAfterTouch;
fits.shapeParam_dropPellet=shapeParam;
fits.shapeParam_reachButNoTouch=shapeParam_noTouch;
fits.shapeParam_successTouch=shapeParam_successTouch;
fits.rateParam_dropPellet=rateParam;
fits.rateParam_reachButNoTouch=rateParam_noTouch;
fits.rateParam_successTouch=rateParam_successTouch;

end

function [nSecondsAfterTouch,touch_rates,touch_shapes,touch_scales]=getRateAndShapeForWaiting(nSecondsAfterTouch,count_waits,timeStep)

touch_scales=nan(1,length(nSecondsAfterTouch));
touch_rates=nan(1,length(nSecondsAfterTouch));
touch_shapes=nan(1,length(nSecondsAfterTouch));
disp(['Counting to ' num2str(length(nSecondsAfterTouch))]);
for i=1:length(nSecondsAfterTouch)
    if mod(i,10)==0
        disp(i);
    end
    currTimeDelay=nSecondsAfterTouch(i);
    f=find(count_waits(i,end:-1:1)~=0,1,'first');
    [touch_scales(i),touch_rates(i),touch_shapes(i)]=fitGamma(timeStep:timeStep:(length(count_waits(i,:))-(f-1))*timeStep,count_waits(i,1:length(count_waits(i,:))-(f-1)));
end

end

function [best_scale,best_rate,best_shape]=fitGamma(waittimes,countwaits)

% try different scales for count of wait times
% try different rates for gamma pdf
% try a small range of shape parameters
% but note that if the shape parameter is 1, this is the exponential pdf,
% which represents the probability distribution of the time between events
% in a Poisson process (also exponential random variable is memoryless -- pure rate)
% however, reach batches may be more grouped in time or more periodic than pure Poisson?

% measure fit quality as the sum of squared differences divided by the
% integral of count of wait times

% try_scaling=[0:0.1:1];
try_scaling=[1]; % chosen from previous runs of grid search
try_rate_params=[0.1:0.1:100]; % actually, according to Matlab paramaterization of gamma pdf, rate is 1/try_rate_params
% try_shape_params=[0.5 0.75 1 1.25 1.5 2 3 5 10];
try_shape_params=[0.5:0.01:1.3];
% try_shape_params=[0.8817]; % chosen from previous runs of grid search

fitqual=nan(length(try_scaling),length(try_rate_params),length(try_shape_params));
normFactor=nansum(countwaits);
for i=1:length(try_scaling)
    curr_scale=try_scaling(i);
    for j=1:length(try_rate_params)
        curr_rate_param=try_rate_params(j);
        for k=1:length(try_shape_params)
            curr_shape_param=try_shape_params(k);
            y=curr_scale.*gampdf(waittimes,curr_shape_param,curr_rate_param);
            fitqual(i,j,k)=nansum((countwaits-y).^2)./normFactor;
        end
    end
end

[fitqualval,fitind]=nanmin(fitqual(1:end)); 
[a,b,c]=ind2sub(size(fitqual),fitind);
best_scale=try_scaling(a);
best_rate=try_rate_params(b);
best_shape=try_shape_params(c);

% figure();
% y=gampdf(waittimes,try_shape_params(c),try_rate_params(b));
% plot(waittimes,try_scaling(a).*countwaits,'Color','k'); 
% hold on;
% plot(waittimes,y,'Color','r');

end

function waitTimes=getWaitTimes(reaches,timeStepForInd)

waitTimes=nan(1,size(reaches,1));
for i=1:size(reaches,1)
    f=find(reaches(i,:)>0.5);
    if ~isempty(f)
        waitTimes(i)=nanmin(f).*timeStepForInd;
    end
end

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

function [x_backup,counts1,counts2]=plotHist(data1,data2,bins,tit,xlab,c1,c2)

[n,x]=histcounts(data1,bins);
x_backup=x;
[n,x]=cityscape_hist(n,x);
figure();
plot(x,n./nansum(n),'Color',c1);
counts1=n;
xlabel(xlab);
ylabel('Count');
title(tit);

[n,x]=histcounts(data2,x_backup);
[n,x]=cityscape_hist(n,x);
hold on;
plot(x,n./nansum(n),'Color',c2);
counts2=n;
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


