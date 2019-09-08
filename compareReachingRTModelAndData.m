function [prediction,a,b,lse,err,n,best_predict,differ]=compareReachingRTModelAndData(data,dataForModel,fits,behaviorEvent,ITI,suppressOutput)

% Model
% Rate-based term is fit from data
% Depends on behavioral event and ITI
% RPE-based term is built from principle
% Need to fit alpha, the learning rate

% Choose parameters for model build and analysis
bins=0:0.035:9;
rts=dataForModel.allTrialsSequence_RT_trial1InSeq{1};
rts_firstTrialInSequence=data.event_RT_trial1InSeq{1};
shortestCuedRT=0.125; % in seconds, accounts for sensory detection + motor planning/execution delay
preemptCue_secBeforeCue=0.1; % earliest that mouse can detect preemptive cue, in seconds before cue

% Mouse is not doing a cued reach if reaches before shortestCuedRT, e.g., 125 ms
rts(rts<=shortestCuedRT)=nan;
rts_firstTrialInSequence(rts_firstTrialInSequence<=shortestCuedRT)=nan;

% For default (i.e., comparing consecutive trials), do not change these
sameRTforEachTrial=false; % if true, will assume same reaction time for all subsequent trials, otherwise will sample randomly from current rt pdf
n_update_steps=1; % How many update steps (e.g., trials)


% Rate term
rt_pdf = @(rts,bins) histcounts(rts,bins)./nansum(histcounts(rts,bins));
bin_centers = @(bins) nanmean([bins(1:end-1); bins(2:end)],1);
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
    [rt_change_pdfs(i,:),rt_change_bins]=getRTchange_fromCurrAndUpdatedRTpdfs(rt_pdf(rts_firstTrialInSequence,bins),rt_pdf_out,bins,true,[bins(i) bins(i+1)]); % if last argument is empty, will plot change distribution for all current reaction times
end
% Plot imagesc(rt_pdf_outs) to see how pdf of RTs on next trial is expected
% to change, given a change in rate only
prediction.rate_term.rt_change_pdfs=rt_change_pdfs;
prediction.rate_term.rt_change_bins=rt_change_bins;
prediction.rate_term.rt_pdf_outs=rt_pdf_outs;


% RPE update
% For consecutive trials, this is just the RT pdf pushing subsequent RTs
% forward in time
% Note that positive RPE update is proportional to 1-cdf of RTs
smearSize=60;
temp=zeros(size(rt_pdf_outs));
% Use the followed code (commented out) to convince yourself that the
% RPE update to the PDF of reaction times is just rt_pdf
% curr_cdf=accumulateDistribution(rt_pdf(rts_firstTrialInSequence,bins));
% curr_updates=1-curr_cdf;
% curr_updates=[0 diff(curr_updates)];
curr_updates=rt_pdf(rts_firstTrialInSequence,bins);
curr_updates(curr_updates==0)=0.00001;
curr_updates=curr_updates./nansum(curr_updates);
for i=1:size(temp,1)
    if i-(smearSize-1)<1
        useInd=1;
    else
        useInd=i-(smearSize-1);
    end
    temp(i,useInd:i)=curr_updates(i)/smearSize;
end
temp(temp(1:end)<=0)=0.00001;
temp=temp./nansum(nansum(temp));
rt_pdf_outs=temp;
prediction.rpe_only_consec_update.rt_pdf_outs=rt_pdf_outs;

rt_change_pdfs=nan(length(try_curr_rts),length(bin_centers(unique([-fliplr(bins) bins]))));
for i=1:length(try_curr_rts)
    rt_pdf_out=rt_pdf_outs(i,:);
    [rt_change_pdfs(i,:),rt_change_bins]=getRTchange_fromCurrAndUpdatedRTpdfs(rt_pdf(rts_firstTrialInSequence,bins),rt_pdf_out,bins,true,[bins(i) bins(i+1)]); % if last argument is empty, will plot change distribution for all current reaction times
end
prediction.rpe_only_consec_update.rt_change_pdfs=rt_change_pdfs;
prediction.rpe_only_consec_update.rt_change_bins=rt_change_bins;

% prediction.rpe_only_consec_update.rt_change_pdfs=zeros(size(rt_change_pdfs));
% % prediction.rpe_only_consec_update.rt_change_pdfs(:,find(bin_centers(rt_change_bins)>0,1,'first'))=rt_pdf(rts_firstTrialInSequence,bins);
% for i=1:length(try_curr_rts)
%     f=find(bin_centers(rt_change_bins)>0 & bin_centers(rt_change_bins)<try_curr_rts(i));
%     if length(f)>smearSize % 160
%         f=f(1:smearSize);
%     end
%     if isempty(f)
%         continue
%     end
%     temp=rt_pdf(rts_firstTrialInSequence,bins);
%     prediction.rpe_only_consec_update.rt_change_pdfs(i,f)=repmat(temp(i)./length(f),1,length(f));
% end
% prediction.rpe_only_consec_update.rt_change_bins=rt_change_bins;

% RPE effects over multiple trials
alpha=1;
try_curr_rts=bin_centers(bins);
rt_change_pdfs=nan(length(try_curr_rts),length(bin_centers(unique([-fliplr(bins) bins]))));
rt_pdf_outs=nan(length(try_curr_rts),length(bin_centers(bins)));
for i=1:length(try_curr_rts)
    rt_pdf_out=prediction.rpe_only_consec_update.rt_pdf_outs(i,:);
    if sameRTforEachTrial==true
        for j=1:n_update_steps-1 % if RT is same for each trial
            rt_pdf_out=update_RPE_term(rt_pdf_out,bins,try_curr_rts(i),alpha,true);
        end
    else
        for j=1:n_update_steps-1 % if RT is randomly sampled from RT distribution for each trial
            if j==1
                rt_pdf_out=update_RPE_term(rt_pdf_out,bins,try_curr_rts(i),alpha,true);
            else
                rt_pdf_out=update_RPE_term(rt_pdf_out,bins,randsample(try_curr_rts,1,true,rt_pdf(rts,bins)),alpha,true);
            end
        end
    end
    rt_pdf_outs(i,:)=rt_pdf_out;
    [rt_change_pdfs(i,:),rt_change_bins]=getRTchange_fromCurrAndUpdatedRTpdfs(rt_pdf(rts,bins),rt_pdf_out,bins,true,[bins(i) bins(i+1)]);
end
prediction.rpe_term.rt_change_pdfs=rt_change_pdfs;
prediction.rpe_term.rt_change_bins=rt_change_bins;
prediction.rpe_term.rt_pdf_outs=rt_pdf_outs;

% Fit terms to data
i=n_update_steps;
temp1_seq2=data.event_RT_trial1InSeq{i};
temp1_seq2=temp1_seq2(data.realrtpair_seq2{i}==1);
temp2_seq2=data.event_RT_trialiInSeq{i};
temp2_seq2=temp2_seq2(data.realrtpair_seq2{i}==1);
rtchanges_seq2=temp1_seq2-temp2_seq2;
[n,x,y]=make2Dpdf(temp1_seq2,rtchanges_seq2,{bin_centers(bins), bin_centers(rt_change_bins)});
% Smooth all before fitting
K=ones(7);
n=conv2(n,K,'same');
prediction.rate_term.rt_change_pdfs=conv2(prediction.rate_term.rt_change_pdfs,K,'same');
prediction.rpe_term.rt_change_pdfs=conv2(prediction.rpe_term.rt_change_pdfs,K,'same');
prediction.rpe_only_consec_update.rt_change_pdfs=conv2(prediction.rpe_only_consec_update.rt_change_pdfs,K,'same');
prediction.rate_term.rt_change_pdfs=prediction.rate_term.rt_change_pdfs./nansum(nansum(prediction.rate_term.rt_change_pdfs));
prediction.rpe_term.rt_change_pdfs=prediction.rpe_term.rt_change_pdfs./nansum(nansum(prediction.rpe_term.rt_change_pdfs));
prediction.rpe_only_consec_update.rt_change_pdfs=prediction.rpe_only_consec_update.rt_change_pdfs./nansum(nansum(prediction.rpe_only_consec_update.rt_change_pdfs));
n=n./nansum(nansum(n));
% fitMask=repmat(bin_centers(rt_change_bins)>0,length(bin_centers(bins)),1);
fitMask=ones(size(prediction.rpe_term.rt_change_pdfs));
% fitMask(bin_centers(bins)>2,:)=0;
[a,b,lse]=getFitCoefficients(n,prediction.rate_term.rt_change_pdfs,prediction.rpe_term.rt_change_pdfs,suppressOutput,fitMask);
% rateMask=repmat(bin_centers(rt_change_bins)>-0.5,length(bin_centers(bins)),1);
% rateMask(bin_centers(bins)<0.5,:)=1;
% rpeMask=repmat(bin_centers(rt_change_bins)<=0 | bin_centers(rt_change_bins)>=0.25,length(bin_centers(bins)),1);
% [a,b,lse]=fitRateThenRPE(n,prediction.rate_term.rt_change_pdfs,prediction.rpe_only_consec_update.rt_change_pdfs,rateMask,rpeMask);
if suppressOutput==0
    disp('ratio of rpe to rate term');
    disp(b/a);
    disp('rate term coefficient');
    disp(a);
    disp('rpe term coefficient');
    disp(b);
    disp('least squares error');
    disp(lse);
end

best_predict=a.*prediction.rate_term.rt_change_pdfs+b.*prediction.rpe_only_consec_update.rt_change_pdfs;
best_predict=best_predict./nansum(nansum(best_predict));

if suppressOutput==0
    plot_deltaRT_as_func_of_RT(x,y,[n best_predict],'Real Data then Fit to Data');
end

K=ones(10);
smoothMat=conv2([n best_predict],K,'same');
if suppressOutput==0
    plot_deltaRT_as_func_of_RT(x,y,smoothMat,'Real Data then Fit to Data -- smoothed');
end

% Plot difference between real data and fit
% differ=conv2(n,K,'same')-conv2(best_predict,K,'same');
differ=n-best_predict;
err=nansum(nansum(differ));
differ=differ-nanmin(nanmin(differ));
if suppressOutput==0
    plot_deltaRT_as_func_of_RT(x,y,differ,'Real Data Minus Fit');
end



% Preemptive reaching term
prediction.preempt_term.preempt_details=fitPreemptiveReaching(x,y,differ,0,rts,bins,shortestCuedRT,preemptCue_secBeforeCue);
process_a=gampdf(x,prediction.preempt_term.preempt_details.shape,prediction.preempt_term.preempt_details.rate);
process_b=-prediction.preempt_term.preempt_details.subtractPDF;
for i=1:length(try_curr_rts)
%     rt_pdf_outs(i,:)=prediction.preempt_term.preempt_details.process_a_coefficients(i).*process_a+prediction.preempt_term.preempt_details.process_b_coefficients(i).*process_b;
    rt_pdf_outs(i,:)=prediction.preempt_term.preempt_details.process_a_coefficients(i).*process_a;
end
rt_pdf_outs=rt_pdf_outs-nanmin(rt_pdf_outs(1:end));
rt_pdf_outs(rt_pdf_outs<=0)=0.000000001;
prediction.preempt_term.rt_pdf_outs=rt_pdf_outs;
rt_change_pdfs=nan(length(try_curr_rts),length(bin_centers(unique([-fliplr(bins) bins]))));
for i=1:length(try_curr_rts)
    rt_pdf_out=rt_pdf_outs(i,:);
    [rt_change_pdfs(i,:),rt_change_bins]=getRTchange_fromCurrAndUpdatedRTpdfs(rt_pdf(data.event_RT_trial1InSeq{1},bins),rt_pdf_out,bins,true,[bins(i) bins(i+1)]); % if last argument is empty, will plot change distribution for all current reaction times
end
prediction.preempt_term.rt_change_pdfs=rt_change_pdfs;
prediction.preempt_term.rt_change_bins=rt_change_bins;
if suppressOutput==0
    plot_deltaRT_as_func_of_RT(x,y,rt_change_pdfs,'Preemptive reaching RT change PDF');
end

end

function preempt=fitPreemptiveReaching(x,y,differ,suppressOutput,rts,bins,shortestCuedRT,preemptCue_secBeforeCue)

rt_pdf = @(rts,bins) histcounts(rts,bins)./nansum(histcounts(rts,bins));
bin_centers = @(bins) nanmean([bins(1:end-1); bins(2:end)],1);

% Model preemptive reaching as a switch between 2 strategies (or states)
% Cued reaching and pre-emptive reaching
% Fit emission and transition probabilities
% Cued reaching (> shortestCuedRT) follow the RT pdf
% Preemptive reaching + counting as a reaction time combines the PDF of
% emitting a preemptive reach and the requirement of RT (must be >0 sec, step function)
% Fit preemptive reaching as gamma PDF
preemptFalloff=getPDFat(x,y,differ,[-shortestCuedRT shortestCuedRT]);
preemptFalloff=preemptFalloff-preemptFalloff(end);
[best_a,best_rate,best_shape,fitqualval,subtractPDF,x_corrected,fitToPreempt,nadded]=fitGamma_w_transition(x,preemptFalloff,rt_pdf(rts,bins),preemptCue_secBeforeCue,bins);
if suppressOutput==0
    figure(); 
    preemptFalloff=preemptFalloff./nansum(abs(preemptFalloff));
    plot(x,preemptFalloff,'Color','k');
    hold on; 
    plot(x_corrected,fitToPreempt,'Color','r');
    legend({'Data', 'Fit'});
    title('Fit preemptive reaching');
end
preempt.scale=best_a;
preempt.rate=best_rate;
preempt.shape=best_shape;
preempt.updateToRTpdf=fitToPreempt(1:end-nadded);
preempt.subtractPDF=subtractPDF(1:end-nadded);
preempt.process_a_coefficients=best_a.*gampdf(x_corrected,best_shape,best_rate);
preempt.process_a_coefficients=preempt.process_a_coefficients(1:end-nadded);
preempt.process_b_coefficients=(1-best_a).*subtractPDF;
preempt.process_b_coefficients=preempt.process_b_coefficients(1:end-nadded);

end

function [best_a,best_rate,best_shape,fitqualval,subtractPDF,waittimes,y_w_subtract,n_to_add]=fitGamma_w_transition(waittimes,countwaits,subtractPDF,preemptCue_secBeforeCue,bins)

% Range chosen from previous runs of grid search
% actually, according to Matlab paramaterization of gamma pdf, rate is 1/try_rate_params
try_rate_params=[0.01:0.01:0.2]; % from previous run, best rate is 0.1
try_shape_params=[0.5:0.01:1.5]; % from previous run, best shape is 1
try_a=[0.01:0.01:0.1]; % from previous run, best a is 0.02

% Preemptive reach cue is BEFORE real cue
bin_step=bins(2)-bins(1);
n_to_add=floor(preemptCue_secBeforeCue./bin_step);
waittimes=[waittimes-n_to_add*bin_step:bin_step:waittimes-bin_step waittimes];
backup_waittimes=waittimes;
waittimes=waittimes-min(waittimes); % start at preempt cue
% countwaits=[zeros(n_to_add,1); countwaits];
subtractPDF=[zeros(1,n_to_add) subtractPDF];
countwaits=[countwaits; zeros(n_to_add,1)];
% Am in preemptive reaching cue reference frame now (in terms of time)

fitqual=nan(length(try_a),length(try_rate_params),length(try_shape_params));
countwaits=countwaits./nansum(abs(countwaits));
subtractPDF=subtractPDF./nansum(subtractPDF);
for i=1:length(try_a)
    curr_a=try_a(i);
    for j=1:length(try_rate_params)
        curr_rate_param=try_rate_params(j);
        for k=1:length(try_shape_params)
            curr_shape_param=try_shape_params(k);
            y=gampdf(waittimes,curr_shape_param,curr_rate_param);
            y_w_subtract=curr_a.*y-(1-curr_a).*subtractPDF; % Probability of cued reaching is 1 - probability of preemptive reaching
            y_w_subtract=(y_w_subtract./nanmax(y_w_subtract)).*nanmax(countwaits);
            fitqual(i,j,k)=nansum((countwaits-y_w_subtract').^2);
        end
    end
end

[fitqualval,fitind]=nanmin(fitqual(1:end)); 
[a,b,c]=ind2sub(size(fitqual),fitind);
best_a=try_a(a);
disp(['try_a min and max ' num2str([nanmin(try_a) nanmax(try_a)])]);
disp('best a');
disp(best_a);
best_rate=try_rate_params(b);
disp(['try_rate min and max ' num2str([nanmin(try_rate_params) nanmax(try_rate_params)])]);
disp('best rate');
disp(best_rate);
best_shape=try_shape_params(c);
disp(['try_shape min and max ' num2str([nanmin(try_shape_params) nanmax(try_shape_params)])]);
disp('best shape');
disp(best_shape);

curr_a=best_a; 
y=gampdf(waittimes,best_shape,best_rate);
y_w_subtract=curr_a.*y-(1-curr_a).*subtractPDF;

end

function [best_rate,best_shape,fitqualval]=fitGamma(waittimes,countwaits)

% Range chosen from previous runs of grid search
try_rate_params=[0.01:0.01:0.5]; % actually, according to Matlab paramaterization of gamma pdf, rate is 1/try_rate_params
try_shape_params=[0.9:0.01:1.5];

fitqual=nan(length(try_rate_params),length(try_shape_params));
countwaits=countwaits-nanmin(countwaits);
countwaits=countwaits./nansum(countwaits);
for j=1:length(try_rate_params)
    curr_rate_param=try_rate_params(j);
    for k=1:length(try_shape_params)
        curr_shape_param=try_shape_params(k);
        y=gampdf(waittimes,curr_shape_param,curr_rate_param);
        y=y./nansum(y); % this is a pdf
        fitqual(j,k)=nansum((countwaits-y').^2);
    end
end

[fitqualval,fitind]=nanmin(fitqual(1:end)); 
[b,c]=ind2sub(size(fitqual),fitind);
best_rate=try_rate_params(b);
best_shape=try_shape_params(c);

end

function plot_deltaRT_as_func_of_RT(x,y,pdf2D,tit)

figure();
% Note intensity is plotted on log scale!
imagesc(x,y,log(pdf2D'));
set(gca,'YDir','normal');
xlabel('Reaction time (seconds)');
ylabel('Change in reaction time (seconds)');
title(tit);

end

function slice=getPDFat(x,y,pdf2D,y_bin)

slice=nanmean(pdf2D(:,y>=y_bin(1) & y<=y_bin(2)),2);

end

function [a,b,fitqualval]=fitRateThenRPE(data,rateTerm,rpeTerm,rateMask,rpeMask)

% a is coefficient for rateTerm
% b is coefficient for rpeTerm

try_a=[0.001:0.05:10];
cost=nan(1,length(try_a));
masked_data=data;
masked_data(rateMask==1)=0;
masked_rateTerm=rateTerm;
masked_rateTerm(rateMask==1)=0;
masked_data=masked_data./nansum(nansum(masked_data));
for i=1:length(try_a)
    curr_a=try_a(i);
    cost(i)=nansum(nansum((masked_data-(curr_a.*masked_rateTerm)./nansum(nansum(curr_a.*masked_rateTerm))).^2));
end
[fitqualval1,fitind]=nanmin(cost);
a=try_a(fitind);

try_b=[0.001:0.05:10];
cost=nan(1,length(try_b));
masked_data=data;
masked_data=masked_data-a.*rateTerm./nansum(nansum(a.*rateTerm));
masked_data(rpeMask==1)=0;
masked_rpeTerm=rpeTerm;
masked_rpeTerm(rpeMask==1)=0;
masked_data=masked_data./nansum(nansum(masked_data));
for i=1:length(try_b)
    curr_b=try_b(i);
    cost(i)=nansum(nansum((masked_data-(curr_b.*masked_rpeTerm)./nansum(nansum(curr_b.*masked_rpeTerm))).^2));
end
[fitqualval2,fitind]=nanmin(cost);
b=try_b(fitind);
fitqualval=fitqualval1+fitqualval2;

end

function [a,b,fitqualval,coefficients]=getFitCoefficients(data,rateTerm,rpeTerm,suppressOutput,fitMask)

% a is coefficient for rateTerm
% b is coefficient for rpeTerm

data=data./nansum(nansum(data));

try_a=[0.001:0.05:3];
try_b=[0.001:0.05:3];
% Looking for ratio of a and b that best fits data
cost=nan(length(try_a),length(try_b));
for i=1:length(try_a)
    curr_a=try_a(i);
    if mod(i,100)==0
        disp(i);
    end
    for j=1:length(try_b)
        curr_b=try_b(j);
        temp=(data-(curr_a.*rateTerm+curr_b.*rpeTerm)./nansum(nansum(curr_a.*rateTerm+curr_b.*rpeTerm))).^2;
        cost(i,j)=nansum(nansum(temp(fitMask==1)));
    end
end
% Best ratio of a and b defines a line in cost space
% Find this line
nSamples=20;
minPairs=nan(nSamples,2);
maxout=nanmax(cost(1:end));
maxout=maxout+1000;
backup_cost=cost;
for i=1:nSamples
    [~,fitind]=nanmin(cost(1:end)); 
    [a_i,b_i]=ind2sub(size(cost),fitind);
    minPairs(i,:)=[a_i b_i];
    cost(a_i,b_i)=maxout;
end
cost=backup_cost;

if suppressOutput==0
    figure();
    imagesc(try_a,try_b,cost');
    set(gca,'YDir','normal');
    xlabel('Rate coefficient');
    ylabel('RPE coefficient');
end
    
coefficients=polyfit(try_a(minPairs(:,1)),try_b(minPairs(:,2)),1);
if suppressOutput==0
    figure();
    scatter(try_a(minPairs(:,1)),try_b(minPairs(:,2)));
    xFit=linspace(min(try_a(minPairs(:,1))),max(try_a(minPairs(:,1))),1000);
    yFit=polyval(coefficients,xFit);
    hold on;
    plot(xFit,yFit);
end

[fitqualval,fitind]=nanmin(cost(1:end)); 
[a_i,b_i]=ind2sub(size(cost),fitind);
a=try_a(a_i);
b=try_b(b_i);

end

function [n,x,y]=make2Dpdf(x1,y1,nBinsPerDim)

% make 2D histogram
if size(x1,2)>1
    % make column vector
    x1=x1';
    y1=y1';
end
if length(nBinsPerDim)>1
    [n,bin_c]=hist3([x1 y1],'Ctrs',nBinsPerDim);
else
    [n,bin_c]=hist3([x1 y1],[nBinsPerDim nBinsPerDim]);
end

% Normalize 
n=n./nansum(nansum(n,1),2);
x=bin_c{1};
y=bin_c{2};

end

function invertRTchange_to_rt_pdf_update(rt_change_pdf,)

% For each reaction time, sample changes from rt_change_pdf
% Apply these changes to curr_rt_pdf to get updated_rt_pdf
sample_n=30000;
bin_centers = @(bins) nanmean([bins(1:end-1); bins(2:end)],1);
rt_changes=randsample(

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