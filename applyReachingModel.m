function applyReachingModel(use_dataset_dir,fits_file,preempt_file,alltbt,out,metadata,dprime_range,event_type,override_file,n_update_steps)

% Load reaching data
a=load([use_dataset_dir '\pdf.mat']);
dataset=a.dataset;

% Filter reaching data for sessions with d-prime in dprime_range
[dprime_alltbt]=filtTbt(alltbt,out,'dprime',dprime_range,metadata);

% Get base-subtracted cued reaching
cued_process=getCuedProcess(dprime_alltbt,false);

% Load preemptive reaching term
a=load([preempt_file]);
rt_pdf_for_preemptModel=a.rt_pdf_for_preemptModel;
preemptive_process=a.preemptive_process;
ss=a.ss;

% Load fits for rate term
a=load([fits_file]);
fits=a.fits;

% Load prediction for preemptive reaching
if ~isempty(override_file)
    a=load(override_file);
    override_prediction=a.prediction;
else
    override_prediction=[];
end

% Get model fits
compareReachingRTModelAndData(dataset.realDistributions,dataset.realDistributions,fits,event_type,9.5,0,rt_pdf_for_preemptModel,preemptive_process,ss,cued_process,override_prediction,n_update_steps);

end

function cued_process=getCuedProcess(dprime_alltbt,suppressOutput)

% shortestCuedRT=0.125; % in seconds, accounts for sensory detection + motor planning/execution delay
shortestCuedRT=0.2; % in seconds, accounts for sensory detection + motor planning/execution delay
removeShortRTs=1;

bin_centers = @(bins) nanmean([bins(1:end-1); bins(2:end)],1);
bins=0:0.035:9;

% Use baseline-subtracted reaching for cue pdf
rts=getBaseSubtractedRTs(dprime_alltbt,'cueZone_onVoff','all_reachBatch');
% Mouse is not doing a cued reach if reaches before shortestCuedRT, e.g., 125 ms
if removeShortRTs==1
    rts(rts<=shortestCuedRT)=nan;
end
n=histcounts(rts,bins);
cued_process.rt_pdf=n;
cued_process.rt_pdf=cued_process.rt_pdf./nansum(cued_process.rt_pdf);
cued_process.rt_pdf(cued_process.rt_pdf==0)=0.002;
cued_process.rt_pdf_outs=repmat(cued_process.rt_pdf,length(bin_centers(bins)),1);
try_curr_rts=bin_centers(bins);
for i=1:length(try_curr_rts)
    rt_pdf_out=cued_process.rt_pdf_outs(i,:);
    if all(rt_pdf_out==0)
        rt_pdf_out(1:end)=0.000000001;
    end
    [rt_change_pdfs(i,:),rt_change_bins]=getRTchange_fromCurrAndUpdatedRTpdfs(cued_process.rt_pdf,rt_pdf_out,bins,true,[bins(i) bins(i+1)]);
end
cued_process.rt_change_pdfs=rt_change_pdfs;
cued_process.rt_change_bins=rt_change_bins;
cued_process.rt_change_pdfs=cued_process.rt_change_pdfs./nansum(nansum(abs(cued_process.rt_change_pdfs)));

if suppressOutput==false
    figure();
    plot(bin_centers(bins),cued_process.rt_pdf,'Color','k');
    xlabel('RT (sec)');
    ylabel('Count');
    title('Baseline reaching-subtracted cued process');
    
    figure(); 
    imagesc(bin_centers(bins),bin_centers(rt_change_bins),rt_change_pdfs'); 
    set(gca,'YDir','normal');
end
    
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

function reactionTimes=getBaseSubtractedRTs(tbt,useAsCue,whichReach)

% cue ind
avcue=nanmean(tbt.(useAsCue),1);
[~,maind]=max(avcue);

% which reaches to get
temp=tbt.(whichReach);
newtemp=nan(size(tbt.(whichReach)));

% get only first reaches
firstreaches=nan(1,size(temp,1));
baselineWindows=zeros(size(temp,1),maind);
% note that everything has already been cue-aligned
for i=1:size(temp,1)
    fi=find(temp(i,maind+1:end)>0,1,'first')+maind;
    % subtract off baseline reaching
    % first, calculate baseline reaching
    fi_base=find(temp(i,maind:-1:1)>0,1,'first');
    if ~isempty(fi_base)
        baselineWindows(i,fi_base)=1;
    end
    if ~isempty(fi)
        fi=fi(1);
        firstreaches(i)=fi;
    end
    if ~isempty(fi)
        newtemp(i,:)=zeros(size(newtemp(i,:)));
        newtemp(i,fi)=1;
    end
end
reactionTimes=(firstreaches-maind).*mode(diff(nanmean(tbt.times,1)));
binToRTmapping=nan(1,size(temp,2));
binToRTmapping(maind+1:end)=[1:length(binToRTmapping(maind+1:end))].*mode(diff(nanmean(tbt.times,1)));
tbt.firstReachesAfterCueOnset=newtemp;
% find time bins with reach # below baseline reach rate
% zero out these reaction times
baseThresh=nanmean(nansum(baselineWindows,1)./size(baselineWindows,1));
RTsummary=smooth(nansum(newtemp,1)./size(newtemp,1),floor(0.3/mode(diff(nanmean(tbt.times,1)))));
whichBinsBelow=RTsummary<=baseThresh;
throwOutAllAfter1stBaseline=0;
if throwOutAllAfter1stBaseline==1
    firstAtBase=find(whichBinsBelow(maind+1:end)==1,1,'first');
    whichBinsBelow(maind+firstAtBase:end)=1;
end
tbt.firstReachesAfterCueOnset(:,whichBinsBelow)=0;
zeroTheseRT=binToRTmapping(whichBinsBelow);
reactionTimes(ismember(reactionTimes,zeroTheseRT))=nan;

end