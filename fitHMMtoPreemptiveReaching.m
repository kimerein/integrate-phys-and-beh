function [rt_pdf_outs,best_p_cc,best_p_pp,best_switchingTime,ss,preemptive_process,cued_process]=fitHMMtoPreemptiveReaching(x,y,differ,dataset,alltbt)

bins=x;
bin_centers = @(bins) nanmean([bins(1:end-1); bins(2:end)],1);

% preemptive_process.shape=1;
% preemptive_process.rate=0.1;
% preemptive_process.shape=1.2; % best ss assumption
% preemptive_process.rate=0.1; % best ss assumption
% preemptive_process.shape=1.5; % best ss assumption
% preemptive_process.rate=0.11; % best ss assumption
preemptive_process.shape=1.5; % best ss assumption
preemptive_process.rate=0.15; % best ss assumption
prepdf=gampdf(x,preemptive_process.shape,preemptive_process.rate);
preemptive_process.preempt_pdf=prepdf;
preemptive_process.preempt_pdf=preemptive_process.preempt_pdf./nansum(preemptive_process.preempt_pdf);

bins=0:0.035:9;
% Use baseline-subtracted reaching for cue pdf
rts=getBaseSubtractedRTs(alltbt,'cueZone_onVoff','all_reachBatch');
% rts=dataset.realDistributions.allTrialsSequence_RT_trial1InSeq{1};
n=histcounts(rts,0:0.035:9);
cued_process.rt_pdf=n;
cued_process.rt_pdf=cued_process.rt_pdf./nansum(cued_process.rt_pdf);

% fitMask=repmat(y>-0.3 & y<0.3,length(bin_centers(bins)),1);
% fitMask(bin_centers(bins)>0.45,:)=0;
fitMask=repmat(y>-0.4 & y<0.4,length(bin_centers(bins)),1);
fitMask(bin_centers(bins)>0.5,:)=0;
% fitMask=ones(length(x),length(y));

% try_p_cc=[0.1:0.1:1-0.1];
% try_p_pp=try_p_cc+0.1;
% try_switchingTimeCost=[0:0.25:0.8];

% try_p_cc=[0.15:0.01:0.25];
% try_p_pp=[0.35:0.01:0.45];
% try_switchingTimeCost=[0:0.025:0.4];

% try_p_cc=[0.18:0.01:0.2];
% try_p_pp=[0.43:0.01:0.46];
% try_switchingTimeCost=[0.18:0.01:0.22];

% try_p_cc=[0.17]; % best
% try_p_pp=[0.45]; % best
% try_switchingTimeCost=[0.18]; % best

% try_p_cc=[0.15:0.01:0.19]; % best
% try_p_pp=[0.43:0.01:0.47]; % best
% try_switchingTimeCost=[0.16:0.01:0.2]; % best

% try_p_cc=[0.2]; % best ss assumption
% try_p_pp=[0.39]; % best ss assumption
% try_switchingTimeCost=[0.18]; % best ss assumption

% try_p_cc=[0.2]; % best ss assumption
% try_p_pp=[0.39]; % best ss assumption
% try_switchingTimeCost=[0.18]; % best ss assumption

% try_p_cc=[0.45:0.01:0.55]; % best ss assumption
% try_p_pp=[0.25:0.01:0.35]; % best ss assumption
% try_switchingTimeCost=[0.05:0.01:0.15]; % best ss assumption

try_p_cc=[0.5]; % best ss assumption
try_p_pp=[0.31]; % best ss assumption
try_switchingTimeCost=[0.09]; % best ss assumption

differ=differ-nanmean(nanmean(differ(repmat(y<-6,length(bin_centers(bins)),1))));
differ=differ./nansum(nansum(abs(differ(fitMask==1))));

subtractCued=1;
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
cued_process.rt_change_pdfs=cued_process.rt_change_pdfs./nansum(nansum(abs(cued_process.rt_change_pdfs(fitMask==1))));

% Search space for best fit
disp(['Counting to ' num2str(length(try_p_cc))]);
cost=nan(length(try_p_cc),length(try_p_pp),length(try_switchingTimeCost));
for i=1:length(try_p_cc)
    disp(i);
    for j=1:length(try_p_pp)
        disp(j);
        for k=1:length(try_switchingTimeCost)
            [rt_pdf_outs,ss]=makeHMMforPreemptiveReaching(preemptive_process,cued_process,bins,try_p_cc(i),try_p_pp(j),try_switchingTimeCost(k));
            rt_change_pdfs=getRTchangeFromRTpdfs(bins,rt_pdf_outs,rts,preemptive_process,ss.cued,ss.preempt);
            rt_change_pdfs=rt_change_pdfs-nanmean(nanmean(rt_change_pdfs(repmat(y<-6,length(bin_centers(bins)),1))));
            rt_change_pdfs=rt_change_pdfs./nansum(nansum(abs(rt_change_pdfs(fitMask==1))));
            preemptiveOnly=rt_change_pdfs;
            if subtractCued==1
                rt_change_pdfs=rt_change_pdfs-cued_process.rt_change_pdfs;
                rt_change_pdfs=rt_change_pdfs-nanmean(nanmean(rt_change_pdfs(repmat(y<-6,length(bin_centers(bins)),1))));
                rt_change_pdfs=rt_change_pdfs./nansum(nansum(abs(rt_change_pdfs(fitMask==1))));
            end
            % Get difference with data
            temp=abs(differ-rt_change_pdfs);
            cost(i,j,k)=nansum(nansum(temp(fitMask==1)));
%             figure(); 
%             imagesc(x,y,[differ rt_change_pdfs]'); 
%             set(gca,'YDir','normal');
%             figure(); 
%             imagesc(x,y,(differ-rt_change_pdfs)'); 
%             set(gca,'YDir','normal');
        end
    end
end
[fitqualval,fitind]=nanmin(cost(1:end)); 
disp('fit cost');
disp(fitqualval);
[a_i,b_i,c_i]=ind2sub(size(cost),fitind);
best_p_cc=try_p_cc(a_i);
best_p_pp=try_p_pp(b_i);
best_switchingTime=try_switchingTimeCost(c_i);

[rt_pdf_outs,ss]=makeHMMforPreemptiveReaching(preemptive_process,cued_process,bins,best_p_cc,best_p_pp,best_switchingTime);
rt_change_pdfs=getRTchangeFromRTpdfs(bins,rt_pdf_outs,rts,preemptive_process,ss.cued,ss.preempt);
rt_change_pdfs=rt_change_pdfs-nanmean(nanmean(rt_change_pdfs(repmat(y<-6,length(bin_centers(bins)),1))));
rt_change_pdfs=rt_change_pdfs./nansum(nansum(abs(rt_change_pdfs(fitMask==1))));
preemptiveOnly=rt_change_pdfs;
if subtractCued==1
    rt_change_pdfs=rt_change_pdfs-cued_process.rt_change_pdfs;
    rt_change_pdfs=rt_change_pdfs-nanmean(nanmean(rt_change_pdfs(repmat(y<-6,length(bin_centers(bins)),1))));
    rt_change_pdfs=rt_change_pdfs./nansum(nansum(abs(rt_change_pdfs(fitMask==1))));
end

figure(); imagesc(x,x,rt_pdf_outs); title('RT pdf outs');
figure(); imagesc(x,y,conv2(rt_change_pdfs,ones(7),'same')'); set(gca,'YDir','normal');
figure(); imagesc(x,y,differ'); set(gca,'YDir','normal');
figure(); imagesc(x,y,[differ rt_change_pdfs]'); set(gca,'YDir','normal');

end

function reactionTimes=getBaseSubtractedRTs(tbt,useAsCue,whichReach)

%function [reactionTimes,tbt,reactionTimes_preCue]=plotOnlyFirstReach(tbt,ds,whichReach,useAsCue,out,outfield,doPlot)

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

function rt_change_pdfs=getRTchangeFromRTpdfs(bins,rt_pdf_outs,rts,preemptive_process,cue_coeff,preempt_coeff)

% initial_input_rt_pdf='uniform';
initial_input_rt_pdf='steady state';

rt_pdf_outs=rt_pdf_outs./nansum(nansum(rt_pdf_outs));

rt_pdf = @(rts,bins) histcounts(rts,bins)./nansum(histcounts(rts,bins));
bin_centers = @(bins) nanmean([bins(1:end-1); bins(2:end)],1);
try_curr_rts=bin_centers(bins);
rt_change_pdfs=nan(length(try_curr_rts),length(bin_centers(unique([-fliplr(bins) bins]))));

switch initial_input_rt_pdf
    case 'uniform'
        init_rt_pdf=ones(size(bin_centers(bins)));
        init_rt_pdf=init_rt_pdf./nansum(init_rt_pdf);
        for i=1:length(try_curr_rts)
            rt_pdf_out=rt_pdf_outs(i,:);
            if all(rt_pdf_out==0)
                rt_pdf_out(1:end)=0.000000001;
            end
            [rt_change_pdfs(i,:),rt_change_bins]=getRTchange_fromCurrAndUpdatedRTpdfs(init_rt_pdf,rt_pdf_out,bins,true,[bins(i) bins(i+1)]); % if last argument is empty, will plot change distribution for all current reaction times
        end
    case 'steady state'
        init_rt_pdf=rt_pdf(rts,bins);
        init_rt_pdf=init_rt_pdf./nansum(init_rt_pdf);
        for i=1:length(try_curr_rts)
            rt_pdf_out=rt_pdf_outs(i,:);
            if all(rt_pdf_out==0)
                rt_pdf_out(1:end)=0.000000001;
            end
            [rt_change_pdfs1(i,:),rt_change_bins]=getRTchange_fromCurrAndUpdatedRTpdfs(init_rt_pdf,rt_pdf_out,bins,true,[bins(i) bins(i+1)]); % if last argument is empty, will plot change distribution for all current reaction times
        end
        init_rt_pdf=gampdf(bin_centers(bins),preemptive_process.shape,preemptive_process.rate); 
        init_rt_pdf=init_rt_pdf./nansum(init_rt_pdf);
        for i=1:length(try_curr_rts)
            rt_pdf_out=rt_pdf_outs(i,:);
            if all(rt_pdf_out==0)
                rt_pdf_out(1:end)=0.000000001;
            end
            [rt_change_pdfs2(i,:),rt_change_bins]=getRTchange_fromCurrAndUpdatedRTpdfs(init_rt_pdf,rt_pdf_out,bins,true,[bins(i) bins(i+1)]); % if last argument is empty, will plot change distribution for all current reaction times
        end
        rt_change_pdfs=cue_coeff*rt_change_pdfs1+preempt_coeff*rt_change_pdfs2;
end

end

function [rt_change_pdf,rt_change_bins]=getRTchange_fromCurrAndUpdatedRTpdfs(curr_rt_pdf,updated_rt_pdf,bins,suppressOutput,curr_rt)

% Get the distribution of reaction time changes expected given the input
% reaction time distribution (curr_rt_pdf) and the reaction time
% distribution on the next trial (updated_rt_pdf)

% This is the difference of these two random variables
sample_n=30000;
% sample_n=5000;
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

function [rt_pdf_outs,ss]=makeHMMforPreemptiveReaching(preemptive_process,cued_process,rt_bins,p_cc,p_pp,switchingTimeCost)

p_cp=1-p_cc;
p_pc=1-p_pp;

% Transition probabilities from each state must sum to 1
temp1=p_pp/(p_pp+p_pc);
temp2=p_pc/(p_pp+p_pc);
p_pp=temp1;
p_pc=temp2;

temp1=p_cc/(p_cc+p_cp);
temp2=p_cp/(p_cc+p_cp);
p_cc=temp1;
p_cp=temp2;

bin_centers = @(bins) nanmean([bins(1:end-1); bins(2:end)],1);

% PDF of reaction times for preemptive reaching
preempt_pdf=preemptive_process.preempt_pdf;

% PDF of reaction times for cued reaching
cued_pdf=cued_process.rt_pdf;
cued_pdf=cued_pdf./nansum(cued_pdf);
cued_pdf(cued_pdf<=0)=0.00000001;

% Simulate
nSamples=30000;
% nSamples=5000;
% Choose initial state randomly
state=rand<0.5; % if state is 0, cued reaching; if state is 1, preemptive reaching
all_rts=nan(1,nSamples);
all_pdfs=zeros(length(bin_centers(rt_bins)),length(bin_centers(rt_bins)));
justSwitched=1;
timeStep=bin_centers(rt_bins);
timeStep=timeStep(2)-timeStep(1);
whichState=nan(1,nSamples);
for i=1:nSamples
    whichState(i)=state;
    if state==0
        % cued reaching
        % draw curr_rt from cued rt pdf
        if justSwitched==1
            all_rts(i)=randsample(switchingTimeCost+bin_centers(rt_bins),1,true,cued_pdf);
            if all_rts(i)>bin_centers(rt_bins)
                f=length(bin_centers(rt_bins));
            else
                f=find(bin_centers(rt_bins)>=all_rts(i),1,'first');
            end
            all_pdfs(f,:)=all_pdfs(f,:)+[zeros(1,floor(switchingTimeCost/timeStep)) cued_pdf(1:end-floor(switchingTimeCost/timeStep))];
        else
            all_rts(i)=randsample(bin_centers(rt_bins),1,true,cued_pdf);
            %         all_pdfs(find(bin_centers(rt_bins)>=all_rts(i),1,'first'),:)=all_pdfs(find(bin_centers(rt_bins)>=all_rts(i),1,'first'),:)-cued_pdf;
            all_pdfs(find(bin_centers(rt_bins)>=all_rts(i),1,'first'),:)=all_pdfs(find(bin_centers(rt_bins)>=all_rts(i),1,'first'),:)+cued_pdf;
        end
        % test for transition
        % also transition in proportion to how certain I am that in this
        % state
        if rand<p_cc.*cued_pdf(find(bin_centers(rt_bins)>=all_rts(i),1,'first'))/(p_cc.*cued_pdf(find(bin_centers(rt_bins)>=all_rts(i),1,'first'))+p_cp.*preempt_pdf(find(bin_centers(rt_bins)>=all_rts(i),1,'first')))
            % stay in this state
            justSwitched=0;
        else
            state=1;
            justSwitched=1;
        end
    elseif state==1
        % preemptive reaching
        % draw curr_rt from preemptive rt pdf
        if justSwitched==1
            all_rts(i)=randsample(switchingTimeCost+bin_centers(rt_bins),1,true,preempt_pdf);
            if all_rts(i)>bin_centers(rt_bins)
                f=length(bin_centers(rt_bins));
            else
                f=find(bin_centers(rt_bins)>=all_rts(i),1,'first');
            end
            all_pdfs(f,:)=all_pdfs(f,:)+[zeros(1,floor(switchingTimeCost/timeStep)) preempt_pdf(1:end-floor(switchingTimeCost/timeStep))];
        else
            all_rts(i)=randsample(bin_centers(rt_bins),1,true,preempt_pdf);
            all_pdfs(find(bin_centers(rt_bins)>=all_rts(i),1,'first'),:)=all_pdfs(find(bin_centers(rt_bins)>=all_rts(i),1,'first'),:)+preempt_pdf;
        end
        % test for transition
        % also transition in proportion to how certain I am that in this
        % state
        if rand<p_pp.*preempt_pdf(find(bin_centers(rt_bins)>=all_rts(i),1,'first'))/(p_pp.*preempt_pdf(find(bin_centers(rt_bins)>=all_rts(i),1,'first'))+p_pc*cued_pdf(find(bin_centers(rt_bins)>=all_rts(i),1,'first')))
            % stay in this state
            justSwitched=0;
        else
            state=0;
            justSwitched=1;
        end
    end
end
ss.cued=nansum(whichState==0)./nSamples;
ss.preempt=nansum(whichState==1)./nSamples;

% Make plot of output rt pdfs as a function of previous trial rt
rt_pdf_outs=nan(length(bin_centers(rt_bins)),length(bin_centers(rt_bins)));
for i=1:length(bin_centers(rt_bins))
    curr_rt_bin=[rt_bins(i) rt_bins(i+1)];
    f=find(all_rts>=curr_rt_bin(1) & all_rts<curr_rt_bin(2));
    f=f+1;
    f=f(f<=length(all_rts));
    next_rts=all_rts(f);
    [n,x]=histcounts(next_rts,rt_bins);
    rt_pdf_outs(i,:)=n;
end
% figure();
% imagesc(bin_centers(rt_bins),bin_centers(rt_bins),rt_pdf_outs);
% 
% figure();
% imagesc(bin_centers(rt_bins),bin_centers(rt_bins),all_pdfs);

end