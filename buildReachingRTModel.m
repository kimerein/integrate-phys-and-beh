function [dataset,correctedDistributions]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir)

% For the following behavioral events, get how delta_RT is changed
% 1 trial forward, 2 trials forward, 3 trials forward, etc.
% Try to get functions for the pdfs(delta_RT) for each of these behavioral
% events, independent of mouse and d-prime
% These pdfs(delta_RT) then serve as model components

% Questions:
% 1. How good is the model?
% How much of the animal's behavior across many trials can be predicted
% from these estimates of the effects of single behavioral events?
% How good is the model across individual mice?
% How good is the model across stages of learning?
%
% 2. What do the forms of the model components tell us about what these
% components represent?

% Behavioral Events:
% Mouse touches pellet
% Mouse consumes pellet
% Mouse reaches but no pellet touch
% Mouse touches pellet and opto stim on
% Mouse consumes pellet and opto stim on
% Mouse reaches but no pellet touch and opto stim on 

% Test Hypothesis 1: 
% If opto silencing of DMS does not change underlying RT distribution but
% only impacts CHANGES in RT distribution, then the sequence trial type A, led
% silencing, and trial type B will look the same as trial type A followed
% immediately by trial type B, in terms of RT change

% Note that, if underlying RT distributions are very different, it will be
% next to impossible to compare delta_RTs as a result of different
% behavioral events
% Thus, need to "correct input distributions" and bootstrap to get matching
% initial RT distributions, in order to focus on delta_RT changes
% For comparisons, resample all RT distributions to match the average RT
% distribution across all trials and all mice



templateSequence1_cond=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
templateSequence1_end=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;


% Get effects of, e.g., "Mouse touches pellet" on RT
nInSequence=[2 3 4];
nRepsForBootstrap=20;
% templateSequence2_cond=any(alltbt.all_reachBatch>0.5,2) & trialTypes.touched_pellet==0 & trialTypes.led==0 & trialTypes.paw_during_wheel==0;
templateSequence2_cond=trialTypes.touched_pellet==1 & trialTypes.led==1 & trialTypes.consumed_pellet==0;
templateSequence2_end=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
% Get raw reaching (all reaches)
% Get reaction times
% Get change in reaction times (non-corrected input distributions)
% Get change in reaction times (corrected input distributions)
dataset.realDistributions.event_name='dprime_1to2point5_mouse_touches_but_does_not_consume_pellet_str_led_on';
dataset.realDistributions.nInSequence=nInSequence;
dataset.realDistributions.nRepsForBootstrap=nRepsForBootstrap;
dataset.realDistributions.templateSequence1_cond=templateSequence1_cond;
dataset.realDistributions.templateSequence1_end=templateSequence1_end;
dataset.realDistributions.templateSequence2_cond=templateSequence2_cond;
dataset.realDistributions.templateSequence2_end=templateSequence2_end;

saveDir=[saveDir '\' dataset.realDistributions.event_name];
mkdir(saveDir);

doCorrectInputDistribution=false;
[dim1_seq_rtchanges_cond1,dim1_seq_rtchanges_cond2,dim2_seq_rtchanges_cond1,dim2_seq_rtchanges_cond2,ns,cond1_isSeq,cond2_isSeq,rr_cond1_trial1,rr_cond1_triali,rr_cond2_trial1,rr_cond2_triali,cond1_rts_trial1,cond1_rts_triali,cond2_rts_trial1,cond2_rts_triali,rr_cond1_trial1_se,rr_cond1_triali_se,rr_cond2_trial1_se,rr_cond2_triali_se,realrtpair_seq1,realrtpair_seq2,alldim_seq_rtchanges_cond1,alldim_seq_rtchanges_cond2]=getFxOfBehaviorEvent(templateSequence1_cond,templateSequence1_end,templateSequence2_cond,templateSequence2_end,nInSequence,1,doCorrectInputDistribution,alltbt,trialTypes,metadata,fakeCueInd);
dataset.realDistributions.(['alldim_rtchanges_allTrialsSequence'])=alldim_seq_rtchanges_cond1;
dataset.realDistributions.(['alldim_rtchanges_event'])=alldim_seq_rtchanges_cond2;
dataset.realDistributions.(['dim1_rtchanges_allTrialsSequence'])=dim1_seq_rtchanges_cond1;
dataset.realDistributions.(['dim1_rtchanges_event'])=dim1_seq_rtchanges_cond2;
dataset.realDistributions.(['dim2_rtchanges_allTrialsSequence'])=dim2_seq_rtchanges_cond1;
dataset.realDistributions.(['dim2_rtchanges_event'])=dim2_seq_rtchanges_cond2;
dataset.realDistributions.(['rawReaching_allTrialsSequence_trial1InSeq'])=rr_cond1_trial1;
dataset.realDistributions.(['rawReaching_allTrialsSequence_trialiInSeq'])=rr_cond1_triali;
dataset.realDistributions.(['rawReaching_event_trial1InSeq'])=rr_cond2_trial1;
dataset.realDistributions.(['rawReaching_event_trialiInSeq'])=rr_cond2_triali;
dataset.realDistributions.(['se_rawReaching_allTrialsSequence_trial1InSeq'])=rr_cond1_trial1_se;
dataset.realDistributions.(['se_rawReaching_allTrialsSequence_trialiInSeq'])=rr_cond1_triali_se;
dataset.realDistributions.(['se_rawReaching_event_trial1InSeq'])=rr_cond2_trial1_se;
dataset.realDistributions.(['se_rawReaching_event_trialiInSeq'])=rr_cond2_triali_se;
dataset.realDistributions.(['allTrialsSequence_isSeq'])=cond1_isSeq;
dataset.realDistributions.(['event_isSeq'])=cond2_isSeq;
dataset.realDistributions.(['allTrialsSequence_RT_trial1InSeq'])=cond1_rts_trial1;
dataset.realDistributions.(['allTrialsSequence_RT_trialiInSeq'])=cond1_rts_triali;
dataset.realDistributions.(['event_RT_trial1InSeq'])=cond2_rts_trial1;
dataset.realDistributions.(['event_RT_trialiInSeq'])=cond2_rts_triali;
dataset.realDistributions.(['ns'])=ns;
dataset.realDistributions.(['realrtpair_seq1'])=realrtpair_seq1;
dataset.realDistributions.(['realrtpair_seq2'])=realrtpair_seq2;

mkdir([saveDir '\real_distributions']); 
save([saveDir '\real_distributions\pdf.mat'],'dataset');

doCorrectInputDistribution=true;
[dim1_seq_rtchanges_cond1,dim1_seq_rtchanges_cond2,dim2_seq_rtchanges_cond1,dim2_seq_rtchanges_cond2,ns,cond1_isSeq,cond2_isSeq,rr_cond1_trial1,rr_cond1_triali,rr_cond2_trial1,rr_cond2_triali,cond1_rts_trial1,cond1_rts_triali,cond2_rts_trial1,cond2_rts_triali,rr_cond1_trial1_se,rr_cond1_triali_se,rr_cond2_trial1_se,rr_cond2_triali_se,realrtpair_seq1,realrtpair_seq2,alldim_seq_rtchanges_cond1,alldim_seq_rtchanges_cond2]=getFxOfBehaviorEvent(templateSequence1_cond,templateSequence1_end,templateSequence2_cond,templateSequence2_end,nInSequence,nRepsForBootstrap,doCorrectInputDistribution,alltbt,trialTypes,metadata,fakeCueInd);
correctedDistributions.event_name=dataset.realDistributions.event_name;
correctedDistributions.nInSequence=nInSequence;
correctedDistributions.nRepsForBootstrap=nRepsForBootstrap;
correctedDistributions.templateSequence1_cond=templateSequence1_cond;
correctedDistributions.templateSequence1_end=templateSequence1_end;
correctedDistributions.templateSequence2_cond=templateSequence2_cond;
correctedDistributions.templateSequence2_end=templateSequence2_end;
correctedDistributions.(['alldim_rtchanges_allTrialsSequence'])=alldim_seq_rtchanges_cond1;
correctedDistributions.(['alldim_rtchanges_event'])=alldim_seq_rtchanges_cond2;
correctedDistributions.(['dim1_rtchanges_allTrialsSequence'])=dim1_seq_rtchanges_cond1;
correctedDistributions.(['dim1_rtchanges_event'])=dim1_seq_rtchanges_cond2;
correctedDistributions.(['dim2_rtchanges_allTrialsSequence'])=dim2_seq_rtchanges_cond1;
correctedDistributions.(['dim2_rtchanges_event'])=dim2_seq_rtchanges_cond2;
correctedDistributions.(['rawReaching_allTrialsSequence_trial1InSeq'])=rr_cond1_trial1;
correctedDistributions.(['rawReaching_allTrialsSequence_trialiInSeq'])=rr_cond1_triali;
correctedDistributions.(['rawReaching_event_trial1InSeq'])=rr_cond2_trial1;
correctedDistributions.(['rawReaching_event_trialiInSeq'])=rr_cond2_triali;
correctedDistributions.(['se_rawReaching_allTrialsSequence_trial1InSeq'])=rr_cond1_trial1_se;
correctedDistributions.(['se_rawReaching_allTrialsSequence_trialiInSeq'])=rr_cond1_triali_se;
correctedDistributions.(['se_rawReaching_event_trial1InSeq'])=rr_cond2_trial1_se;
correctedDistributions.(['se_rawReaching_event_trialiInSeq'])=rr_cond2_triali_se;
correctedDistributions.(['allTrialsSequence_isSeq'])=cond1_isSeq;
correctedDistributions.(['event_isSeq'])=cond2_isSeq;
correctedDistributions.(['allTrialsSequence_RT_trial1InSeq'])=cond1_rts_trial1;
correctedDistributions.(['allTrialsSequence_RT_trialiInSeq'])=cond1_rts_triali;
correctedDistributions.(['event_RT_trial1InSeq'])=cond2_rts_trial1;
correctedDistributions.(['event_RT_trialiInSeq'])=cond2_rts_triali;
correctedDistributions.(['ns'])=ns;
correctedDistributions.(['realrtpair_seq1'])=realrtpair_seq1;
correctedDistributions.(['realrtpair_seq2'])=realrtpair_seq2;

% Save 
f=fieldnames(correctedDistributions);
mkdir([saveDir '\corrected_distributions']); 
for i=1:length(f)
    temp=correctedDistributions.(f{i});
    save([saveDir '\corrected_distributions\' f{i} '.mat'],'temp');
end

end

function [dim1_seq_rtchanges_cond1,dim1_seq_rtchanges_cond2,dim2_seq_rtchanges_cond1,dim2_seq_rtchanges_cond2,ns,cond1_isSeq,cond2_isSeq,rr_cond1_trial1,rr_cond1_triali,rr_cond2_trial1,rr_cond2_triali,cond1_rts_trial1,cond1_rts_triali,cond2_rts_trial1,cond2_rts_triali,rr_cond1_trial1_se,rr_cond1_triali_se,rr_cond2_trial1_se,rr_cond2_triali_se,realrtpair_seq1,realrtpair_seq2,alldim_seq_rtchanges_cond1,alldim_seq_rtchanges_cond2]=getFxOfBehaviorEvent(templateSequence1_cond,templateSequence1_end,templateSequence2_cond,templateSequence2_end,nInSequence,nRepsForBootstrap,doCorrectInputDistribution,alltbt,trialTypes,metadata,fakeCueInd)

allDims_rt_change_cond1=nan(length(nInSequence),nRepsForBootstrap);
allDims_rt_change_cond2=nan(length(nInSequence),nRepsForBootstrap);
allDims_p=nan(length(nInSequence),nRepsForBootstrap);
dim1_rt_change_cond1=nan(length(nInSequence),nRepsForBootstrap);
dim1_rt_change_cond2=nan(length(nInSequence),nRepsForBootstrap);
dim1_p=nan(length(nInSequence),nRepsForBootstrap);
dim2_rt_change_cond1=nan(length(nInSequence),nRepsForBootstrap);
dim2_rt_change_cond2=nan(length(nInSequence),nRepsForBootstrap);
dim2_p=nan(length(nInSequence),nRepsForBootstrap);
dim2_closeup_rt_change_cond1=nan(length(nInSequence),nRepsForBootstrap);
dim2_closeup_rt_change_cond2=nan(length(nInSequence),nRepsForBootstrap);
dim2_closeup_p=nan(length(nInSequence),nRepsForBootstrap);
ns.cond1=nan(length(nInSequence),nRepsForBootstrap);
ns.cond2=nan(length(nInSequence),nRepsForBootstrap);
dim1_seq_rtchanges_cond1=cell(1,length(nInSequence));
dim1_seq_rtchanges_cond2=cell(1,length(nInSequence));
dim2_seq_rtchanges_cond1=cell(1,length(nInSequence));
dim2_seq_rtchanges_cond2=cell(1,length(nInSequence));
cond1_isSeq=cell(1,length(nInSequence));
cond2_isSeq=cell(1,length(nInSequence));
cond1_rts_trial1=cell(1,length(nInSequence));
cond1_rts_triali=cell(1,length(nInSequence));
cond2_rts_trial1=cell(1,length(nInSequence));
cond2_rts_triali=cell(1,length(nInSequence));
rr_cond1_trial1=cell(1,length(nInSequence));
rr_cond1_triali=cell(1,length(nInSequence));
rr_cond2_trial1=cell(1,length(nInSequence));
rr_cond2_triali=cell(1,length(nInSequence));
rr_cond1_trial1_se=cell(1,length(nInSequence));
rr_cond1_triali_se=cell(1,length(nInSequence));
rr_cond2_trial1_se=cell(1,length(nInSequence));
rr_cond2_triali_se=cell(1,length(nInSequence));
realrtpair_seq1=cell(1,length(nInSequence));
realrtpair_seq2=cell(1,length(nInSequence));
for i=1:length(nInSequence)
    nSeq=nInSequence(i);
    disp(nSeq);
    for j=1:nSeq
        if j==1
            templateSequence1{j}=templateSequence1_cond;
            templateSequence2{j}=templateSequence2_cond;
        elseif j==nSeq
            templateSequence1{j}=templateSequence1_end;
            templateSequence2{j}=templateSequence2_end;
        else
            templateSequence1{j}=templateSequence1_cond;
            templateSequence2{j}=templateSequence2_cond;
        end
    end
    allVals1_alldim=[];
    allVals2_alldim=[];
    allVals1_dim1=[];
    allVals2_dim1=[];
    allVals1_dim2=[];
    allVals2_dim2=[];
    rawR_cond1_trial1=[];
    rawR_cond1_triali=[];
    rawR_cond2_trial1=[];
    rawR_cond2_triali=[];
    rawR_cond1_trial1_se=[];
    rawR_cond1_triali_se=[];
    rawR_cond2_trial1_se=[];
    rawR_cond2_triali_se=[];
    allRTs_cond1_trial1=[];
    allRTs_cond1_triali=[];
    allRTs_cond2_trial1=[];
    allRTs_cond2_triali=[];
    realpair_seq1=[];
    realpair_seq2=[];
    for j=1:nRepsForBootstrap
        disp(j);
        [allDims,dim1,dim2,dim2_closeup,rawReaching_cond1_trial1,rawReaching_cond1_triali,rawReaching_cond2_trial1,rawReaching_cond2_triali,isSeq1,isSeq2,rts_seq1_trial1,rts_seq1_triali,rts_seq2_trial1,rts_seq2_triali,rtpairs_real_seq1,rtpairs_real_seq2]=separateCueDep_and_Ind_Components(alltbt,trialTypes,metadata,fakeCueInd,templateSequence1,templateSequence2,nSeq-1,nSeq-1,doCorrectInputDistribution);
        allDims_rt_change_cond1(i,j)=allDims.med1;
        allDims_rt_change_cond2(i,j)=allDims.med2;
        allDims_p(i,j)=allDims.p;
        ns.cond1(i,j)=allDims.n_cond1;
        ns.cond2(i,j)=allDims.n_cond2;
        dim1_rt_change_cond1(i,j)=dim1.med1;
        dim1_rt_change_cond2(i,j)=dim1.med2;
        dim1_p(i,j)=dim1.p;
        dim2_rt_change_cond1(i,j)=dim2.med1;
        dim2_rt_change_cond2(i,j)=dim2.med2;
        dim2_p(i,j)=dim2.p;
        dim2_closeup_rt_change_cond1(i,j)=dim2_closeup.med1;
        dim2_closeup_rt_change_cond2(i,j)=dim2_closeup.med2;
        dim2_closeup_p(i,j)=dim2_closeup.p;
        allVals1_alldim=[allVals1_alldim allDims.vals1];
        allVals2_alldim=[allVals2_alldim allDims.vals2];
        allVals1_dim2=[allVals1_dim2 dim2.vals1];
        allVals2_dim2=[allVals2_dim2 dim2.vals2];
        allVals1_dim1=[allVals1_dim1 dim1.vals1];
        allVals2_dim1=[allVals2_dim1 dim1.vals2];
        % save means of bootstrap
        rawR_cond1_trial1=[rawR_cond1_trial1; nanmean(rawReaching_cond1_trial1,1)];
        rawR_cond1_triali=[rawR_cond1_triali; nanmean(rawReaching_cond1_triali,1)];
        rawR_cond2_trial1=[rawR_cond2_trial1; nanmean(rawReaching_cond2_trial1,1)];
        rawR_cond2_triali=[rawR_cond2_triali; nanmean(rawReaching_cond2_triali,1)];
        % uncertainties add as sqrt of sum of squares
        rawR_cond1_trial1_se=[rawR_cond1_trial1_se; std(rawReaching_cond1_trial1,[],1)./sqrt(size(rawReaching_cond1_trial1,1))];
        rawR_cond1_triali_se=[rawR_cond1_triali_se; std(rawReaching_cond1_triali,[],1)./sqrt(size(rawReaching_cond1_triali,1))];
        rawR_cond2_trial1_se=[rawR_cond2_trial1_se; std(rawReaching_cond2_trial1,[],1)./sqrt(size(rawReaching_cond2_trial1,1))];
        rawR_cond2_triali_se=[rawR_cond2_triali_se; std(rawReaching_cond2_triali,[],1)./sqrt(size(rawReaching_cond2_triali,1))];
        allRTs_cond1_trial1=[allRTs_cond1_trial1 rts_seq1_trial1];
        allRTs_cond1_triali=[allRTs_cond1_triali rts_seq1_triali];
        allRTs_cond2_trial1=[allRTs_cond2_trial1 rts_seq2_trial1];
        allRTs_cond2_triali=[allRTs_cond2_triali rts_seq2_triali];
        realpair_seq1=[realpair_seq1 rtpairs_real_seq1];
        realpair_seq2=[realpair_seq2 rtpairs_real_seq2];
    end
    alldim_seq_rtchanges_cond1{i}=allVals1_alldim;
    alldim_seq_rtchanges_cond2{i}=allVals2_alldim;
    dim1_seq_rtchanges_cond1{i}=allVals1_dim1;
    dim1_seq_rtchanges_cond2{i}=allVals2_dim1;
    dim2_seq_rtchanges_cond1{i}=allVals1_dim2;
    dim2_seq_rtchanges_cond2{i}=allVals2_dim2;
    cond1_isSeq{i}=isSeq1;
    cond2_isSeq{i}=isSeq2;
    rr_cond1_trial1{i}=rawR_cond1_trial1;
    rr_cond1_triali{i}=rawR_cond1_triali;
    rr_cond2_trial1{i}=rawR_cond2_trial1;
    rr_cond2_triali{i}=rawR_cond2_triali;
    rr_cond1_trial1_se{i}=rawR_cond1_trial1_se;
    rr_cond1_triali_se{i}=rawR_cond1_triali_se;
    rr_cond2_trial1_se{i}=rawR_cond2_trial1_se;
    rr_cond2_triali_se{i}=rawR_cond2_triali_se;
    cond1_rts_trial1{i}=allRTs_cond1_trial1;
    cond1_rts_triali{i}=allRTs_cond1_triali;
    cond2_rts_trial1{i}=allRTs_cond2_trial1;
    cond2_rts_triali{i}=allRTs_cond2_triali;
    realrtpair_seq1{i}=realpair_seq1;
    realrtpair_seq2{i}=realpair_seq2;
end

% i=1;
% temp1=forHists_cond1{i};
% temp2=forHists_cond2{i};
% plotHist(temp1,temp2,200,'Histo','y-vals');
% plotCDF(temp1,temp2,200,'CDF');
% range_learning=2; % in seconds
% plotHist(temp1(temp1>-range_learning & temp1<range_learning),temp2(temp2>-range_learning & temp2<range_learning),200,'Close-up Histo','y-vals');
% plotCDF(temp1(temp1>-range_learning & temp1<range_learning),temp2(temp2>-range_learning & temp2<range_learning),200,'Close-up CDF');
% disp('CLOSE-UP P-VAL');
% testRanksum(temp1(temp1>-range_learning & temp1<range_learning),temp2(temp2>-range_learning & temp2<range_learning),1);
% 
% 
% for i=1:length(nInSequence)
%     disp([num2str(nanmean(ns.cond1(i,:),2)) ' trials for cond1 of nInSequence ' num2str(nInSequence(i))]);
%     disp([num2str(nanmean(ns.cond2(i,:),2)) ' trials for cond2 of nInSequence ' num2str(nInSequence(i))]);
% end

end

function plotOutput(nInSequence,allDims_rt_change_cond1,allDims_rt_change_cond2,allDims_p,tit,nRepsForBootstrap)

figure();
plot(nInSequence,nanmean(allDims_rt_change_cond1-allDims_rt_change_cond2,2),'Color','k');
hold on;
plot(nInSequence,nanmean(allDims_p,2),'Color','b');
legend({'Median shift in RT change','p-value'});
for i=1:length(nInSequence)
    temp1=allDims_rt_change_cond1(i,:);
    temp2=allDims_rt_change_cond2(i,:);
    scatter(ones(size(temp1)).*nInSequence(i),temp1,[],'k');
    scatter(ones(size(temp2)).*nInSequence(i),temp2,[],'r');
    line([nInSequence(i) nInSequence(i)],[nanmean(temp1-temp2,2)-nanstd(temp1-temp2,[],2)./sqrt(nRepsForBootstrap) nanmean(temp1-temp2,2)+(nanstd(temp1-temp2,[],2)./sqrt(nRepsForBootstrap))],'Color','k');
end
for i=1:length(nInSequence)
    temp1=allDims_p(i,:);
    scatter(ones(size(temp1)).*nInSequence(i),temp1,[],'b');
    line([nInSequence(i) nInSequence(i)],[nanmean(temp1,2)-(nanstd(temp1,[],2)./sqrt(nRepsForBootstrap)) nanmean(temp1,2)+(nanstd(temp1,[],2)./sqrt(nRepsForBootstrap))],'Color','b');
end

title(tit);

end

function [allDims,dim1,dim2,dim2_closeup,rawReaching_cond1_trial1,rawReaching_cond1_triali,rawReaching_cond2_trial1,rawReaching_cond2_triali,isSeq1,isSeq2,rts_seq1_trial1,rts_seq1_triali,rts_seq2_trial1,rts_seq2_triali,rtpairs_real_seq1,rtpairs_real_seq2]=separateCueDep_and_Ind_Components(alltbt,trialTypes,metadata,fakeCueInd,templateSequence1,templateSequence2,nNext1,nNext2,doCorrectInputDistribution)

% separate cue-dependent and cue-independent components of RT change

% plotting settings
jitterStep=0.01;
bins=200;

% Get RT change after real cue

% Format 3 (10 args); tbt, trialTypes, metadata, templateSequence1,
%   nNext1, templateSequence2, nNext2, useAsCue, whichReach, doPlot
% disp('Plotting results for real post-cue reaction times');
[~,~,rt_pairs]=compareTrialCombos(alltbt,trialTypes,metadata,templateSequence1,nNext1,templateSequence2,nNext2,'cueZone_onVoff','all_reachBatch',0);
% show where these template sequences come from as a function of fraction
% through session
% [n,x]=hist(metadata.nth_session(rt_pairs.sequenceMatchStarts1==1)+metadata.fractionThroughSess(rt_pairs.sequenceMatchStarts1==1),bins);
% figure();
% plot(x,n,'Color','k');
% hold on;
% [n,x]=hist(metadata.nth_session(rt_pairs.sequenceMatchStarts2==1)+metadata.fractionThroughSess(rt_pairs.sequenceMatchStarts2==1),bins);
% plot(x,n,'Color','r');
% title('Template match from nth session');

temp=find(rt_pairs.sequenceMatchStarts1==1);
rawReaching_cond1_trial1=alltbt.all_reachBatch(temp,:);
rawReaching_cond1_triali=alltbt.all_reachBatch(temp+nNext1,:);
temp=find(rt_pairs.sequenceMatchStarts2==1);
rawReaching_cond2_trial1=alltbt.all_reachBatch(temp,:);
rawReaching_cond2_triali=alltbt.all_reachBatch(temp+nNext2,:);

% Match reaction time distributions for conditions 1 and 2
if doCorrectInputDistribution==true
    [rt_pairs,resampInds]=correctInputDistribution(rt_pairs,1,[]);
end
isSeq1=rt_pairs.sequenceMatchStarts1;
isSeq2=rt_pairs.sequenceMatchStarts2;
rtpairs_real_seq1=rt_pairs.real_rt_pair1(rt_pairs.sequenceMatchStarts1==1);
rtpairs_real_seq2=rt_pairs.real_rt_pair2(rt_pairs.sequenceMatchStarts2==1);
rts_seq1_trial1=rt_pairs.all_rt1(rt_pairs.sequenceMatchStarts1==1);
rts_seq1_triali=rt_pairs.all_rt1_triali(rt_pairs.sequenceMatchStarts1==1);
rts_seq2_trial1=rt_pairs.all_rt2(rt_pairs.sequenceMatchStarts2==1);
rts_seq2_triali=rt_pairs.all_rt2_triali(rt_pairs.sequenceMatchStarts2==1);

% rt_pairs.rt_pairs1_contingent contains reaction time pairs sub-sampled to
% include only real pairs (not separated by a file split) and belonging to
% templateSequence1
real_rt_pairs1=rt_pairs.rt_pairs1_contingent;
real_rt_pairs2=rt_pairs.rt_pairs2_contingent;

% plotHist(rt_pairs.rt_pairs1_contingent,rt_pairs.rt_pairs2_contingent,bins,'RT change','Count');
% plotCDF(rt_pairs.rt_pairs1_contingent,rt_pairs.rt_pairs2_contingent,bins,'CDF');
[p_allDims,med1_allDims,med2_allDims]=testRanksum(rt_pairs.rt_pairs1_contingent,rt_pairs.rt_pairs2_contingent,0);
allDims.p=p_allDims;
allDims.med1=med1_allDims;
allDims.med2=med2_allDims;
allDims.vals1=real_rt_pairs1;
allDims.vals2=real_rt_pairs2;
allDims.n_cond1=nansum(rt_pairs.sequenceMatchStarts1==1);
allDims.n_cond2=nansum(rt_pairs.sequenceMatchStarts2==1);

% Get RT change after fake cue

% Make fake cue for measuring cue-independent shift in reaction time
% figure(); plot(nanmean(alltbt.all_reachBatch,1),'Color','k');
% if isempty(fakeCueInd)
%     fakeCueInd=input('Enter index of first pre-cue reaching across all trials. Will put fake cue here.');
% end
% hold all;
% line([fakeCueInd fakeCueInd],[0 nanmax(nanmean(alltbt.all_reachBatch,1))]);

alltbt.fakeCue_for_preCueRT=zeros(size(alltbt.cue));
alltbt.fakeCue_for_preCueRT(:,fakeCueInd)=1;

% disp('Plotting results for fake pre-cue reaction times');
[~,~,rt_pairs]=compareTrialCombos(alltbt,trialTypes,metadata,templateSequence1,nNext1,templateSequence2,nNext2,'fakeCue_for_preCueRT','all_reachBatch',0);

% Match reaction time distributions for conditions 1 and 2
if doCorrectInputDistribution==true
    [rt_pairs]=correctInputDistribution(rt_pairs,1,resampInds);
end

fake_rt_pairs1=rt_pairs.rt_pairs1_contingent;
fake_rt_pairs2=rt_pairs.rt_pairs2_contingent;

% Calculate PCAed matrix that orthogonalizes cue-dependent and cue-independent components

% Plot RT change before PCA
% jitter1=rand(size(real_rt_pairs1)).*jitterStep;
% jitter2=rand(size(real_rt_pairs2)).*jitterStep;
% plotScatter(real_rt_pairs1,fake_rt_pairs1,real_rt_pairs2,fake_rt_pairs2,jitter1,jitter2,'RT change before rotation');

coeff=[0.6110    0.7885   -0.0709; ...
       0.5554   -0.4908   -0.6713; ...
       0.5641   -0.3707    0.7378];
% disp(coeff);
pc1=coeff(:,1);
pc2=coeff(:,2);
pc3=coeff(:,3);

% hold on; line([0 pc1(1)],[0 pc1(2)],'Color','b'); 
% hold on; line([0 pc2(1)],[0 pc2(2)],'Color','g');
% hold on; line([0 pc3(1)],[0 pc3(2)],'Color','c');

% Project onto each PC
Z_input=[real_rt_pairs1 real_rt_pairs2; fake_rt_pairs1 fake_rt_pairs2; fake_rt_pairs1 fake_rt_pairs2];
Z_onto_pc1=nan(1,size(Z_input,2));
Z_onto_pc2=nan(1,size(Z_input,2));
Z_onto_pc3=nan(1,size(Z_input,2));
for i=1:size(Z_input,2)
    currData=Z_input(:,i);
    % Z_onto_pc1(:,i)=bsxfun(@times, dot(currData,pc1), pc1);
    Z_onto_pc1(i)=dot(currData,pc1);
    Z_onto_pc2(i)=dot(currData,pc2);
    Z_onto_pc3(i)=dot(currData,pc3);    
end

% Plot RT change after rotation
% plotScatter(Z_onto_pc1(1:length(real_rt_pairs1)),Z_onto_pc2(1:length(real_rt_pairs1)),Z_onto_pc1(length(real_rt_pairs1)+1:end),Z_onto_pc2(length(real_rt_pairs1)+1:end),jitter1,jitter2,'RT change after rotation');

Zpca=[Z_onto_pc1; Z_onto_pc2];
dim=1;
% plotHist(Zpca(dim,1:length(real_rt_pairs1)),Zpca(dim,length(real_rt_pairs1)+1:end),bins,'Histo of dim 1','x-vals');
% plotCDF(Zpca(dim,1:length(real_rt_pairs1)),Zpca(dim,length(real_rt_pairs1)+1:end),bins,'CDF of dim 1');
[p,med1,med2]=testRanksum(Zpca(dim,1:length(real_rt_pairs1)),Zpca(dim,length(real_rt_pairs1)+1:end),0);
dim1.p=p;
dim1.med1=med1;
dim1.med2=med2;
temp1=Zpca(1,1:length(real_rt_pairs1));
temp2=Zpca(1,length(real_rt_pairs1)+1:end);
dim1.vals1=temp1;
dim1.vals2=temp2;

dim=2;
% plotHist(Zpca(dim,1:length(real_rt_pairs1)),Zpca(dim,length(real_rt_pairs1)+1:end),bins,'Histo of dim 2','y-vals');
% plotCDF(Zpca(dim,1:length(real_rt_pairs1)),Zpca(dim,length(real_rt_pairs1)+1:end),bins,'CDF of dim 2');
[p,med1,med2]=testRanksum(Zpca(dim,1:length(real_rt_pairs1)),Zpca(dim,length(real_rt_pairs1)+1:end),0);
dim2.p=p;
dim2.med1=med1;
dim2.med2=med2;
temp1=Zpca(dim,1:length(real_rt_pairs1));
temp2=Zpca(dim,length(real_rt_pairs1)+1:end);
dim2.vals1=temp1;
dim2.vals2=temp2;

temp1=Zpca(dim,1:length(real_rt_pairs1));
temp2=Zpca(dim,length(real_rt_pairs1)+1:end);
% plotCDF(temp1(temp1>-2 & temp1<2),temp2(temp2>-2 & temp2<2),bins,'CDF of dim 2 -- only changes within 2 sec faster or slower');
% disp('CLOSE-UP P-VAL');
[p,med1,med2]=testRanksum(temp1(temp1>-4 & temp1<4),temp2(temp2>-4 & temp2<4),0);
dim2_closeup.p=p;
dim2_closeup.med1=med1;
dim2_closeup.med2=med2;

end

function plotScatter(RT_pairs1_x,RT_pairs1_y,RT_pairs2_x,RT_pairs2_y,jitter1,jitter2,tit)

figure();
if length(jitter1)==1
    useX=RT_pairs1_x+rand(size(RT_pairs1_x)).*jitter;
    useY=RT_pairs1_y+rand(size(RT_pairs1_y)).*jitter;
else
    useX=RT_pairs1_x+jitter1;
    useY=RT_pairs1_y+jitter1;
end
    
s=scatter(useX,useY,100,'k','filled');
set(gcf,'position',[10,10,1000,1000]);
pause;
m=get(s,'MarkerHandle');
alpha=0.3;
m.FaceColorData=uint8(255*[0;0;0;alpha]);
xlabel('real cue');
ylabel('fake cue');
title(tit);
% xlim([0 9.5]);
% ylim([0 9.5]);

hold on;
if length(jitter2)==1
    useX=RT_pairs2_x+rand(size(RT_pairs2_x)).*jitter;
    useY=RT_pairs2_y+rand(size(RT_pairs2_y)).*jitter;
else
    useX=RT_pairs2_x+jitter2;
    useY=RT_pairs2_y+jitter2;
end
s=scatter(useX,useY,100,'r','filled');
pause;
m=get(s,'MarkerHandle');
alpha=0.3;
m.FaceColorData=uint8(255*[1;0;0;alpha]);

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

function plotCDF(data1,data2,bins,tit)

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

function plotHist(data1,data2,bins,tit,xlab)

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