function [forHists_cond1,forHists_cond2]=showHowLearningIncrements(alltbt,trialTypes,metadata,fakeCueInd)

templateSequence1_cond=trialTypes.touched_pellet==1 & trialTypes.led==0;
templateSequence1_end=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;

templateSequence2_cond=trialTypes.touched_pellet==0 & trialTypes.led==0;
templateSequence2_end=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;

% templateSequence1_cond=trialTypes.touched_pellet==1 & trialTypes.led==0;
% templateSequence1_end=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
% 
% templateSequence2_cond=trialTypes.touched_pellet==1 & trialTypes.led==1;
% templateSequence2_end=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;

% nInSequence=[2 3 4];
nInSequence=[2 4];

nRepsForBootstrap=20;

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
forHists_cond1=cell(1,length(nInSequence));
forHists_cond2=cell(1,length(nInSequence));
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
    allVals1=[];
    allVals2=[];
    for j=1:nRepsForBootstrap
        disp(j);
        [allDims,dim1,dim2,dim2_closeup]=separateCueDep_and_Ind_Components(alltbt,trialTypes,metadata,fakeCueInd,templateSequence1,templateSequence2,nSeq-1,nSeq-1);
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
        allVals1=[allVals1 dim2_closeup.vals1];
        allVals2=[allVals2 dim2_closeup.vals2];
%         allVals1=[allVals1 dim1.vals1];
%         allVals2=[allVals2 dim1.vals2];
    end
    forHists_cond1{i}=allVals1;
    forHists_cond2{i}=allVals2;
end

plotOutput(nInSequence,allDims_rt_change_cond1,allDims_rt_change_cond2,allDims_p,'All dimensions',nRepsForBootstrap);
plotOutput(nInSequence,dim1_rt_change_cond1,dim1_rt_change_cond2,dim1_p,'Dimension 1 -- cue-independent',nRepsForBootstrap);
plotOutput(nInSequence,dim2_rt_change_cond1,dim2_rt_change_cond2,dim2_p,'Dimension 2 -- cue-dependent',nRepsForBootstrap);
plotOutput(nInSequence,dim2_closeup_rt_change_cond1,dim2_closeup_rt_change_cond2,dim2_closeup_p,'Dimension 2 -- cue-dependent (close-up)',nRepsForBootstrap);

i=1;
temp1=forHists_cond1{i};
temp2=forHists_cond2{i};
plotHist(temp1,temp2,200,'Histo','y-vals');
plotCDF(temp1,temp2,200,'CDF');
range_learning=2; % in seconds
plotHist(temp1(temp1>-range_learning & temp1<range_learning),temp2(temp2>-range_learning & temp2<range_learning),200,'Close-up Histo','y-vals');
plotCDF(temp1(temp1>-range_learning & temp1<range_learning),temp2(temp2>-range_learning & temp2<range_learning),200,'Close-up CDF');
disp('CLOSE-UP P-VAL');
testRanksum(temp1(temp1>-range_learning & temp1<range_learning),temp2(temp2>-range_learning & temp2<range_learning),1);


for i=1:length(nInSequence)
    disp([num2str(nanmean(ns.cond1(i,:),2)) ' trials for cond1 of nInSequence ' num2str(nInSequence(i))]);
    disp([num2str(nanmean(ns.cond2(i,:),2)) ' trials for cond2 of nInSequence ' num2str(nInSequence(i))]);
end

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

function [allDims,dim1,dim2,dim2_closeup,ns]=separateCueDep_and_Ind_Components(alltbt,trialTypes,metadata,fakeCueInd,templateSequence1,templateSequence2,nNext1,nNext2)

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

% Match reaction time distributions for conditions 1 and 2
[rt_pairs,resampInds]=correctInputDistribution(rt_pairs,1,[]);

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
[rt_pairs]=correctInputDistribution(rt_pairs,1,resampInds);

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
% plotCDF(temp1(temp1>-2 & temp1<2),temp2(temp2>-2 & temp2<2),bins,'CDF of dim 2 -- only changes within 2 sec faster or slower');
% disp('CLOSE-UP P-VAL');
[p,med1,med2]=testRanksum(temp1(temp1>-4 & temp1<4),temp2(temp2>-4 & temp2<4),0);
dim2_closeup.p=p;
dim2_closeup.med1=med1;
dim2_closeup.med2=med2;
dim2_closeup.vals1=temp1;
dim2_closeup.vals2=temp2;

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