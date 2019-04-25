function separateCueDep_and_Ind_Components_of_RTchange(alltbt,trialTypes,metadata,fakeCueInd)

% separate cue-dependent and cue-independent components of RT change

%% plotting settings
jitterStep=0.01;
bins=200;

%% choose trial types to analyze
% template sequence 1 
% templateSequence1{1}=trialTypes.consumed_pellet==1 & trialTypes.led==0;
% templateSequence1{2}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
% templateSequence1{2}=trialTypes.chewing_at_trial_start==0 & trialTypes.paw_during_wheel==0;
templateSequence1{1}=trialTypes.touch_in_cued_window==1 & trialTypes.led==0 & any(alltbt.isChewing(:,200:486),2)==0;
templateSequence1{2}=trialTypes.chewing_at_trial_start==0;
nNext1=1;

% template sequence 2 
% templateSequence2{1}=trialTypes.consumed_pellet==1 & trialTypes.led==1;
% templateSequence2{2}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
% templateSequence2{2}=trialTypes.chewing_at_trial_start==0 & trialTypes.paw_during_wheel==0;
templateSequence2{1}=trialTypes.touch_in_cued_window==1 & trialTypes.led==1 & any(alltbt.isChewing(:,200:486),2)==0;
templateSequence2{2}=trialTypes.chewing_at_trial_start==0;
nNext2=1;

%% Get RT change after real cue

% Format 3 (10 args); tbt, trialTypes, metadata, templateSequence1,
%   nNext1, templateSequence2, nNext2, useAsCue, whichReach, doPlot
disp('Plotting results for real post-cue reaction times');
[~,~,rt_pairs]=compareTrialCombos(alltbt,trialTypes,metadata,templateSequence1,nNext1,templateSequence2,nNext2,'cueZone_onVoff','all_reachBatch',1);
% rt_pairs.rt_pairs1_contingent contains reaction time pairs sub-sampled to
% include only real pairs (not separated by a file split) and belonging to
% templateSequence1
real_rt_pairs1=rt_pairs.rt_pairs1_contingent;
real_rt_pairs2=rt_pairs.rt_pairs2_contingent;
pause;

%% Get RT change after fake cue

% Make fake cue for measuring cue-independent shift in reaction time
figure(); plot(nanmean(alltbt.all_reachBatch,1),'Color','k');
if isempty(fakeCueInd)
    fakeCueInd=input('Enter index of first pre-cue reaching across all trials. Will put fake cue here.');
end
hold all;
line([fakeCueInd fakeCueInd],[0 nanmax(nanmean(alltbt.all_reachBatch,1))]);

alltbt.fakeCue_for_preCueRT=zeros(size(alltbt.cue));
alltbt.fakeCue_for_preCueRT(:,fakeCueInd)=1;

disp('Plotting results for fake pre-cue reaction times');
[~,~,rt_pairs]=compareTrialCombos(alltbt,trialTypes,metadata,templateSequence1,nNext1,templateSequence2,nNext2,'fakeCue_for_preCueRT','all_reachBatch',1);
fake_rt_pairs1=rt_pairs.rt_pairs1_contingent;
fake_rt_pairs2=rt_pairs.rt_pairs2_contingent;

%% Calculate rotated matrix that orthogonalizes cue-dependent and cue-independent components
U_inv=[-0.4308   -1.2762; 0.9024   -0.6093];
mu=[-0.0631; -0.2280];
Z_input=[real_rt_pairs1 real_rt_pairs2; fake_rt_pairs1 fake_rt_pairs2];

Zpca=U_inv*(Z_input-repmat(mu,1,size(Z_input,2)));

% Plot RT change before rotation
jitter1=rand(size(real_rt_pairs1)).*jitterStep;
jitter2=rand(size(real_rt_pairs2)).*jitterStep;
plotScatter(real_rt_pairs1,fake_rt_pairs1,real_rt_pairs2,fake_rt_pairs2,jitter1,jitter2,'RT change before rotation');

% Plot RT change after rotation
plotScatter(Zpca(1,1:length(real_rt_pairs1)),Zpca(2,1:length(real_rt_pairs1)),Zpca(1,length(real_rt_pairs1)+1:end),Zpca(2,length(real_rt_pairs1)+1:end),jitter1,jitter2,'RT change after rotation');

% Plot summaries of RT changes in orthogonalized dimensions
% Note that for dim 1, shift to the LEFT in CDF is RT speeding up
% Note that for dim 2, shift to the right in CDF is RT speeding up
% 
% Flipping sign in histograms for dim 1
% For consistency, so in histograms
% Note that for dim 1, shift to the right in CDF is RT speeding up
% Note that for dim 2, shift to the right in CDF is RT speeding up
dim=1;
plotHist(-Zpca(dim,1:length(real_rt_pairs1)),-Zpca(dim,length(real_rt_pairs1)+1:end),bins,'Histo of dim 1','x-vals');
plotCDF(-Zpca(dim,1:length(real_rt_pairs1)),-Zpca(dim,length(real_rt_pairs1)+1:end),bins,'CDF of dim 1');
testRanksum(-Zpca(dim,1:length(real_rt_pairs1)),-Zpca(dim,length(real_rt_pairs1)+1:end),1);
dim=2;
plotHist(Zpca(dim,1:length(real_rt_pairs1)),Zpca(dim,length(real_rt_pairs1)+1:end),bins,'Histo of dim 2','y-vals');
plotCDF(Zpca(dim,1:length(real_rt_pairs1)),Zpca(dim,length(real_rt_pairs1)+1:end),bins,'CDF of dim 2');
testRanksum(Zpca(dim,1:length(real_rt_pairs1)),Zpca(dim,length(real_rt_pairs1)+1:end),1);

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
    
s=scatter(useX,useY,[],'k','filled');
pause;
m=get(s,'MarkerHandle');
alpha=0.3;
m.FaceColorData=uint8(255*[0;0;0;alpha]);
xlabel('real cue');
ylabel('fake cue');
title(tit);

hold on;
if length(jitter2)==1
    useX=RT_pairs2_x+rand(size(RT_pairs2_x)).*jitter;
    useY=RT_pairs2_y+rand(size(RT_pairs2_y)).*jitter;
else
    useX=RT_pairs2_x+jitter2;
    useY=RT_pairs2_y+jitter2;
end
s=scatter(useX,useY,[],'r','filled');
pause;
m=get(s,'MarkerHandle');
alpha=0.3;
m.FaceColorData=uint8(255*[1;0;0;alpha]);

end

function p=testRanksum(data1,data2,dispStuff)

if dispStuff==1
    disp('median of data1');
    disp(nanmedian(data1,2));
    disp('median of data2');
    disp(nanmedian(data2,2));
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





