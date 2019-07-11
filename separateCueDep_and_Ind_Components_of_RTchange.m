function separateCueDep_and_Ind_Components_of_RTchange(alltbt,trialTypes,metadata,fakeCueInd)

% separate cue-dependent and cue-independent components of RT change

%% plotting settings
jitterStep=0.01;
bins=200;

%% choose trial types to analyze
% template sequence 1 
% templateSequence1{1}=trialTypes.success_in_cued_window==1 & trialTypes.led==0 & any(alltbt.isChewing(:,200:486),2)==0;
% templateSequence1{2}=trialTypes.chewing_at_trial_start==0 & trialTypes.paw_during_wheel==0 & metadata.mouseid==1;
% templateSequence1{2}=trialTypes.chewing_at_trial_start==0 & trialTypes.paw_during_wheel==0;
% templateSequence1{1}=trialTypes.touch_in_cued_window==1 & trialTypes.led==0 & any(alltbt.isChewing(:,200:486),2)==0;
% templateSequence1{1}=trialTypes.touched_pellet==1 & trialTypes.cued_reach==1 & trialTypes.led==0 & any(alltbt.isChewing(:,200:486),2)==0;
templateSequence1{1}=trialTypes.touched_pellet==1 & any(alltbt.all_reachBatch(:,94:151),2)==1 & trialTypes.led==0 & metadata.dprimes>0 & metadata.dprimes<2;
templateSequence1{2}=trialTypes.touched_pellet==1 & any(alltbt.all_reachBatch(:,94:151),2)==1 & trialTypes.consumed_pellet==0 & trialTypes.led==0 & metadata.dprimes>0 & metadata.dprimes<2;
templateSequence1{3}=trialTypes.chewing_at_trial_start==0 & metadata.dprimes>0 & metadata.dprimes<2;
% templateSequence1{2}=trialTypes.chewing_at_trial_start==0;
% templateSequence1{3}=trialTypes.chewing_at_trial_start==0; any(alltbt.all_reachBatch(:,94:151),2)==1


templateSequence1{1}=trialTypes.touched_pellet==1 & trialTypes.led==0;
templateSequence1{2}=trialTypes.touched_pellet==1 & trialTypes.led==0;
% templateSequence1{3}=trialTypes.touched_pellet==1 & trialTypes.led==0;
% templateSequence1{4}=trialTypes.touched_pellet==1 & trialTypes.led==0;
% templateSequence1{5}=trialTypes.touched_pellet==1 & trialTypes.led==0;
% templateSequence1{6}=trialTypes.touched_pellet==1 & trialTypes.led==0;
% templateSequence1{7}=trialTypes.touched_pellet==1 & trialTypes.led==0;
% templateSequence1{8}=trialTypes.touched_pellet==1 & trialTypes.led==0;
templateSequence1{3}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
nNext1=2;

% template sequence 2 
% templateSequence2{1}=trialTypes.success_in_cued_window==0 & trialTypes.after_cue_drop==1 & trialTypes.led==0 & any(alltbt.isChewing(:,200:486),2)==0;
% templateSequence2{2}=trialTypes.chewing_at_trial_start==0 & trialTypes.paw_during_wheel==0 & metadata.mouseid==1;
% templateSequence2{2}=trialTypes.chewing_at_trial_start==0 & trialTypes.paw_during_wheel==0;
% templateSequence2{1}=trialTypes.touch_in_cued_window==0 & trialTypes.led==0 & any(alltbt.isChewing(:,200:486),2)==0;
% templateSequence2{1}=trialTypes.touched_pellet==1 & trialTypes.cued_reach==1 & trialTypes.led==1 & any(alltbt.isChewing(:,200:486),2)==0;
templateSequence2{1}=trialTypes.touched_pellet==0 & any(alltbt.all_reachBatch(:,94:151),2)==1 & trialTypes.led==0 & metadata.dprimes>0 & metadata.dprimes<2;
templateSequence2{2}=trialTypes.touched_pellet==0 & trialTypes.consumed_pellet==0 & trialTypes.led==0 & metadata.dprimes>0 & metadata.dprimes<2;
templateSequence2{3}=trialTypes.chewing_at_trial_start==0 & metadata.dprimes>0 & metadata.dprimes<2;
% templateSequence2{2}=trialTypes.chewing_at_trial_start==0; any(alltbt.all_reachBatch(:,94:151),2)==1 
% templateSequence2{3}=trialTypes.chewing_at_trial_start==0;
nNext2=2;

templateSequence2{1}=trialTypes.touched_pellet==0 & trialTypes.led==0;
templateSequence2{2}=trialTypes.touched_pellet==0 & trialTypes.led==0;
% templateSequence2{3}=trialTypes.touched_pellet==0 & trialTypes.led==0;
% templateSequence2{4}=trialTypes.touched_pellet==0 & trialTypes.led==0;
% templateSequence2{5}=trialTypes.touched_pellet==0 & trialTypes.led==0;
% templateSequence2{6}=trialTypes.touched_pellet==0 & trialTypes.led==0;
% templateSequence2{7}=trialTypes.touched_pellet==0 & trialTypes.led==0;
% templateSequence2{8}=trialTypes.touched_pellet==0 & trialTypes.led==0;
templateSequence2{3}=trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1;
nNext2=2;



%% Get RT change after real cue

% Format 3 (10 args); tbt, trialTypes, metadata, templateSequence1,
%   nNext1, templateSequence2, nNext2, useAsCue, whichReach, doPlot
disp('Plotting results for real post-cue reaction times');
[~,~,rt_pairs]=compareTrialCombos(alltbt,trialTypes,metadata,templateSequence1,nNext1,templateSequence2,nNext2,'cueZone_onVoff','all_reachBatch',1);
% show where these template sequences come from as a function of fraction
% through session
[n,x]=hist(metadata.nth_session(rt_pairs.sequenceMatchStarts1==1)+metadata.fractionThroughSess(rt_pairs.sequenceMatchStarts1==1),bins);
figure();
plot(x,n,'Color','k');
hold on;
[n,x]=hist(metadata.nth_session(rt_pairs.sequenceMatchStarts2==1)+metadata.fractionThroughSess(rt_pairs.sequenceMatchStarts2==1),bins);
plot(x,n,'Color','r');
title('Template match from nth session');

% Match reaction time distributions for conditions 1 and 2
[rt_pairs,resampInds]=correctInputDistribution(rt_pairs,1,[]);

% rt_pairs.rt_pairs1_contingent contains reaction time pairs sub-sampled to
% include only real pairs (not separated by a file split) and belonging to
% templateSequence1
real_rt_pairs1=rt_pairs.rt_pairs1_contingent;
real_rt_pairs2=rt_pairs.rt_pairs2_contingent;

plotHist(rt_pairs.rt_pairs1_contingent,rt_pairs.rt_pairs2_contingent,bins,'RT change','Count');
plotCDF(rt_pairs.rt_pairs1_contingent,rt_pairs.rt_pairs2_contingent,bins,'CDF');
testRanksum(rt_pairs.rt_pairs1_contingent,rt_pairs.rt_pairs2_contingent,1);

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

% Match reaction time distributions for conditions 1 and 2
[rt_pairs]=correctInputDistribution(rt_pairs,1,resampInds);

fake_rt_pairs1=rt_pairs.rt_pairs1_contingent;
fake_rt_pairs2=rt_pairs.rt_pairs2_contingent;

%% Calculate rotated matrix that orthogonalizes cue-dependent and cue-independent components
% U_inv=[-0.4308   -1.2762; 0.9024   -0.6093];
% mu=[-0.0631; -0.2280];
Z_input=[real_rt_pairs1 real_rt_pairs2; fake_rt_pairs1 fake_rt_pairs2];

% Zpca=U_inv*(Z_input-repmat(mu,1,size(Z_input,2)));

% Plot RT change before rotation
jitter1=rand(size(real_rt_pairs1)).*jitterStep;
jitter2=rand(size(real_rt_pairs2)).*jitterStep;
plotScatter(real_rt_pairs1,fake_rt_pairs1,real_rt_pairs2,fake_rt_pairs2,jitter1,jitter2,'RT change before rotation');

coeff=[0.6110    0.7885   -0.0709; ...
       0.5554   -0.4908   -0.6713; ...
       0.5641   -0.3707    0.7378];
disp(coeff);
pc1=coeff(:,1);
pc2=coeff(:,2);
pc3=coeff(:,3);

hold on; line([0 pc1(1)],[0 pc1(2)],'Color','b'); 
hold on; line([0 pc2(1)],[0 pc2(2)],'Color','g');
hold on; line([0 pc3(1)],[0 pc3(2)],'Color','c');

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

% plotScatter(Z_onto_pc1(1:length(real_rt_pairs1)),Z_onto_pc2(1:length(real_rt_pairs1)),Z_onto_pc1(length(real_rt_pairs1)+1:end),Z_onto_pc2(length(real_rt_pairs1)+1:end),jitter1,jitter2,'RT change after rotation');

% theta=-pi/4;
% R=[cos(theta) -sin(theta); sin(theta) cos(theta)];
% Zpca=[Z_onto_pc1; Z_onto_pc2];
% Zpca=R*Zpca;
% Z_onto_pc1=Zpca(1,:);
% Z_onto_pc2=Zpca(2,:);
% theta=pi/3;
% R=[cos(theta) -sin(theta); sin(theta) cos(theta)];
% Zpca=[Z_onto_pc1; Z_onto_pc2];
% Zpca=R*Zpca;
% Z_onto_pc1=Zpca(1,:);
% Z_onto_pc2=Zpca(2,:);
% theta=(10/10)*(2*pi);
% R=[cos(theta) -sin(theta); sin(theta) cos(theta)];
% Zpca=[Z_onto_pc1; Z_onto_pc2];
% Zpca=R*Zpca;
% Z_onto_pc1=Zpca(1,:);
% Z_onto_pc2=Zpca(2,:);

% Plot RT change after rotation
plotScatter(Z_onto_pc1(1:length(real_rt_pairs1)),Z_onto_pc2(1:length(real_rt_pairs1)),Z_onto_pc1(length(real_rt_pairs1)+1:end),Z_onto_pc2(length(real_rt_pairs1)+1:end),jitter1,jitter2,'RT change after rotation');

% Plot summaries of RT changes in orthogonalized dimensions
% Note that for dim 1, shift to the LEFT in CDF is RT speeding up
% Note that for dim 2, shift to the right in CDF is RT speeding up
% 
% Flipping sign in histograms for dim 1
% For consistency, so in histograms
% Note that for dim 1, shift to the right in CDF is RT speeding up
% Note that for dim 2, shift to the right in CDF is RT speeding up
% dim=1;
% plotHist(-Zpca(dim,1:length(real_rt_pairs1)),-Zpca(dim,length(real_rt_pairs1)+1:end),bins,'Histo of dim 1','x-vals');
% plotCDF(-Zpca(dim,1:length(real_rt_pairs1)),-Zpca(dim,length(real_rt_pairs1)+1:end),bins,'CDF of dim 1');
% testRanksum(-Zpca(dim,1:length(real_rt_pairs1)),-Zpca(dim,length(real_rt_pairs1)+1:end),1);
% dim=2;
% plotHist(Zpca(dim,1:length(real_rt_pairs1)),Zpca(dim,length(real_rt_pairs1)+1:end),bins,'Histo of dim 2','y-vals');
% plotCDF(Zpca(dim,1:length(real_rt_pairs1)),Zpca(dim,length(real_rt_pairs1)+1:end),bins,'CDF of dim 2');
% testRanksum(Zpca(dim,1:length(real_rt_pairs1)),Zpca(dim,length(real_rt_pairs1)+1:end),1);

Zpca=[Z_onto_pc1; Z_onto_pc2];
dim=1;
plotHist(Zpca(dim,1:length(real_rt_pairs1)),Zpca(dim,length(real_rt_pairs1)+1:end),bins,'Histo of dim 1','x-vals');
plotCDF(Zpca(dim,1:length(real_rt_pairs1)),Zpca(dim,length(real_rt_pairs1)+1:end),bins,'CDF of dim 1');
testRanksum(Zpca(dim,1:length(real_rt_pairs1)),Zpca(dim,length(real_rt_pairs1)+1:end),1);
dim=2;
plotHist(Zpca(dim,1:length(real_rt_pairs1)),Zpca(dim,length(real_rt_pairs1)+1:end),bins,'Histo of dim 2','y-vals');
plotCDF(Zpca(dim,1:length(real_rt_pairs1)),Zpca(dim,length(real_rt_pairs1)+1:end),bins,'CDF of dim 2');
testRanksum(Zpca(dim,1:length(real_rt_pairs1)),Zpca(dim,length(real_rt_pairs1)+1:end),1);

temp1=Zpca(dim,1:length(real_rt_pairs1));
temp2=Zpca(dim,length(real_rt_pairs1)+1:end);
plotCDF(temp1(temp1>-4 & temp1<4),temp2(temp2>-4 & temp2<4),bins,'CDF of dim 2 -- only changes within 2 sec faster or slower');
disp('CLOSE-UP P-VAL');
testRanksum(temp1(temp1>-4 & temp1<4),temp2(temp2>-4 & temp2<4),1);

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





