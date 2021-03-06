function comparePreCueAndCuedRTs(alltbt,trialTypes,metadata,fakeCueInd1,fakeCueInd2)

jitter=0.01;
bins=200;

% template sequence 1 
templateSequence1{1}=trialTypes.touched_pellet==1 & trialTypes.cued_reach==1 & trialTypes.led==0;
templateSequence1{2}=trialTypes.chewing_at_trial_start==0 & trialTypes.paw_during_wheel==0;
% templateSequence1{1}=trialTypes.touch_in_cued_window==1 & trialTypes.led==0;
% templateSequence1{2}=trialTypes.touch_in_cued_window==1 & trialTypes.led==0;
% templateSequence1{3}=trialTypes.touch_in_cued_window==1 & trialTypes.led==0;
% templateSequence1{4}=trialTypes.chewing_at_trial_start==0 & trialTypes.paw_during_wheel==0;
nNext1=1;

% template sequence 2 
templateSequence2{1}=trialTypes.touched_pellet==1 & trialTypes.cued_reach==1 & trialTypes.led==1;
templateSequence2{2}=trialTypes.chewing_at_trial_start==0 & trialTypes.paw_during_wheel==0;
% templateSequence2{1}=trialTypes.touched_pellet==1 & trialTypes.cued_reach==1 & trialTypes.led==0;
% templateSequence2{2}=trialTypes.touched_pellet==1 & trialTypes.cued_reach==1 & trialTypes.led==0;
% templateSequence2{3}=trialTypes.touched_pellet==1 & trialTypes.cued_reach==1 & trialTypes.led==0;
% templateSequence2{4}=trialTypes.chewing_at_trial_start==0 & trialTypes.paw_during_wheel==0;
nNext2=1;

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


% Make fake cue for measuring cue-independent shift in reaction time
figure(); plot(nanmean(alltbt.all_reachBatch,1),'Color','k');
if isempty(fakeCueInd1)
    fakeCueInd1=input('Enter index of first pre-cue reaching across all trials. Will put fake cue here.');
end
hold all;
line([fakeCueInd1 fakeCueInd1],[0 nanmax(nanmean(alltbt.all_reachBatch,1))]);

alltbt.fakeCue_for_preCueRT=zeros(size(alltbt.cue));
alltbt.fakeCue_for_preCueRT(:,fakeCueInd1)=1;

disp('Plotting results for fake pre-cue reaction times');
[~,~,rt_pairs]=compareTrialCombos(alltbt,trialTypes,metadata,templateSequence1,nNext1,templateSequence2,nNext2,'fakeCue_for_preCueRT','all_reachBatch',1);
fake_rt_pairs1=rt_pairs.rt_pairs1_contingent;
fake_rt_pairs2=rt_pairs.rt_pairs2_contingent;


% Make fake cue 2 for measuring cue-independent shift in reaction time
% figure(); plot(nanmean(alltbt.all_reachBatch,1),'Color','k');
% if isempty(fakeCueInd2)
%     fakeCueInd2=input('Enter index of later pre-cue reaching across all trials. Will put fake cue here.');
% end
% hold all;
% line([fakeCueInd2 fakeCueInd2],[0 nanmax(nanmean(alltbt.all_reachBatch,1))]);
% 
% alltbt.fakeCue_for_preCueRT=zeros(size(alltbt.cue));
% alltbt.fakeCue_for_preCueRT(:,fakeCueInd2)=1;
% 
% disp('Plotting results for fake pre-cue reaction times');
% [~,~,rt_pairs_postcue]=compareTrialCombos(alltbt,trialTypes,metadata,templateSequence1,nNext1,templateSequence2,nNext2,'fakeCue_for_preCueRT','all_reachBatch',1);
% fake_rt_pairs1_postcue=rt_pairs_postcue.rt_pairs1_contingent;
% fake_rt_pairs2_postcue=rt_pairs_postcue.rt_pairs2_contingent;
% 
% figure(); 
% scatter3(real_rt_pairs1,fake_rt_pairs1,fake_rt_pairs1_postcue);
% % fake_rt_pairs1=fake_rt_pairs1_postcue;
% % fake_rt_pairs2=fake_rt_pairs2_postcue;

% 3D PCA
plotScatter(real_rt_pairs1,fake_rt_pairs1,real_rt_pairs2,fake_rt_pairs2,jitter,'RT change dependent on cue');
fake_rt_pairs1_postcue=fake_rt_pairs1;
fake_rt_pairs2_postcue=fake_rt_pairs2;
n_PCs=3;
temp1=[real_rt_pairs1 real_rt_pairs2];
temp2=[fake_rt_pairs1 fake_rt_pairs2];
temp3=[fake_rt_pairs1_postcue fake_rt_pairs2_postcue];
% [ica_out1,W,T,mu]=fastICA([temp1; temp2; temp3],n_PCs);
temp_together=[temp1; temp2; temp3];
[Zpca,U,mu,eigVecs] = PCA([temp1; temp2; temp3],n_PCs);
ica_out1=Zpca;
% plotScatter(ica_out1(1,1:length(real_rt_pairs1)),ica_out1(3,1:length(real_rt_pairs1)),ica_out1(1,length(real_rt_pairs1)+1:end),ica_out1(3,length(real_rt_pairs1)+1:end),jitter,'RT change dependent on cue');
plotScatter(ica_out1(1,1:length(real_rt_pairs1)),ica_out1(2,1:length(real_rt_pairs1)),ica_out1(1,length(real_rt_pairs1)+1:end),ica_out1(2,length(real_rt_pairs1)+1:end),jitter,'RT change dependent on cue');
title('Best dims from PCA');
% plotScatter(ica_out1(3,1:length(real_rt_pairs1)),ica_out1(2,1:length(real_rt_pairs1)),ica_out1(3,length(real_rt_pairs1)+1:end),ica_out1(2,length(real_rt_pairs1)+1:end),jitter,'RT change dependent on cue');
figure();
scatter3(ica_out1(1,1:length(real_rt_pairs1)),ica_out1(2,1:length(real_rt_pairs1)),ica_out1(3,1:length(real_rt_pairs1)),'MarkerFaceColor','k');
hold on;
scatter3(ica_out1(1,length(real_rt_pairs1)+1:end),ica_out1(2,length(real_rt_pairs1)+1:end),ica_out1(3,length(real_rt_pairs1)+1:end),'MarkerFaceColor','r');

figure(); 
title('Input samples and scaled eigenvectors from PCA');
scatter3(real_rt_pairs1,fake_rt_pairs1,fake_rt_pairs1_postcue,'MarkerFaceColor','k');
hold on;
scatter3(real_rt_pairs2,fake_rt_pairs2,fake_rt_pairs2_postcue,'MarkerFaceColor','r');
x=@(i,j) temp_together(i)+[0 eigVecs(i,j)];
plot3(x(1,1),x(2,1),x(3,1),'b','Linewidth',3);
plot3(x(1,2),x(2,2),x(3,2),'g','Linewidth',3);
plot3(x(1,3),x(2,3),x(3,3),'m','Linewidth',3);
% for j=1:n_PCs % this many PCs
%     plot3(x(1,j),x(2,j),x(3,j),'b','Linewidth',3);
% end

figure();
scatter(real_rt_pairs1,fake_rt_pairs1,'MarkerFaceColor','k');
hold on; 
scatter(real_rt_pairs2,fake_rt_pairs2,'MarkerFaceColor','r');
plot(x(1,1),x(2,1),'b','Linewidth',3);
plot(x(1,2),x(2,2),'g','Linewidth',3);
figure();
Zr = U * Zpca + repmat(mu,1,length(temp1));
scatter(Zr(1,1:length(real_rt_pairs1)),Zr(2,1:length(real_rt_pairs1)),'MarkerFaceColor','k');
hold on;
scatter(Zr(1,length(real_rt_pairs1)+1:end),Zr(2,length(real_rt_pairs1)+1:end),'MarkerFaceColor','r');
title('Zr');

dim=1;
plotHist(ica_out1(dim,1:length(real_rt_pairs1)),ica_out1(dim,length(real_rt_pairs1)+1:end),bins,'Histo of ICA dim 1','y-vals');
plotCDF(ica_out1(dim,1:length(real_rt_pairs1)),ica_out1(dim,length(real_rt_pairs1)+1:end),bins,'CDF of ICA dim 1');
dim=2;
plotHist(ica_out1(dim,1:length(real_rt_pairs1)),ica_out1(dim,length(real_rt_pairs1)+1:end),bins,'Histo of ICA dim 2','y-vals');
plotCDF(ica_out1(dim,1:length(real_rt_pairs1)),ica_out1(dim,length(real_rt_pairs1)+1:end),bins,'CDF of ICA dim 2');
dim=3;
plotHist(ica_out1(dim,1:length(real_rt_pairs1)),ica_out1(dim,length(real_rt_pairs1)+1:end),bins,'Histo of ICA dim 3','y-vals');
plotCDF(ica_out1(dim,1:length(real_rt_pairs1)),ica_out1(dim,length(real_rt_pairs1)+1:end),bins,'CDF of ICA dim 3');
testRanksum(ica_out1(1,1:length(real_rt_pairs1)),ica_out1(1,length(real_rt_pairs1)+1:end),1);
% testRanksum(ica_out1(2,1:length(real_rt_pairs1)),ica_out1(2,length(real_rt_pairs1)+1:end),1);
% testRanksum(ica_out1(3,1:length(real_rt_pairs1)),ica_out1(3,length(real_rt_pairs1)+1:end),1);
return






% n_PCs=2;
% temp1=[real_rt_pairs1 real_rt_pairs2];
% temp2=[fake_rt_pairs1 fake_rt_pairs2];
% temp_together=[temp1; temp2];
% [Zpca,~,~,eigVecs] = PCA([temp1; temp2],n_PCs);
% ica_out1=Zpca;
% plotScatter(ica_out1(1,1:length(real_rt_pairs1)),ica_out1(2,1:length(real_rt_pairs1)),ica_out1(1,length(real_rt_pairs1)+1:end),ica_out1(2,length(real_rt_pairs1)+1:end),jitter,'RT change dependent on cue');
% 
% figure(); 
% title('Input samples and scaled eigenvectors from PCA');
% scatter(real_rt_pairs1,fake_rt_pairs1,'MarkerFaceColor','k');
% hold on;
% scatter(real_rt_pairs2,fake_rt_pairs2,'MarkerFaceColor','r');
% x=@(i,j) temp_together(i)+[0 eigVecs(i,j)];
% plot(x(1,1),x(2,1),'b','Linewidth',3);
% plot(x(1,2),x(2,2),'g','Linewidth',3);
% 
% dim=1;
% plotHist(ica_out1(dim,1:length(real_rt_pairs1)),ica_out1(dim,length(real_rt_pairs1)+1:end),bins,'Histo of ICA dim 1','y-vals');
% plotCDF(ica_out1(dim,1:length(real_rt_pairs1)),ica_out1(dim,length(real_rt_pairs1)+1:end),bins,'CDF of ICA dim 1');
% dim=2;
% plotHist(ica_out1(dim,1:length(real_rt_pairs1)),ica_out1(dim,length(real_rt_pairs1)+1:end),bins,'Histo of ICA dim 2','y-vals');
% plotCDF(ica_out1(dim,1:length(real_rt_pairs1)),ica_out1(dim,length(real_rt_pairs1)+1:end),bins,'CDF of ICA dim 2');
% testRanksum(ica_out1(1,1:length(real_rt_pairs1)),ica_out1(1,length(real_rt_pairs1)+1:end),1);
% testRanksum(ica_out1(2,1:length(real_rt_pairs1)),ica_out1(2,length(real_rt_pairs1)+1:end),1);
% 
% return












% % Plot real vs fake RT pairs
plotScatter(real_rt_pairs1,fake_rt_pairs1,real_rt_pairs2,fake_rt_pairs2,jitter,'RT change dependent on cue');
% plotScatter(real_rt_pairs1,fake_rt_pairs1_postcue,real_rt_pairs2,fake_rt_pairs2_postcue,jitter,'RT change dependent on cue');
x=real_rt_pairs1;
y=fake_rt_pairs1;
% Three modes in this graph
% mode 1: y=x
% mode 2: x=0
% mode 3: y=0
% Remove y=x

% Remove trials that do not vary with real cue (lobes above and below x
% axis)
x2_1=real_rt_pairs1;
y2_1=fake_rt_pairs1;
x2_2=real_rt_pairs2;
y2_2=fake_rt_pairs2;
use1=zeros(1,length(x2_1));
use2=zeros(1,length(x2_2));
for i=1:length(x2_1)
%     if ~(abs(y2_1(i))>0.2 & abs(x2_1(i))<0.1)
    if ~(abs(x2_1(i))<0.1)    
        use1(i)=1;
    end
end
for i=1:length(x2_2)
%     if ~(abs(y2_2(i))>0.2 & abs(x2_2(i))<0.1)
    if ~(abs(x2_2(i))<0.1)
        use2(i)=1;
    end
end        
x2_1(use1==0)=nan;
y2_1(use1==0)=nan;
x2_2(use2==0)=nan;
y2_2(use2==0)=nan;
plotScatter(x2_1,y2_1,x2_2,y2_2,jitter,'Removed lobes above and below origin');
real_rt_pairs1=x2_1;
fake_rt_pairs1=y2_1;
real_rt_pairs2=x2_2;
fake_rt_pairs2=y2_2;

% Remove trials along y=x
% use1=zeros(1,length(x2_1));
% use2=zeros(1,length(x2_2));
% for i=1:length(x2_1)
%     if ~(y2_1(i)<=x2_1(i)+0.06 & y2_1(i)>=x2_1(i)-0.06)
%         use1(i)=1;
%     end
% end
% for i=1:length(x2_2)
%     if ~(y2_2(i)<=x2_2(i)+0.06 & y2_2(i)>=x2_2(i)-0.06)
%         use2(i)=1;
%     end
% end
% x2_1(use1==0)=nan;
% y2_1(use1==0)=nan;
% x2_2(use2==0)=nan;
% y2_2(use2==0)=nan;
% plotScatter(x2_1,y2_1,x2_2,y2_2,jitter,'Removed y=x');



% Plot rotated scatter to isolate cue-dependent and cue-independent
% components
% If x and y are same, make this y=0 (y is cue-dependent component)
angle=-(1.5/4)*pi;
x2_1=real_rt_pairs1.*cos(angle)-fake_rt_pairs1.*sin(angle);
y2_1=fake_rt_pairs1.*cos(angle)+real_rt_pairs1.*sin(angle);
x2_2=real_rt_pairs2.*cos(angle)-fake_rt_pairs2.*sin(angle);
y2_2=fake_rt_pairs2.*cos(angle)+real_rt_pairs2.*sin(angle);
plotScatter(x2_1,y2_1,x2_2,y2_2,jitter,'RT change dependent on cue -- rotated');
% % Find trials along y axis
% x_cutoff=0.5;
% use1=abs(x2_1)<=x_cutoff;
% use2=abs(x2_2)<=x_cutoff;
% x2_1(use1==0)=nan;
% y2_1(use1==0)=nan;
% x2_2(use2==0)=nan;
% y2_2(use2==0)=nan;
% plotScatter(x2_1,y2_1,x2_2,y2_2,jitter,'Removed x not along y axis');
% Find trials where "real cue" effect is greater than "fake cue" effect
% use1=zeros(1,length(x2_1));
% for i=1:length(x2_1)
%     if abs(y2_1(i))>abs(x2_1(i)) % real cue effect stronger
%         use1(i)=1;
%     end
% end
% use2=zeros(1,length(x2_2));
% for i=1:length(x2_2)
%     if abs(y2_2(i))>abs(x2_2(i)) % real cue effect stronger
%         use2(i)=1;
%     end
% end
% x2_1(use1==0)=nan;
% y2_1(use1==0)=nan;
% x2_2(use2==0)=nan;
% y2_2(use2==0)=nan;
% plotScatter(x2_1,y2_1,x2_2,y2_2,jitter,'Removed x>y');
plotHist(y2_1,y2_2,bins,'Histo of y','y-vals');
plotCDF(y2_1,y2_2,bins,'CDF of y');
testRanksum(y2_1,y2_2,1);

% angle=-(1/4)*pi;
% x2_1=real_rt_pairs1.*cos(angle)-fake_rt_pairs1.*sin(angle);
% y2_1=fake_rt_pairs1.*cos(angle)+real_rt_pairs1.*sin(angle);
% x2_2=real_rt_pairs2.*cos(angle)-fake_rt_pairs2.*sin(angle);
% y2_2=fake_rt_pairs2.*cos(angle)+real_rt_pairs2.*sin(angle);
% plotScatter(x2_1,y2_1,x2_2,y2_2,jitter,'RT change dependent on cue -- rotated');
% y_dev=0.16;
% y2_1_backup=y2_1;
% y2_2_backup=y2_2;
% x2_1(abs(y2_1_backup)<y_dev)=nan;
% y2_1(abs(y2_1_backup)<y_dev)=nan;
% x2_2(abs(y2_2_backup)<y_dev)=nan;
% y2_2(abs(y2_2_backup)<y_dev)=nan;
% plotScatter(x2_1,y2_1,x2_2,y2_2,jitter,'RT change dependent on cue -- removed cue-independent component');
% 
% angle=(1/4)*pi;
% x2_1=x2_1.*cos(angle)-y2_1.*sin(angle);
% y2_1=y2_1.*cos(angle)+x2_1.*sin(angle);
% x2_2=x2_2.*cos(angle)-y2_2.*sin(angle);
% y2_2=y2_2.*cos(angle)+x2_2.*sin(angle);
% plotScatter(x2_1,y2_1,x2_2,y2_2,jitter,'Rotated 2nd time');
% x_dev=0.06;
% x2_1_backup=x2_1;
% x2_2_backup=x2_2;
% x2_1(abs(x2_1_backup)<x_dev)=nan;
% y2_1(abs(x2_1_backup)<x_dev)=nan;
% x2_2(abs(x2_2_backup)<x_dev)=nan;
% y2_2(abs(x2_2_backup)<x_dev)=nan;
% plotScatter(x2_1,y2_1,x2_2,y2_2,jitter,'Rotated 2nd time -- then removed unchanging component');
% plotHist(x2_1,x2_2,bins,'Histo of x','x-vals');
% plotCDF(x2_1,x2_2,bins,'CDF of x');
% testRanksum(x2_1,x2_2,1);

% Remove component of x variance that can be explained by y=x
% [x_1,y_1]=regressOut(x,y,y,y);
% x=real_rt_pairs2;
% y=fake_rt_pairs2;
% [x_2,y_2]=regressOut(x,y,y,y);
% plotScatter(x_1,y_1,x_2,y_2,jitter,'Regressed out y=x');
% plotHist(x_1,x_2,bins,'Histo of x','x-vals');
% plotCDF(x_1,x_2,bins,'CDF of x');
% testRanksum(x_1,x_2,1);

% [x,y]=regressOut(x,y,explainedX)
% regressOut(x,y,explainedX);
% 
% 
% 
% Plot rotated scatter to isolate cue-dependent and cue-independent
% components
% % If x and y are same, make this y=0 (y is cue-dependent component)
% angle=-(1.5/4)*pi;
% x2_1=real_rt_pairs1.*cos(angle)-fake_rt_pairs1.*sin(angle);
% y2_1=fake_rt_pairs1.*cos(angle)+real_rt_pairs1.*sin(angle);
% x2_2=real_rt_pairs2.*cos(angle)-fake_rt_pairs2.*sin(angle);
% y2_2=fake_rt_pairs2.*cos(angle)+real_rt_pairs2.*sin(angle);
% plotScatter(x2_1,y2_1,x2_2,y2_2,jitter,'RT change dependent on cue -- rotated');
% plotHist(y2_1,y2_2,bins,'Histo of y from rotated plot','y-vals');
% plotCDF(y2_1,y2_2,bins,'CDF of y from rotated plot');
% testRanksum(y2_1,y2_2,1);


% n_PCs=3;
% ica(x2_1,y2_1,x2_2,y2_2,n_PCs);

end

function [x,y]=regressOut(x,y,explainedX,explainedY)

% Remove component of x that can be explained by explainedX
x=x-explainedX;
y=y-explainedY;

end

function ica(x2_1,y2_1,x2_2,y2_2,n_PCs)

% Inputs:       Z is a d x n matrix containing n samples of d-dimensional
% %               data
% x2_1=turnToRowVector(x2_1);
% y2_1=turnToRowVector(y2_1);
% x2_2=turnToRowVector(x2_2);
% y2_2=turnToRowVector(y2_2);
% 
% % make A
% A=zeros(length(x2_1),length(y2_1));
% for i=1:length(y2_1)
%     % y2_1_ith=diag_val*x2_1_ith
%     diag_val=y2_1(i)/x2_1(i);
%     A(i,i)=diag_val;
% end
% A(isnan(A))=nanmax(A(1:end));
% [eig_vecs,eig_val_D]=eig(A);
% eig_vals=nan(1,length(y2_1));
% for i=1:length(y2_1)
%     eig_vals(i)=eig_val_D(i,i);
% end
%     
%     
% 
% disp('hi');

% [ica_out1,W,T,mu]=fastICA([x2_1; y2_1],n_PCs);
% Zpca = PCA(A,n_PCs);
% ica_out1=Zpca;
% % 
% % % Plot output
% if n_PCs==3
%     figure();
%     scatter3(ica_out1(1,:),ica_out1(2,:),ica_out1(3,:));
% else
%     disp('too many PCs to plot as 3D scatter');
% end

% for i = 1:r
%     plot(Zfica(i,:),'-','Color',cm(i,:)); hold on;
% end
% title('Independent components [Fast ICA]');
% axis tight;

end

function x2_1=turnToRowVector(x2_1)

if size(x2_1,1)>1
    % make row vector
    x2_1=x2_1';
end
if size(x2_1,1)>1
    error('improper data structure size in ica');
end
    
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

function plotScatter(RT_pairs1_x,RT_pairs1_y,RT_pairs2_x,RT_pairs2_y,jitter,tit)

figure();
useX=RT_pairs1_x+rand(size(RT_pairs1_x)).*jitter;
useY=RT_pairs1_y+rand(size(RT_pairs1_y)).*jitter;
s=scatter(useX,useY,[],'k','filled');
pause;
m=get(s,'MarkerHandle');
alpha=0.3;
m.FaceColorData=uint8(255*[0;0;0;alpha]);
xlabel('real cue');
ylabel('fake cue');
title(tit);

hold on;
useX=RT_pairs2_x+rand(size(RT_pairs2_x)).*jitter;
useY=RT_pairs2_y+rand(size(RT_pairs2_y)).*jitter;
s=scatter(useX,useY,[],'r','filled');
pause;
m=get(s,'MarkerHandle');
alpha=0.3;
m.FaceColorData=uint8(255*[1;0;0;alpha]);

end



