function trialTypeDecode(tensor,allLabels,timepoints_for_tensor)

disp(['Rows of tensor must be ']);
disp('Cued grp 1 cells');
disp('Cued grp 2 cells');
disp('Uncued grp 1 cells');
disp('Uncued grp 2 cells');

timebin1=[0 8*0.06*5]; % in secs
timebin2=[8*0.06*5 18*0.06*5]; % in secs
timebin3=[5*0.06*5 18*0.06*5]; % in secs
overweightCueNeurons=1;

% Fill nans w zeros
tensor(isnan(tensor))=0;

% Get inds in tensor corresponding to these time bins
indbin1=find(timepoints_for_tensor>=timebin1(1),1,'first'):find(timepoints_for_tensor<=timebin1(2),1,'last');
indbin2=find(timepoints_for_tensor>=timebin2(1),1,'first'):find(timepoints_for_tensor<=timebin2(2),1,'last');
indbin3=find(timepoints_for_tensor>=timebin3(1),1,'first'):find(timepoints_for_tensor<=timebin3(2),1,'last');

% Get cue axis values
cue_axis_vals_part1=overweightCueNeurons*mean(mean(tensor(1:2,indbin1,:),1,'omitnan'),2,'omitnan')-mean(mean(tensor(3:4,indbin1,:),1,'omitnan'),2,'omitnan');
cue_axis_vals_part2=overweightCueNeurons*mean(mean(tensor(1:2,indbin2,:),1,'omitnan'),2,'omitnan')-mean(mean(tensor(3:4,indbin2,:),1,'omitnan'),2,'omitnan');
cue_axis_vals_part2=-cue_axis_vals_part2;
cue_axis_vals=cue_axis_vals_part1-cue_axis_vals_part2;

% Get outcome axis values
out_axis_vals=mean(mean(tensor([1 3],indbin3,:),1,'omitnan'),2,'omitnan')-mean(mean(tensor([2 4],indbin3,:),1,'omitnan'),2,'omitnan');

% % % Drop zeros trials, assumption is that when all cells don't spike, we have no info
todrop=abs(out_axis_vals)<0.1 | abs(cue_axis_vals)<0.1;
% % % Remove outliers
% % [~,rm1]=rmoutliers(squeeze(cue_axis_vals),"median","ThresholdFactor",4);
% % [~,rm2]=rmoutliers(squeeze(out_axis_vals),"median","ThresholdFactor",4);
% % todrop(rm1)=1; todrop(rm2)=1;
% % % Drop
cue_axis_vals=cue_axis_vals(:,:,todrop==0);
out_axis_vals=out_axis_vals(:,:,todrop==0);
allLabels=allLabels(todrop==0);

% Demix
[U,S,V]=svd([squeeze(out_axis_vals) squeeze(cue_axis_vals)]);
out_axis_vals(1,1,:)=U(:,1)*S(1,1);
cue_axis_vals(1,1,:)=U(:,2)*S(2,2);
% I'm assuming Matlab always produces positive singular values
if S(1,1)<0 || S(2,2)<0
    error(['Negative singular values']);
end

% Reorient
% According to my hypothesis, grp1 > grp2 defines outcome axis
% and cued > uncued defines cue axis
% To orient in consistent way across sessions, project example vecs
% This does nothing to the data, just changes axis labels and produces
% consistent plot orientation
example_vec_outPos=[1 0]'; % x axis of input was outcome, y axis of input was cue
example_vec_cuePos=[0 1]';
projXaxis_ontoS1=example_vec_outPos.*V(:,1);
projXaxis_ontoS2=example_vec_outPos.*V(:,2);
projYaxis_ontoS1=example_vec_cuePos.*V(:,1);
projYaxis_ontoS2=example_vec_cuePos.*V(:,2);
if abs(V(2,1))<=abs(V(1,1)) % SV1 is X axis
    % good, leave alone
    if projXaxis_ontoS1(1)>0 % good, leave alone
    else
        out_axis_vals=-out_axis_vals; % flip
    end
    if projYaxis_ontoS2(2)>0 % good, leave alone
    else
        cue_axis_vals=-cue_axis_vals; % flip
    end
else % SV1 is Y axis
    % exchange X and Y for plot
    temp=cue_axis_vals;
    cue_axis_vals=out_axis_vals;
    out_axis_vals=temp;
    if projXaxis_ontoS2(1)>0 % good, leave alone
    else
        cue_axis_vals=-cue_axis_vals;
    end
    if projYaxis_ontoS1(2)>0 % good, leave alone
    else
        out_axis_vals=-out_axis_vals; % flip
    end
end

% Plot
c{1}='b'; c{2}='g'; c{3}='r'; c{4}='k';
figure();
unique_allLabels=unique(allLabels);
meanOfAll=[nanmean(squeeze(out_axis_vals(:,:,:))) nanmean(squeeze(cue_axis_vals(:,:,:)))];
for i=1:length(unique_allLabels)
    r=rand(size(squeeze(out_axis_vals(:,:,allLabels==unique_allLabels(i)))))*0.2;
    scatter(squeeze(out_axis_vals(:,:,allLabels==unique_allLabels(i)))+r,squeeze(cue_axis_vals(:,:,allLabels==unique_allLabels(i)))+r,[],c{i}); hold on;
    scatter(nanmean(squeeze(out_axis_vals(:,:,allLabels==unique_allLabels(i)))),nanmean(squeeze(cue_axis_vals(:,:,allLabels==unique_allLabels(i)))),[],c{i},'filled');
    line([meanOfAll(1) nanmean(squeeze(out_axis_vals(:,:,allLabels==unique_allLabels(i))))],[meanOfAll(2) nanmean(squeeze(cue_axis_vals(:,:,allLabels==unique_allLabels(i))))],'LineWidth',2,'Color',c{i});
end
xlabel('Outcome axis'); ylabel('Cue axis');
legend({'cued succ','','','uncued succ','','','cued fail','','','uncued fail','',''});

% Plot trial label shuffle
c{1}='b'; c{2}='g'; c{3}='r'; c{4}='k';
figure();
unique_allLabels=unique(allLabels);
meanOfAll=[nanmean(squeeze(out_axis_vals(:,:,:))) nanmean(squeeze(cue_axis_vals(:,:,:)))];
backup_allLabels=allLabels;
allLabels=allLabels(randperm(length(allLabels)));
for i=1:length(unique_allLabels)
    r=rand(size(squeeze(out_axis_vals(:,:,allLabels==unique_allLabels(i)))))*0.2;
    scatter(squeeze(out_axis_vals(:,:,allLabels==unique_allLabels(i)))+r,squeeze(cue_axis_vals(:,:,allLabels==unique_allLabels(i)))+r,[],c{i}); hold on;
    scatter(nanmean(squeeze(out_axis_vals(:,:,allLabels==unique_allLabels(i)))),nanmean(squeeze(cue_axis_vals(:,:,allLabels==unique_allLabels(i)))),[],c{i},'filled');
    line([meanOfAll(1) nanmean(squeeze(out_axis_vals(:,:,allLabels==unique_allLabels(i))))],[meanOfAll(2) nanmean(squeeze(cue_axis_vals(:,:,allLabels==unique_allLabels(i))))],'LineWidth',2,'Color',c{i});
end
title('TRIAL LABEL SHUFFLE');
xlabel('Outcome axis'); ylabel('Cue axis');
legend({'cued succ','','','uncued succ','','','cued fail','','','uncued fail','',''});

end